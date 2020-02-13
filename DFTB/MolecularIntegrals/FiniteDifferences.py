#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
finite differences on a spherical grid
"""
from DFTB.MolecularIntegrals.MulticenterIntegration import select_angular_grid, number_of_radial_points, atomlist2arrays, solve_radial_dgl
#from DFTB.MolecularIntegrals.MulticenterIntegration import multicenter_grids, multicenter_interpolation
from DFTB.MolecularIntegrals.LebedevQuadrature import spherical_harmonics_it, outerN
from DFTB.MolecularIntegrals.Ints1e import integral
from DFTB.MolecularIntegrals import settings
from DFTB.AtomicData import atom_names, slater_radii

import numpy as np
import numpy.linalg as la
import scipy.linalg as sla

from scipy.sparse.linalg import LinearOperator


def angular_laplacian(lebedev_order):
    """
    representation of the angular part L^2 of the Laplacian operator

       2                 1      d                      1     d^2 f
      L  f(th,ph)  = - ------- --- ( sin(th) f ) - --------- -------
                       sin(th) dth                 sin(th)^2 d ph^2

    as a matrix that operates on functions defined on a Lebdev grid:
        2
      (L . f)  = sum   lap     f
             j      j     i,j   j

    Parameters
    ----------
    lebdev_order  :   int, order of Lebedev grid

    Returns
    -------
    lap_ang       :   n x n matrix, representation of angular part of the
                      Laplacian operator on the Lebedev grid with n points
    """
    
    # load Lebedev grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    # number of angular grid points
    n = len(th)

    # matrix elements of Laplacian
    #
    #   lap
    #      i,j
    lap_ang = np.zeros((n,n), dtype=complex)

    sph_it = spherical_harmonics_it(th,ph)

    for Ylm,l,m in sph_it:
        lap_ang += l*(l+1) * 4.0*np.pi * outerN(Ylm, Ylm.conjugate()*angular_weights)

        if m == -Lmax/2:
            break

    return lap_ang

def get_radial_grid(Nr, rm):
    k = np.array(range(1,Nr+1))
    # grid points on interval [-1,1]
    zr = k/(Nr+1.0)
    xr = np.cos(zr * np.pi)
    # radial grid points on interval [0,infinity]
    r = rm * (1+xr)/(1-xr)

    return r


def radial_laplacian(Nr, rm):
    """
    matrix representation of radial part of Laplacian

      __2        1  d^2
      \/ f(r) = --- ---- (r f)
                 r  dr^2

    on a non-equidistant radial grid as described in Becke's articles. 

    Central finite difference formula for 1st and 2nd derivatives are taken from Bickley [1].
    The radial Laplacian matrix is not symmetric!

    Parameters
    ----------
    Nr            :  int, number of radial grid points
    rm            :  float, scaling factor, 1/2 * Slater radius

    Returns
    -------
    lap_rad       :  nrad x nrad matrix, representation of the
                     Laplacian on the radial grid

    References
    ----------
    [1] W. Bickley, "Formulae for Numerical Differentiation",
        The Mathematical Gazette, vol. 25, no. 263, pp. 19-27 (1941)

    """
    # number of radial grid points
    n = Nr
    # grid is equidistant in z-coordinates
    k = np.array(range(1,n+1))
    # grid points on interval [-1,1]
    zr = k/(n+1.0)
    xr = np.cos(zr * np.pi)
    # radial grid points on interval [0,infinity]
    r = rm * (1+xr)/(1-xr)

    omx = 1-xr
    opx = 1+xr

    # The transformation of partial derivatives from
    # r-coordinates to z-coordinates is accomplished as
    #  d^2 u      d u      d^2 u
    #  ----- = c1 --- + c2 -----
    #  d r^2      d z      d z^2
    c1 = 1.0/(4.0*np.pi*rm**2) * omx**(2.5) * (1+opx) / opx**(1.5)
    c2 = 1.0/(4.0*np.pi**2*rm**2) * omx**3 / opx

    # separation between equidistant points
    h = 1.0/(n+1)

    # operators d/dz and d^2/dz^2
    D1 = np.zeros((n,n))
    D2 = np.zeros((n,n))
    # non-centered seven-point formulae for i=0
    # D^1 u_0
    D1[0,0:7] = np.array([ -1764.0, +4320.0,  -5400.0,  +4800.0, -2700.0,  +864.0, -120.0  ])/(720.0*h)
    # D^2 u_0
    D2[0,0:7] = np.array([ +1624.0, -6264.0, +10530.0, -10160.0, +5940.0, -1944.0, +274.0 ])/(360.0*h**2)
    # non-centered seven-point formulae for i=1
    # D^1 u_1
    D1[1,0:7] = np.array([ -120.0, -924.0, +1800.0, -1200.0, +600.0, -180.0, +24.0 ])/(720.0*h)
    # D^2 u_1
    D2[1,0:7] = np.array([ +274.0, -294.0, -510.0, +940.0, -570.0, +186.0, -26.0 ])/(360.0*h**2)
    # non-centered seven-point formulae for i=2
    # D^1 u_2
    D1[2,0:7] = np.array([ +24.0, -288.0, -420.0, +960.0, -360.0, +96.0, -12.0])/(720.0*h)
    # D^2 u_2
    D2[2,0:7] = np.array([ -26.0, +456.0, -840.0, +400.0,  +30.0, -24.0,  +4.0])/(360.0*h**2)
    # centered seven-point formulae for i=3,...,n-4
    for i in range(3, n-3):
        D1[i,i-3:i+4] = np.array([-12.0, +108.0, -540.0,    0.0, +540.0, -108.0, +12.0])/(720.0*h)
        D2[i,i-3:i+4] = np.array([ +4.0,  -54.0, +540.0, -980.0, +540.0,  -54.0,  +4.0])/(360.0*h**2)
    # non-centered seven-point formulae for i=n-3
    # D^1 u_{n-3}   ~   D^1 u_4
    D1[n-3,n-7:] = np.array([+12.0, -96.0, +360.0, -960.0, +420.0, +288.0, -24.0 ])/(720.0*h)
    # D^2 u_{n-3}   ~   D^2 u_4
    D2[n-3,n-7:] = np.array([ +4.0, -24.0,  +30.0, +400.0, -840.0, +456.0, -26.0 ])/(360.0*h**2)
    # non-centered seven-point formulae for i=n-2
    # D^1 u_{n-2}   ~   D^1 u_5
    D1[n-2,n-7:] = np.array([ -24.0, +180.0, -600.0, +1200.0, -1800.0, +924.0, +120.0])/(720.0*h)
    # D^2 u_{n-2}   ~   D^2 u_5
    D2[n-2,n-7:] = np.array([ -26.0, +186.0, -570.0,  +940.0,  -510.0, -294.0, +274.0])/(360.0*h**2)
    # non-centered seven-point formulae for i=n-1
    # D^1 u_{n-1}   ~   D^1 u_6
    D1[n-1,n-7:] = np.array([ +120.0,  -864.0, +2700.0,  -4800.0,  +5400.0, -4320.0, +1764.0])/(720.0*h)
    # D^2 u_{n-1}   ~   D^2 u_6
    D2[n-1,n-7:] = np.array([ +274.0, -1944.0, +5940.0, -10160.0, +10530.0, -6264.0, +1624.0])/(360.0*h**2)

    # finite difference formula converts differential operators into matrices
    
    # radial part of Laplacian
    #  __2        1  d^2
    #  \/ f(r) = --- ---- (r f)
    #             r  dr^2
    #
    # d^2/dr^2
    d2dr2 = np.dot(np.diag(c1), D1) + np.dot(np.diag(c2), D2)
    #
    lap_rad =  np.dot(np.diag(1.0/r), np.dot(d2dr2, np.diag(r)))

    return lap_rad

def laplacian_matrix_spherical(Zat,
                               lebedev_order=23,
                               radial_grid_factor=1):
    """
    Matrix representation of the Laplacian operator

       __2    1  d^2         L^2
       \/  = --- ---- r  -  ----
              r  dr^2        r^2

    on a spherical grid for atom Z.
    """
    # radial grid
    Nr = number_of_radial_points(Zat)
    # increase number of grid points is requested
    Nr *= radial_grid_factor
    # scaling factor
    rm = 0.5*slater_radii[atom_names[Zat-1]]
    r = get_radial_grid(Nr, rm)

    # radial and angular parts of Laplacian
    lap_rad = radial_laplacian(Nr, rm)    
    lap_ang = angular_laplacian(lebedev_order)

    # number of angular grid points
    Nang = lap_ang.shape[0]
    #
    # total number of grid points
    N = Nr*Nang
    lap = np.zeros((N,N), dtype=complex)

    # enumerate columns
    cc = 0
    for ic in range(0, Nr):
        for jc in range(0, Nang):
            # enumerate rows
            rr = 0
            for ir in range(0, Nr):
                for jr in range(0, Nang):
                    # solid angle is the same, but radius may be different
                    if jr == jc:
                        lap[rr,cc] += lap_rad[ir,ic]
                    # radius is the same, but solid angle may be different
                    if ir == ic:
                        lap[rr,cc] -= 1.0/r[ir]**2 * lap_ang[jr,jc]
                    rr += 1
            cc += 1

    assert rr == N
    assert cc == N
    
    return lap

def apply_laplacian_matrix(Zat, x,
                    lebedev_order=23,
                    radial_grid_factor=1):
    """
    compute the matrix-vector product of the Laplacian operator
    with a vector

         y = L.x 
    """
    # radial grid
    Nr = number_of_radial_points(Zat)
    # increase number of grid points is requested
    Nr *= radial_grid_factor
    # scaling factor
    rm = 0.5*slater_radii[atom_names[Zat-1]]
    r = get_radial_grid(Nr, rm)

    # radial and angular parts of Laplacian
    lap_rad = radial_laplacian(Nr, rm)    
    lap_ang = angular_laplacian(lebedev_order)

    # number of angular grid points
    Nang = lap_ang.shape[0]
    #
    # total number of grid points
    N = Nr*Nang

    x = np.reshape(x, (Nr,Nang))
    # L.x
    y = np.zeros((Nr,Nang), dtype=complex)

    # apply angular part of Laplacian operator
    for i in range(0, Nr):
        y[i,:] -= 1.0/r[i]**2 * np.dot(lap_ang, x[i,:])
    # apply radial part of Laplacian operator
    for j in range(0, Nang):
        y[:,j] += np.dot(lap_rad, x[:,j])

    y = y.flatten()
    return y

        
def prune_grids(grid_points, grid_weights, grid_volumes, thresh=1.0e-5):
    """
    remove grid points with zero weights
    """
    pruned_points  = []
    pruned_weights = []
    pruned_volumes = []
    # list of arrays of indices into the original arrays, which identify
    # the points which remain
    active_indices = []
    
    # number of grids
    nc = len(grid_points)
    
    for i in range(0, nc):
        x,y,z = grid_points[i]
        dV = grid_volumes[i]
        w = grid_weights[i]
        # select points with weights > threshold
        pruned_points.append( (x[w > thresh], y[w > thresh], z[w > thresh]) )
        pruned_weights.append( w[w > thresh] )
        pruned_volumes.append( dV[w > thresh] )
        # keep track of indices into original array
        active_indices.append( np.where(w > thresh)[0] )
        
        
    return pruned_points, pruned_weights, pruned_volumes, active_indices

class MulticenterGrid(object):
    def __init__(self, atomlist):
        self.atomlist = atomlist
        # The grid resolution at the time when the creator is invoked
        # has to be save, since variables defined in settings might
        # change.
        self.lebedev_order = settings.lebedev_order
        self.radial_grid_factor = settings.radial_grid_factor
        # generate multicenter grid
        self.points, self.weights, self.volumes = \
            multicenter_grids(atomlist,
                              lebedev_order=self.lebedev_order,
                              radial_grid_factor=self.radial_grid_factor)
        # list with number of points in each grid
        self.num_points = [len(pts) for pts in self.weights]
        # I we want to separate the parts of a vector belonging to each
        # grid, the vector has to be split at these positions
        self.split_positions = np.cumsum(self.num_points)[:-1]
        # total number of points
        self.n = sum(self.num_points)
        self.shape = (self.n, self.n)
        
    def vector2grid(self, f):
        """
        split vector into pieces belonging to different grids and
        insert zeros for inactive points.
        """
        fs = np.split(f, self.split_positions, axis=0)
        return fs

    def grid2vector(self, fs):
        f = np.hstack(fs)
        return f

    def evaluate_function(self, func):
        """
        evaluate a function f(x,y,z) on the multicenter grid
        """
        fs = []
        for (xI,yI,zI) in self.points:
            fI = func(xI,yI,zI)
            fs.append(fI)
        f = self.grid2vector(fs)
        return f

    def scalar_product(self, f1, f2):
        # weights (volume element x fuzzy Voronoi weight function)
        w = np.hstack(self.volumes) * np.hstack(self.weights)
        return np.sum(w * f1*f2)
    
    def interpolation_func(self, f):
        """
        create a function that interpolates the grid values
        """
        fs = self.vector2grid(f)
        atomic_numbers, atomic_coordinates = atomlist2arrays(self.atomlist)
        func = multicenter_interpolation(fs,
                                         atomic_coordinates, atomic_numbers,
                                         lebedev_order=self.lebedev_order,
                                         radial_grid_factor=self.radial_grid_factor,
                                         weighted=True)

        return func

    def combine_subproblems(self, fs):
        atomic_numbers, atomic_coordinates = atomlist2arrays(self.atomlist)
        # create interpolation function
        func = multicenter_interpolation(fs,
                                         atomic_coordinates, atomic_numbers,
                                         lebedev_order=self.lebedev_order,
                                         radial_grid_factor=self.radial_grid_factor,
                                         weighted=False)
        fs_combined = []
        for (xI,yI,zI) in self.points:
            fI = func(xI,yI,zI)
            fs_combined.append(fI)

        return fs_combined
    
class LaplacianOperator(LinearOperator):
    dtype = np.dtype(complex)
    def __init__(self, atomlist):
        self.grid = MulticenterGrid(atomlist)
        # angular Laplacian is the same for all atoms
        self.lap_ang = angular_laplacian(settings.lebedev_order)
        self.Nang = self.lap_ang.shape[0]
        # radial grids are different for each atom type
        self.radial_laplacians = []
        self.radial_grids = []
        self.num_points_rad = []
        
        for Zat,pos in atomlist:
            # radial grid
            Nr = number_of_radial_points(Zat)
            # increase number of grid points is requested
            Nr *= settings.radial_grid_factor
            # scaling factor
            rm = 0.5*slater_radii[atom_names[Zat-1]]
            r = get_radial_grid(Nr, rm)
    
            # radial part of Laplacian
            lap_rad = radial_laplacian(Nr, rm)    

            self.radial_grids.append(r)
            self.radial_laplacians.append(lap_rad)
            self.num_points_rad.append(Nr)
            
        self.shape = self.grid.shape

    def _matvec(self, x):
        """
         y = L.x
        """
        # split vector into arrays belonging to each cell
        xs = self.grid.vector2grid(x)
        ys = []
        # apply Laplacian for each grid
        for I, xI in enumerate(xs):
            Nr = self.num_points_rad[I]
            Nang = self.Nang
            rI = self.radial_grids[I]
            lap_ang = self.lap_ang
            lap_rad = self.radial_laplacians[I]
            wI = self.grid.weights[I]

            # solve the subproblem
            #   (I)           (I)    (I)
            #  y     =  L . (w    * x   )
            xI *= wI
            
            xI = np.reshape(xI, (Nr,Nang))
            # L.x
            yI = np.zeros((Nr,Nang), dtype=complex)

            # apply angular part of Laplacian operator
            for i in range(0, Nr):
                yI[i,:] -= 1.0/rI[i]**2 * np.dot(lap_ang, xI[i,:])
            # apply radial part of Laplacian operator
            for j in range(0, Nang):
                yI[:,j] += np.dot(lap_rad, xI[:,j])

            yI = yI.flatten()
             
            ys.append(yI)

        #                (I)
        # add solutions y    from different subproblems
        ys = self.grid.combine_subproblems(ys)
            
        return self.grid.grid2vector(ys)
                      
class PotentialOperator(LinearOperator):
    dtype = np.dtype(complex)
    
    def __init__(self, atomlist, potential):
        self.grid = MulticenterGrid(atomlist)
        self.shape = self.grid.shape
        # evaluate potential on the grid points
        self.pots = []
        for (xI,yI,zI) in self.grid.points:
            potI = potential(xI,yI,zI)
            self.pots.append( potI )

    def _matvec(self, x):
        # split vector into arrays belonging to each cell
        xs = self.grid.vector2grid(x)
        ys = []
        # apply potential operator on each grid
        for I, xI in enumerate(xs):
            potI = self.pots[I]
            # potential energy is diagonal in position representation
            yI = potI * xI
            
            ys.append(yI)
            
        return self.grid.grid2vector(ys)


     
class HamiltonianOperator(LinearOperator):
    dtype = np.dtype(complex)
    def __init__(self, atomlist, potential):
        """
        create a linear operator that represents the Hamiltonian T + V
        on a multicenter grid

        Parameters
        ==========
        atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                        and positions
        potential    :  callable, potential(x,y,z) evaluates effective potential V
        """
        self.atomlist = atomlist
        # generate multicenter grid
        self.grid_points, self.grid_weights, self.grid_volumes = \
            multicenter_grids(atomlist,
                              lebedev_order=settings.lebedev_order,
                              radial_grid_factor=settings.radial_grid_factor)
        # number of grids
        self.num_grids = len(self.grid_weights)
        self.all_num_points = [len(pts) for pts in self.grid_weights]
        # prune grids, grid points with zero weights are removed
        self.grid_points, self.grid_weights, self.grid_volumes, self.active_indices = prune_grids(self.grid_points, self.grid_weights, self.grid_volumes)
        
        # list with number of points in each grid
        self.num_points = [len(pts) for pts in self.grid_weights]
        # I we want to separate the parts of a vector belonging to each
        # grid, the vector has to be split at these positions
        self.split_positions = np.cumsum(self.num_points)[:-1]
        # total number of points
        self.n = sum(self.num_points)
        self.shape = (self.n, self.n)
        
        # list of square matrices representing the kinetic energy operator
        # on each spherical grid
        self.grid_T_matrices = []
        # list of square matrices representing the potential energy
        # operator on each spherical grid
        self.grid_V_matrices = []
        # list of square matrices representing the Hamiltonian operator T+V
        # on each spherical grid
        self.grid_H_matrices = []
        # 
        self.grid_K_matrices = []
        
        # Now create matrices for T and V operators on each grid
        for I,(Z,pos) in enumerate(atomlist):
            # Laplacian matrix for spherical grid around atom I
            lap = laplacian_matrix_spherical(Z,
                                             lebedev_order=settings.lebedev_order,
                                             radial_grid_factor=settings.radial_grid_factor)
            lap = lap[:,self.active_indices[I]][self.active_indices[I],:]
            # kinetic energy operator on grid I
            T = -0.5 * lap
            # potential energy operator (only diagonal elements are stored)
            xI,yI,zI = self.grid_points[I]
            # evaluate potential function at the grid points belonging to atom I
            V = np.diag( potential(xI,yI,zI) )

            # Hamiltonian operator on grid I
            H = T + V

            # Build matrix for non-Hermitian eigenvalue problem
            #
            # All spherical grids cover the entire space, however with
            # very different resolution for different regions. The multicenter
            # grid is a superposition of all spherical grids, but there is no
            # coupling between the individual grids.
            # The weight factor of a grid point is the product of the volume
            # element and the value of the weight function due to the fuzzy
            # Voronoi decomposition.
            w = self.grid_volumes[I] * self.grid_weights[I]

            # We wish to minimize the expectation value of the energy
            # subject to the constraint that the wavefunction is normalized,
            #     
            #     minimize  <phi|H|phi>   s.t.   <phi|phi> = 1
            #
            # The values of the wavefunction on the grid are f_i = phi(x_i, y_i, z_i).
            # The volume element (weight) of each grid point is w_i. In terms of the
            # values and weights, the minimization problem becomes:
            #                                                           2
            #    minimize  sum    w  f  H    f           s.t.   sum  w f  = 1
            #       f_i       i,j  i  i  i,j  j                    i  i i
            #
            # Minimization problems subject to equality constraints can be solved by
            # introducing Lagrange multipliers. The new objective function becomes
            #                                                  2
            #    L(f) = sum    w  f  H    f   -  E * (sum  w  f   -  1)
            #              i,j  i  i  i,j  j             i  i  i
            #
            # At the minimum we have dL/df_k = 0. When deriving the eigenvalue equation,
            # we have to keep in mind that the Hamiltonian is not symmetric (!), so H_ij != H_ji^*:
            #                       
            #    sum  1/2 (H     +  w / w * H    ) f   =  E f
            #       j       i,j      j   i   j,i    j        i
            #
            # This is a non-Hermitian eigenvalue problem of the form
            #
            #    K.f = E*f
            #

            # Build matrix for non-Hermitian eigenvalue problem
            #  H_ij
            H1 = H
            #  w_i / w_j * H_ij
            H2 = np.multiply.outer(w,1.0/w) * H
            # 1/2 (H_ij + w_j / w_i * H_ji )
            K = 0.5 * (H1 + H2.transpose())
            
            self.grid_T_matrices.append( T )
            self.grid_V_matrices.append( V )
            self.grid_H_matrices.append( H )
            self.grid_K_matrices.append( K )

    def vector2grid(self, x):
        """
        split vector into pieces belonging to different grids and
        insert zeros for pruned points.
        """
        vecs = np.split(x, self.split_positions, axis=0)
        # For grid points with zero weight, the value 0 is filled in.
        vals = [np.zeros(n) for n in self.all_num_points]
        for i in range(0, self.num_grids):
            vals[i][self.active_indices[i]] = vecs[i]

        return vals        
    
    def _matvec(self, x):
        """
         y = K.x
        """
        # split vector into arrays belonging to each cell
        xs = np.split(x, self.split_positions, axis=0)
        ys = []
        for I,xi in enumerate(xs):
            # apply K matrix of I-th Voronoi cell to I-th
            # part of the vector
            K = self.grid_K_matrices[I]
            #
            yi = np.dot(K, xi)
            ys.append(yi)

        # join vectors from all cells
        y = np.hstack(ys)
    
        return y

    def _rmatvec(self, y):
        """
               H
          x = K .y
        """
        # split vector into arrays belonging to each cell
        ys = np.split(y, self.split_positions, axis=0)
        xs = []
        for I,yi in enumerate(ys):
            # apply K^H matrix of I-th Voronoi cell to I-th
            # part of the vector
            Kh = self.grid_K_matrices[I].conjugate().transpose()
            xi = np.dot(Kh, yi)
            xs.append(xi)

        # join vectors from all cells
        x = np.hstack(xs)

        return x
    
    def dense_matrix(self):
        """
        construct a dense matrix representation from the K's
        for each cell. The K matrix is block-diagonal since
        different cells/centers are not coupled.

            (  K1              )
            (     K2     0     ) 
        K = (        K3        )
            (    0      ...    )
            (              Kn  )

        """
        K = sla.block_diag(*self.grid_K_matrices)
        return K
        
    
        
###############################################################################
#
#   Testing
#  
###############################################################################


def test_angular_laplacian():

    lebedev_order = 7
    
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    lap_ang = angular_laplacian(lebedev_order)

    print "laplacian matrix"
    print lap_ang
    
    f = np.sin(th)*np.cos(ph)
    Lf_exact = 2.0*np.sin(th)*np.cos(ph)
    
    # apply Laplacian to f
    Lf = np.dot(lap_ang, f)

    err = la.norm(Lf - Lf_exact)
    
    print Lf
    print Lf_exact
    print "|Lf - Lf (exact)|= %e" % err
    assert err < 1.0e-8

def test_radial_laplacian():
    # radial grid
    Nr = 30
    rm = 0.5 * 0.47
    
    r = get_radial_grid(Nr, rm)

    # compute radial part of Laplacian of this function
    f = np.exp(-r)
    Lf_exact = 1.0/r * np.exp(-r) * (r - 2.0)

    lap_rad = radial_laplacian(Nr, rm)
    Lf = np.dot(lap_rad, f)
    
    print Lf
    print Lf_exact
    
    import matplotlib.pyplot as plt
    plt.plot(r, f, label="f(r)")
    plt.plot(r, Lf, label=r"$\nabla^2 f$ (numerical)")
    plt.plot(r, Lf_exact, ls="-.", label=r"$\nabla^2 f$ (exact)")
    plt.legend()
    plt.show()

def test_radial_laplacian_symmetry():
    """
    Is the radial laplacian matrix symmetric?  
    The answer is no!
    """
    from DFTB import utils
    import matplotlib.pyplot as plt
    
    # radial grid
    Nr = 7
    rm = 0.5 * 0.47
    
    # Laplacian matrix without boundary conditions
    lap_rad = radial_laplacian(Nr, rm)
    # Laplacian matrix with implicit boundary conditions
    #lap_rad_boundary0 = radial_laplacian_boundary0(Nr, rm)

    # print Laplacian matrix
    print "Laplacian matrix (no boundary conditions)"
    labels = ["%3.1d" % i for i in range(0, Nr)]
    print utils.annotated_matrix(lap_rad, labels, labels, format="%+3.3e")

    #print "Laplacian matrix (implicit boundary conditions)"
    ## print symmetric Laplacian matrix (with implicit boundary conditions)
    #labels = ["%3.1d" % i for i in range(0, Nr)]
    #print utils.annotated_matrix(lap_rad_boundary0, labels, labels, format="%+3.3e")

    
    # plot Laplacian matrix

    plt.imshow(lap_rad)
    plt.show()

    
def test_laplacian():
    """
    The Laplacian is computed numerically on a spherical grid
    and compared with the exact result.
    """
    Zat = 1
    atomlist = [(Zat,(0.0, 0.0, 0.0))]

    grid_points, grid_weights, grid_volumes = multicenter_grids(atomlist,
                                                                lebedev_order=settings.lebedev_order,
                                                                radial_grid_factor=settings.radial_grid_factor)
    lap = laplacian_matrix_spherical(Zat,
                                     lebedev_order=settings.lebedev_order,
                                     radial_grid_factor=settings.radial_grid_factor)
    print "size of matrices"
    print lap.shape

    """
    # test function
    def func(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        return np.exp(-r)

    # known Laplacian
    def lap_func_exact(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        return 1.0/r * np.exp(-r) * (r - 2.0)        
    """
    # test function
    def func(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        f = (x+y**2+z**3)*np.exp(-r**2)
        return f

    # known Laplacian
    def lap_func_exact(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        lapf = 2*np.exp(-r**2) * (1.0 + 2*x**3 +2*y**4 + 3*z - 9*z**3 + 2*z**5 + x*(-5+2*y**2+2*z**2) + 2*x**2 * (y**2+z**3) + y**2*(-7+2*z**2*(1+z)))
        return lapf

    x,y,z = grid_points[0]
    # The weight of a point is the product of the weight
    # from the Voronoi decomposition and the volume element.
    w = grid_weights[0] * grid_volumes[0]

    f = func(x,y,z)
    Lf_exact = lap_func_exact(x,y,z)

    Lf = np.dot(lap, f)

    # test implementation of matrix-vector product L.f
    Lf_matvec = apply_laplacian_matrix(Zat, f, 
                                       lebedev_order=settings.lebedev_order,
                                       radial_grid_factor=settings.radial_grid_factor)
    err = la.norm(Lf - Lf_matvec)
    print "|L.f - L.f(matvec)|= %e" % err
    assert err < 1.0e-7
    
    # interpolate Lf
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    lap_func_interp = multicenter_interpolation([Lf], atomic_coordinates, atomic_numbers,
                                          lebedev_order=settings.lebedev_order, radial_grid_factor=settings.radial_grid_factor)

    print Lf
    print Lf_exact

    dif = abs(Lf_exact - Lf)
    print dif
    
    import matplotlib.pyplot as plt
    r = np.linspace(-5.0, 5.0, 1000)
    x = 0*r
    y = 0*r
    z = r

    plt.plot(r, func(x,y,z), label=r"$f$")
    plt.plot(r, lap_func_exact(x,y,z), label=r"$\nabla^2 f$ (exact)")
    plt.plot(r, lap_func_interp(x,y,z), ls="-.", label=r"$\nabla^2 f$ (numerical)")

    plt.legend()
    plt.plot()
    plt.show()

def test_multicenter_laplacian():
    """
    compute the Laplacian of the LCAO wavefunction of the hydrogen
    molecule.
    """
    # bond length in bohr 
    R = 2.0
    # 
    atomlist = [(1, (0,0,-R/2.0)),
                (1, (0,0,+R/2.0))]

    def psi_1s(x,y,z):
        """wavefunction of 1s hydrogen electron"""
        r = np.sqrt(x*x+y*y+z*z)
        psi = 1.0/np.sqrt(np.pi) * np.exp(-r)
        return psi

    def lap_1s(x,y,z):
        """analytical Laplacian for 1s orbital"""
        r = np.sqrt(x*x+y*y+z*z)
        lap = 1.0/np.sqrt(np.pi) * np.exp(-r) * (1 - 2.0/r)
        return lap
    
    def psi_sigma(x,y,z):
        """unnormalized LCAO wavefunction of H-H"""
        # sigma orbital is a linear combination of two
        # 1s orbitals of hydrogen
        psi = psi_1s(x,y,z-R/2.0) + psi_1s(x,y,z+R/2.0)
        return psi

    def lap_sigma_exact(x,y,z):
        """                       __2
        analytical expression for \/  psi_sigma
        """
        return lap_1s(x,y,z-R/2.0) + lap_1s(x,y,z+R/2.0)

    # set radial and angular resolution of grid
    settings.radial_grid_factor=40
    settings.lebedev_order=41

    lap_op = LaplacianOperator(atomlist)

    # evaluate wavefunction on the grid
    psi_sigma_vec = lap_op.grid.evaluate_function(psi_sigma)

    print "matvec L.x"
    lap_sigma_vec = lap_op.matvec(psi_sigma_vec)

    print "interpolation"
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    lap_sigma_grid = multicenter_interpolation(lap_op.grid.vector2grid(lap_sigma_vec),
                                               atomic_coordinates, atomic_numbers,
                                               lebedev_order=settings.lebedev_order,
                                               radial_grid_factor=settings.radial_grid_factor)

    print "Laplacian on Becke grid"
    from DFTB.MolecularIntegrals.MulticenterIntegration import multicenter_laplacian
    lap_sigma_num = multicenter_laplacian(psi_sigma, atomic_coordinates, atomic_numbers,
                                          cusps_separate=True,
                                          radial_grid_factor=settings.radial_grid_factor,
                                          lebedev_order=settings.lebedev_order)

    
    # plot Laplacian along the H-H bond (the z-axis)
    import matplotlib.pyplot as plt
    plt.xlabel("z / bohr", fontsize=17)

    r = np.linspace(-10.0, 10.0, 100000)
    x, y, z = 0*r, 0*r, r

    plt.plot(r, psi_sigma(x,y,z), label=r"wavefunction $\psi_{\sigma}$")
    plt.plot(r, lap_sigma_grid(x,y,z),           label=r"$\nabla^2 \psi$, grid solution")
    plt.plot(r, lap_sigma_num(x,y,z), ls="-.",   label=r"$\nabla^2 \psi$, numerical solution")
    plt.plot(r, lap_sigma_exact(x,y,z), ls="--", label=r"$\nabla^2 \psi$, analytical solution")

    plt.legend()
    plt.show()


from DFTB.MolecularIntegrals import chebyshev
from DFTB.MolecularIntegrals.SphericalCoords import cartesian2spherical
from scipy import interpolate
def multicenter_grids(atomlist,
                      kmax=3,
                      lebedev_order=23, radial_grid_factor=1):
    """
    compute grid points and weights of the multicenter grids for visualization
   
    Parameters
    ----------
    atomlist           : list of tuples (Zat,(xI,yI,zI)) with atomic numbers and 
                         atom positions, which define the multicenter grid
    
    Optional
    --------
    kmax               : How fuzzy should the Voronoi polyhedrons be? Larger kmax
                         means borders are fuzzier.
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    grid_points        : list of tuples (x,y,z) with positions of points in each grid,
                         grid_points[I][0] contains the x-positions of the points
                         belonging to the grid around atom I
    grid_weights       : list of numpy arrays, grid_weights[I][k] contains the weight
                         of the k-th point in the grid around atom I due to the fuzzy
                         Voronoi decomposition.
    grid_volumes       : list of numpy arrays, grid_volumes[I][k] contains the volume
                         element around the k-th point in the grid at atom I.
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=kmax):
        f = mu
        for ik in range(0, k):
            f = 1.5 * f -0.5 * f**3
        return 0.5*(1-f)

    plot_cutoff_profile = False
    if plot_cutoff_profile == True:
        import matplotlib.pyplot as plt
        mu = np.linspace(-1.0,1.0,100)
        for k in range(1,5):
            plt.plot(mu, s(mu,k=k), label=r"$k=%d$" % k)
        plt.legend()
        plt.show()
    
    atomic_names = [atom_names[Z-1] for Z in atomic_numbers]
    
    Nat = atomic_coordinates.shape[1]
    R = np.zeros((Nat,Nat))     # distances between atoms i and j
    a = np.zeros((Nat,Nat))     # scaling factor used in eqn. A2
    for i in range(0, Nat):
        for j in range(i+1, Nat):
            R[i,j] = la.norm(atomic_coordinates[:,i] - atomic_coordinates[:,j])
            R[j,i] = R[i,j]

            # ratio of Slater radii
            chi = slater_radii[atomic_names[i]] / slater_radii[atomic_names[j]]
            uij = (chi-1)/(chi+1)
            a[i,j] = uij/(uij**2 - 1)
            a[j,i] = -a[i,j]

    grid_points = []
    grid_weights = []
    grid_volumes = []
    # atom-centered subintegral
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = (k-0.5)/Nr
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/Nr * np.sin(zr * np.pi)**2
        # from variable transformation
        g = 2 * rm**3 * np.sqrt(((1+xr)/(1-xr)**3)**3)
        radial_weights *= g
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # cartesian coordinates of grid
        x = (outerN(r, sc) + atomic_coordinates[0,I]).flatten()
        y = (outerN(r, ss) + atomic_coordinates[1,I]).flatten()
        z = (outerN(r, c ) + atomic_coordinates[2,I]).flatten()
        weights = outerN(radial_weights, 4.0*np.pi * angular_weights).flatten()
        #
        Npts = Nr*Nang
        # distance between grid points and atom i
        dist = np.zeros((Npts, Nat))
        for i in range(0, Nat):
            dist[:,i] = np.sqrt(    (x - atomic_coordinates[0,i])**2   \
                                   +(y - atomic_coordinates[1,i])**2   \
                                   +(z - atomic_coordinates[2,i])**2 )

        # P_i(r) as defined in eqn. (13)
        P = np.ones((Npts,Nat))
        for i in range(0, Nat):
            for j in range(0, Nat):
                if i==j:
                    continue
                # mu_ij as defined in eqn. (11)
                mu = (dist[:,i]-dist[:,j])/R[i,j]
                nu = mu + a[i,j]*(1-mu**2)
                P[:,i] *= s(nu)
        Ptot = np.sum(P, axis=1)
    
        # weight function due to partitioning of volume
        wr = P[:,I]/Ptot
        
        grid_points.append( (x.flatten(), y.flatten(), z.flatten()) )
        # The weights come from the fuzzy Voronoi partitioning 
        grid_weights.append( wr.flatten() )
        # The naming is a little bit confusing, the `weights` are
        # actually the volume elements dV_i around each point.
        grid_volumes.append( weights.flatten() )

    return grid_points, grid_weights, grid_volumes

def multicenter_interpolation(grid_values,
                              atomic_coordinates, atomic_numbers,
                              kmax=3,
                              lebedev_order=23, radial_grid_factor=1,
                              weighted=True):
    """
    Given the values at the grid points, an interpolation function is created.

    Parameters
    ----------
    grid_values        : list of tuples (fx,fy,fz) with values of function at the points
                         in each grid, 
                         grid_values[I] contains the function values on the grid belonging
                         to the I-th atom, the ordering of the points is the same as in the
                         grids returned by `multicenter_grids`.
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    kmax               : How fuzzy should the Voronoi polyhedrons be? Larger kmax
                         means borders are fuzzier.
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor
    weighted           : flag which determines if the function values
                         should be premultiplied with the weight function

    Returns
    -------
    f                  : callable, f(x,y,z) interpolates the grid data
    """    
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=kmax):
        ff = mu
        for ik in range(0, k):
            ff = 1.5 * ff -0.5 * ff**3
        return 0.5*(1-ff)
    
    atomic_names = [atom_names[Z-1] for Z in atomic_numbers]
    
    Nat = atomic_coordinates.shape[1]
    R = np.zeros((Nat,Nat))     # distances between atoms i and j
    a = np.zeros((Nat,Nat))     # scaling factor used in eqn. A2
    for i in range(0, Nat):
        for j in range(i+1, Nat):
            R[i,j] = la.norm(atomic_coordinates[:,i] - atomic_coordinates[:,j])
            R[j,i] = R[i,j]

            # ratio of Slater radii
            chi = slater_radii[atomic_names[i]] / slater_radii[atomic_names[j]]
            uij = (chi-1)/(chi+1)
            a[i,j] = uij/(uij**2 - 1)
            a[j,i] = -a[i,j]

    radial_functions = []
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        interp = chebyshev.ChebyshevInterpolator(Nr)
        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = (k-0.5)/Nr
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/Nr * np.sin(zr * np.pi)**2
        # from variable transformation
        gg = 2 * rm**3 * np.sqrt(((1+xr)/(1-xr)**3)**3)
        radial_weights *= gg
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # cartesian coordinates of grid
        x = (outerN(r, sc) + atomic_coordinates[0,I])
        y = (outerN(r, ss) + atomic_coordinates[1,I])
        z = (outerN(r, c ) + atomic_coordinates[2,I])
        weights = outerN(radial_weights, 4.0*np.pi * angular_weights)
        #
        Npts = Nr*Nang

        # grid_values[I] should correspond to f(x,y,z), we need to
        # check that the dimensions match.
        assert len(grid_values[I]) == Nr*Nang, "Size of grid data %d  does not match the size of grid  %d*%d = %d!" % (len(grid_values[I]), Nr, Nang, Nr*Nang)
        # We don't need to evaluate the functions fi on the grid,
        # we already have the function values
        fI = np.reshape(grid_values[I], (Nr,Nang))

        # multiply with weight function to get
        #   wI(x,y,z)*f(x,y,z)
        # if desired
        if weighted:
            # distance between grid points and atom i
            dist = np.zeros((Nr,Nang, Nat))
            for i in range(0, Nat):
                dist[:,:,i] = np.sqrt(   (x - atomic_coordinates[0,i])**2   \
                                        +(y - atomic_coordinates[1,i])**2   \
                                        +(z - atomic_coordinates[2,i])**2 )

            # P_i(r) as defined in eqn. (13)
            P = np.ones((Nr,Nang,Nat))
            for i in range(0, Nat):
                for j in range(0, Nat):
                    if i==j:
                        continue
                    # mu_ij as defined in eqn. (11)
                    mu = (dist[:,:,i]-dist[:,:,j])/R[i,j]
                    nu = mu + a[i,j]*(1-mu**2)
                    P[:,:,i] *= s(nu)
            Ptot = np.sum(P, axis=-1)
    
            # weight function
            wr = P[:,:,I]/Ptot

            fI *= wr            
        
        radial_functions.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            fI_lm = 4.0*np.pi * np.sum(fI*wYlm, axis=-1)

            #spline_lm_real = interpolate.splrep(zr, fI_lm.real, s=0)
            #spline_lm_imag = interpolate.splrep(zr, fI_lm.imag, s=0)
            ### DEBUG
            # linear interpolation
            spline_lm_real = interpolate.interp1d(zr, fI_lm.real, fill_value='extrapolate', bounds_error=False, kind=3)
            spline_lm_imag = interpolate.interp1d(zr, fI_lm.imag, fill_value='extrapolate', bounds_error=False, kind=3)
            #spline_lm_real = interpolate.PchipInterpolator(zr, fI_lm.real, extrapolate=True)
            #spline_lm_imag = interpolate.PchipInterpolator(zr, fI_lm.imag, extrapolate=True)
            #spline_lm_real = interp.expansion_coefficients(fI_lm.real)
            #spline_lm_imag = interp.expansion_coefficients(fI_lm.imag)
            ###
            radial_functions[-1][(l,m)] = spline_lm_real, spline_lm_imag

            if m == -(Lmax-1)/2:
                break

    def interpolation_func(x,y,z):
        """
        function for interpolating f
        """
        f = 0j*x
        # sum over centers
        #  f = sum_I  f^(I)
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            Nr = number_of_radial_points(atomic_numbers[I])
            # increase number of grid points is requested
            Nr *= radial_grid_factor
            interp = chebyshev.ChebyshevInterpolator(Nr)
            
            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:
                
                spline_lm_real, spline_lm_imag = radial_functions[I][(l,m)]
                # interpolate
                #fI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                #        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)
                ### DEBUG
                #fI_lm = interp.interpolate(spline_lm_real, xr) + 1.0j*interp.interpolate(spline_lm_imag, xr)
                fI_lm = spline_lm_real(zr) + 1.0j * spline_lm_imag(zr)
                #if l == 0:
                #    import matplotlib.pyplot as plt
                #    plt.cla()
                #plt.plot(zr, fI_lm.real, label=r"l=%d m=%d  real" % (l,m))
                #plt.plot(zr, fI_lm.imag, label=r"l=%d m=%d  imag." % (l,m))
                #if m == -(Lmax-1)/2:
                #    plt.legend()
                #    plt.show()

                ###
                
                f += fI_lm*Ylm

                if m == -(Lmax-1)/2:
                    break

        return f.real
    
    return interpolation_func


    

def test_interpolation():
    """
    check whether repeated cycles of sampling/interpolation on the
    grid change the function, ideally the function should remain
    the same.
    """
    R = 2.0
    atomlist = [(1,(0.0, 0.0, -R/2.0)),
                (1,(0.0, 0.0, +R/2.0))]
    #atomlist = [(1,(0.0, 0.0, 0.0))]

    # set radial and angular resolution of grid
    settings.radial_grid_factor=40 #5 #1 #40
    settings.lebedev_order=41   # 23 # 41
    
    def psi_1s(x,y,z):
        #wavefunction of 1s hydrogen electron
        r = np.sqrt(x*x+y*y+z*z)
        psi = 1.0/np.sqrt(np.pi) * np.exp(-r)
        return psi

    def psi_sigma(x,y,z):
        # unnormalized LCAO wavefunction of H-H
        # sigma orbital is a linear combination of two
        # 1s orbitals of hydrogen
        psi = psi_1s(x,y,z-R/2.0) + psi_1s(x,y,z+R/2.0)
        return psi

    def rho0(x,y,z):
        #return psi_1s(x,y,z)**2
        return 2*psi_sigma(x,y,z)**2
    
    import matplotlib.pyplot as plt
    plt.cla()
    plt.clf()
    # grid for plotting
    r = np.linspace(-15.0, 15.0, 100000)
    x = 0*r
    y = 0*r
    z = r
    
    grid = MulticenterGrid(atomlist)

    rho = rho0
    for i in range(0, 20):
        print "%d x sampling/interpolating" % i
        plt.plot(r, rho(x,y,z), label=r"$\rho$ %d x interp." % i)

        f = grid.evaluate_function(rho)
        # interpolate
        rho = grid.interpolation_func(f)

        # compute deviation from original function
        err = np.sqrt(integral(atomlist, lambda x,y,z: abs(rho(x,y,z) - rho0(x,y,z))**2))

        print "  i= %d     |f_i - f|= %e" % (i, err)
        
    plt.legend()
    plt.show()
    
    
def test_possion_solver_sparse():
    """
    determine the electrostatic potential of an electronic density
    using a sparse solver
    """
    R = 2.0
    atomlist = [(1,(0.0, 0.0, -R/2.0)),
                (1,(0.0, 0.0, +R/2.0))]

    def psi_1s(x,y,z):
        #wavefunction of 1s hydrogen electron
        r = np.sqrt(x*x+y*y+z*z)
        psi = 1.0/np.sqrt(np.pi) * np.exp(-r)
        return psi

    def psi_sigma(x,y,z):
        # unnormalized LCAO wavefunction of H-H
        # sigma orbital is a linear combination of two
        # 1s orbitals of hydrogen
        psi = psi_1s(x,y,z-R/2.0) + psi_1s(x,y,z+R/2.0)
        return psi

    def rho(x,y,z):
        return 2*psi_sigma(x,y,z)**2

    Lop = LaplacianOperator(atomlist)

    from scipy.sparse.linalg import lsmr

    # solve
    #   __2
    #   \/  V = rho
    print "solve Poisson equation..."
    y = Lop.grid.evaluate_function(rho)
    ret = lsmr(Lop, y)
    x = ret[0]

    V_func = Lop.grid.interpolation_func(x)

    import matplotlib.pyplot as plt
    r = np.linspace(-10.0, 10.0, 10000)
    x = 0*r
    y = 0*r
    z = r

    plt.plot(r, rho(x,y,z), label=r"$\rho$")
    plt.plot(r, V_func(x,y,z), label=r"$V$")
    plt.legend()
    plt.show()
    
def test_integration():
    """
    Suppose that f(i) are the function values of f() at the
    grid points. The integral of f is then given by
     
         sum_i weights(i) f(i)
    """
    from DFTB.MolecularIntegrals.Ints1e import integral

    Zat = 1
    atomlist = [(Zat,(0.0, 0.0, 0.0))]

    grid_points, grid_weights, grid_volumes = multicenter_grids(atomlist,
                                                                lebedev_order=settings.lebedev_order,
                                                                radial_grid_factor=settings.radial_grid_factor)

    def f(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        return np.exp(-r)

    # compute integral on the grid
    xi,yi,zi = grid_points[0]
    # The weight of a point is the product of the weight
    # from the Voronoi decomposition and the volume element.
    wi = grid_weights[0] * grid_volumes[0]
    fi = f(xi,yi,zi)

    integ_grid = np.sum(wi*fi)

    #
    integ = integral(atomlist, f)

    print "integral on the grid"
    print " I = %e" % integ_grid
    print "integral using Becke's integration scheme"
    print " I = %e" % integ

    assert abs(integ_grid - integ) < 1.0e-10

def test_hydrogen_spherical_grid():
    """
    compute lowest eigenstate of the hydrogen atom by direct diagonalization
    of the Hamiltonian on a spherical grid
    """
    import matplotlib.pyplot as plt

    # grid for plotting
    r = np.linspace(-5.0, 5.0, 100000)
    x = 0*r
    y = 0*r
    z = r

    # reduce resolution of grid, otherwise the matrices become too large
    settings.radial_grid_factor = 10
    settings.lebedev_order = 3 #11
    
    Zat = 1
    atomlist = [(Zat,(0.0, 0.0, 0.0))]
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    # Coulomb potential
    def potential(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        return -1.0/r
    
    grid_points, grid_weights, grid_volumes = multicenter_grids(atomlist,
                                                                lebedev_order=settings.lebedev_order,
                                                                radial_grid_factor=settings.radial_grid_factor)
    
    # grid points (we have only a single grid here)
    xg,yg,zg = grid_points[0]
    # The weight of a point is the product of the weight
    # from the Voronoi decomposition and the volume element.
    w = grid_weights[0] * grid_volumes[0]

    # evaluate potential on the grid, diagonal elements of matrix
    Vii = potential(xg,yg,zg)
    # potential matrix
    V = np.diag(Vii)
    
    lap = laplacian_matrix_spherical(Zat,
                                     lebedev_order=settings.lebedev_order,
                                     radial_grid_factor=settings.radial_grid_factor)

    #                        __2
    # kinetic energy is -1/2 \/
    T = -0.5*lap

    # Hamiltonian is non-Hermitian, i.e. Hij != Hji^*  
    H = T+V

    # We wish to minimize the expectation value of the energy
    # subject to the constraint that the wavefunction is normalized,
    #     
    #     minimize  <phi|H|phi>   s.t.   <phi|phi> = 1
    #
    # The values of the wavefunction on the grid are f_i = phi(x_i, y_i, z_i).
    # The volume element (weight) of each grid point is w_i. In terms of the
    # values and weights, the minimization problem becomes:
    #                                                           2
    #    minimize  sum    w  f  H    f           s.t.   sum  w f  = 1
    #       f_i       i,j  i  i  i,j  j                    i  i i
    #
    # Minimization problems subject to equality constraints can be solved by
    # introducing Lagrange multipliers. The new objective function becomes
    #                                                  2
    #    L(f) = sum    w  f  H    f   -  E * (sum  w  f   -  1)
    #              i,j  i  i  i,j  j             i  i  i
    #
    # At the minimum we have dL/df_k = 0. When deriving the eigenvalue equation,
    # we have to keep in mind that the Hamiltonian is not symmetric (!), so H_ij != H_ji^*:
    #                       
    #    sum  1/2 (H     +  w / w * H    ) f   =  E f
    #       j       i,j      j   i   j,i    j        i
    #
    # This is a non-Hermitian eigenvalue problem of the form
    #
    #    K.f = E*f

    # Build matrix for non-Hermitian eigenvalue problem
    #  H_ij
    H1 = H
    #  w_i / w_j * H_ij
    H2 = np.multiply.outer(w,1.0/w) * H
    # 1/2 (H_ij + w_j / w_i * H_ji )
    K = 0.5 * (H1 + H2.transpose())

    """
    # solve non-Hermitian eigenvalue problem
    #  K.f = E f
    # For some reason the direct eigenvalue solver produces garbage
    print "diagonalizing %d x %d dimensional non-symmetric matrix" % K.shape
    eigvals, eigvecs = la.eig(K)
    """
    # solve for lowest few eigenvectors iteratively
    from scipy.sparse.linalg import eigs
    print "solve for lowest eigenvectors using iterative algorithm"
    en_guess = -0.5
    eigvals, eigvecs = eigs(K, k=8, sigma=en_guess)
    print eigvals
    # The number of degenerate eigenstates of the hydrogen atom
    # should be
    #
    #   n     En = -1/(2n^2)      degeneracy
    # -----------------------------------------
    #   1       -1/2                 1
    #   2       -1/8                 3+1=4
    #   3       -1/18                5+3+1=9
    #  etc.
    #
    # However, the eigenvalue solver produces vectors that are
    # linearly dependent. For instance, if the order of the Lebedev
    # grid is large enough (11), there are several eigenvectors
    # with eigenvalue -1/2, although the true ground state is not degenerate.
    #
    
    # eigenvectors are sorted in increasing order by eigenvalues
    sort_indx = np.argsort(eigvals)
    eigvals = eigvals[sort_indx]
    eigvecs = eigvecs[:,sort_indx]

    # overlap matrix between eigenvectors
    S = np.dot(eigvecs.transpose(), np.dot(np.diag(w), eigvecs))
    print "overlap matrix between eigenstates"
    print S
    
    
    # lowest eigenvector
    for i in range(0, 6):
        f = eigvecs[:,i]
    
        # The eigenvalue solver produces eigenvectors that are normalized as
        #         2
        #   sum  f  = 1
        #      i  i
        #
        # However, the correct normalization should contain the weights:
        #            2
        #   sum  w  f  = 1
        #      i  i  i
        #
        # Therefore we need to renormalize the f_i's
        f /= np.sqrt( np.sum(w * f**2) )

        print "normalization"
        norm2 = np.sum(w * f**2)
        print np.sqrt(norm2)
    
        en = np.sum(w * f * np.dot(T+V, f)) 
        print "energy expectation value of %d-th eigenvector" % i
        print en

        enK = np.sum(w * f * np.dot(K, f))
        print "(w*f)^T.K.f"
        print enK
    
        # plot function belonging to lowest eigenvalue
        psi0 = multicenter_interpolation([f],
                                         atomic_coordinates, atomic_numbers,
                                         lebedev_order=settings.lebedev_order,
                                         radial_grid_factor=settings.radial_grid_factor)
    
        plt.xlabel("z / bohr")
        plt.plot(r, psi0(x,y,z), label=r"$\psi_0$ ($E=%+5.5f$)" % en)

    def psi_1s(x,y,z):
        """wavefunction of 1s hydrogen electron"""
        r = np.sqrt(x*x+y*y+z*z)
        psi = 1.0/np.sqrt(np.pi) * np.exp(-r)
        return psi

    f_exact = psi_1s(xg,yg,zg)

    print "normalization of exact solution"
    norm2 = np.sum(w * f_exact**2)
    print np.sqrt(norm2)
    
    en = np.sum(w * f_exact * np.dot(T+V, f_exact)) 
    print "energy expectation value of exact solution"
    print en

    enK = np.sum(w * f_exact * np.dot(K, f_exact))
    print "(w*f)^T.K.f (exact)"
    print enK

    
    # plot exact eigenfunctions
    psi_exact = multicenter_interpolation([f_exact],
                                          atomic_coordinates, atomic_numbers,
                                          lebedev_order=settings.lebedev_order,
                                          radial_grid_factor=settings.radial_grid_factor)
    
    plt.plot(r, psi_exact(x,y,z), ls="--", label=r"$\psi_0$ (exact)")
    
    plt.legend()
    plt.show()

def test_linear_operator_hydrogen():
    """
    compute lowest eigenstate of the hydrogen atom by direct diagonalization
    of the Hamiltonian on a spherical grid
    """
    # reduce resolution of grid, otherwise the matrices become too large
    settings.radial_grid_factor = 10
    settings.lebedev_order = 3 #11
    
    Zat = 1
    atomlist = [(Zat,(0.0, 0.0, 0.0))]

    # Coulomb potential
    def potential(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        return -1.0/r

    Hop = HamiltonianOperator(atomlist, potential)

    from scipy.sparse.linalg import eigs
    print "solve for lowest eigenvectors using iterative algorithm"
    en_guess = -0.5

    # check the identity   (K.u)^H.v = u^H.(K^H.v)
    n = Hop.shape[0]
    # random vector u and v
    u = np.random.rand(n)
    v = np.random.rand(n)
    #  (K.u)^H.v
    p1 = np.dot(Hop._matvec(u).conjugate().transpose(), v)
    #  u^H.(K^h.v)
    p2 = np.dot(u.conjugate().transpose(), Hop._rmatvec(v))
    assert abs(p1-p2)/abs(p1) < 1.0e-10

    eigvals, eigvecs = eigs(Hop, k=8, sigma=en_guess)
    
    #K = Hop.dense_matrix()
    #eigvals, eigvecs = eigs(K, k=8, sigma=en_guess)
    print eigvals

def test_multicenter_hmi():
    """
    compute lowest eigenstate of the hydrogen molecular ion (HMI)
    by direct diagonalization of the Hamiltonian on two overlapping
    atom-centered spherical grids
    """
    # reduce resolution of grid, otherwise the matrices become too large
    settings.radial_grid_factor = 20
    settings.lebedev_order = 41

    R = 2.0
    atomlist = [(1,(0.0, 0.0, -R/2.0)),
                (1,(0.0, 0.0, +R/2.0))]
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # only Coulomb potential
    def potential(x,y,z):
        Vnuc = 0.0*x
        for (Zi,posi) in atomlist:            
            Vnuc += - Zi / np.sqrt( (x-posi[0])**2 + (y-posi[1])**2 + (z-posi[2])**2 )
        return Vnuc

    Lop = LaplacianOperator(atomlist)
    grid = Lop.grid    

    # kinetic energy operator
    Top = -0.5 * Lop
    # potential energy operator
    Vop = PotentialOperator(atomlist, potential)
    # Hamiltonian operator
    Hop = Top + Vop

    """
    # compute exact ground state of HMI
    print "compute numerically exact HMI wavefunctions"
    from DFTB.Scattering.hydrogen_molecular_ion import DimerWavefunctions
    wfn = DimerWavefunctions(R,1.0,1.0, plot=False)

    en0, (Rfunc,Sfunc,Pfunc),wavefunction_exact = wfn.getBoundOrbital(0,0,'cos',0)
    print "HMI ground state energy = %e" % en0
    def phi_exact(x,y,z):
        return wavefunction_exact((x,y,z), None)
    """

    def psi_1s(x,y,z):
        #wavefunction of 1s hydrogen electron
        r = np.sqrt(x*x+y*y+z*z)
        psi = 1.0/np.sqrt(np.pi) * np.exp(-r)
        return psi

    def psi_sigma(x,y,z):
        # unnormalized LCAO wavefunction of H-H
        # sigma orbital is a linear combination of two
        # 1s orbitals of hydrogen
        psi = psi_1s(x,y,z-R/2.0) + psi_1s(x,y,z+R/2.0)
        return psi

    phi_exact = psi_sigma
    en0 = -1.0


    # evaluate exact wavefunction on the grid
    f_exact = grid.evaluate_function(phi_exact)

    # sum_i w_i f_i^2
    norm2 = grid.scalar_product(f_exact, f_exact)
    print "normalization norm^2 = %e" % norm2

    # sum_i w_i f_i (H.f)_i
    Hf = Hop.matvec(f_exact)
    en = grid.scalar_product(f_exact, Hf)
    print "energy expectation value <psi|H|psi> of exact ground state"
    print en

    # sum_i w_i f_i ((T+V).f)_i
    TpVf = -0.5 * Lop.matvec(f_exact) + Vop.matvec(f_exact)
    en = grid.scalar_product(f_exact, TpVf)
    print "energy expectation value <psi|T+V|psi> of exact ground state"
    print en

    print "H.f_exact"
    print Hf.tolist()
    print "(T+V).f_exact"
    print TpVf.tolist()

    err = la.norm(Hf - TpVf)
    print "|H.f - (T+V).f|= %e" % err
    
    # plot the exact eigenfunction
    import matplotlib.pyplot as plt
    plt.cla()
    plt.clf()
    # grid for plotting
    r = np.linspace(-15.0, 15.0, 100000)
    x = 0*r
    y = 0*r
    z = r
    
    # plot function belonging to lowest eigenvalue
    phi_exact_interp = grid.interpolation_func(f_exact)

    # check if repeated interpolation changes the function values
    f_interp = np.copy(f_exact)
    for i in range(0, 5):
        f_interp = grid.evaluate_function(grid.interpolation_func(f_interp))
    phi_interp = grid.interpolation_func(f_interp)

    # residual
    r_exact = Hop.matvec(f_exact) - en0 * f_exact
    residual_exact = grid.interpolation_func(r_exact)

    Tphi = Top.matvec(f_exact)
    VmEphi = Vop.matvec(f_exact) - en0*f_exact
    Tphi_func = grid.interpolation_func(Tphi)
    VmEphi_func = grid.interpolation_func(VmEphi)

    potential_interp = grid.interpolation_func(grid.evaluate_function(potential))
    
    plt.xlabel("z / bohr")
    plt.plot(r, phi_exact(x,y,z), label=r"HMI ground state (exact)")
    plt.plot(r, phi_exact_interp(x,y,z), label=r"HMI ground state (exact, interp.)")
    plt.plot(r, phi_interp(x,y,z), label=r"HMI ground state (exact, 5 x interp.)")
    plt.plot(r, residual_exact(x,y,z), label=r"residual $(H-E)\phi$")
    plt.plot(r, Tphi_func(x,y,z), label=r"$T \phi$")
    plt.plot(r, VmEphi_func(x,y,z), label=r"$(V-E) \phi$")
    # potential energy
    plt.plot(r, potential_interp(x,y,z), label=r"$V$ (interpolated)")
    plt.plot(r, potential(x,y,z), ls="--", label=r"$V$")
    plt.legend()

    plt.show()

    
        
if __name__ == "__main__":
    #test_angular_laplacian()
    #test_radial_laplacian()
    #test_radial_laplacian_symmetry()
    #test_laplacian()
    #test_multicenter_laplacian()
    #test_possion_solver_sparse() # not working
    #test_integration()
    #test_hydrogen_spherical_grid()
    #test_interpolation()
    
    #test_linear_operator_hydrogen()
    test_multicenter_hmi()   # not working
