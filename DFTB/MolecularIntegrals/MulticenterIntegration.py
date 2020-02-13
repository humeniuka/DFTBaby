#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Becke's scheme for multicenter integration, solution of Poisson's equation and 
inhomogeneous Schroedinger equation

References
----------
[1] A.Becke, "A multicenter numerical integration scheme for polyatomic molecules",
    J.Chem.Phys. 88, 2547 (1988)
[2] A.Becke, R.Dickson, "Numerical solution of Poisson's equation in polyatomic molecules",
    J.Chem.Phys. 89, 2993 (1988)
[3] A.Becke, R.Dickson, "Numerical solution of Schroedinger's equation in polyatomic molecules",
    J.Chem.Phys. 92, 3610 (1990)

Some useful information is also contained in
[4] T.Shiozaki, S.Hirata, "Grid-based numerical Hartree-Fock solutions of polyatomic molecules",
    Phys.Rev. A 76, 040503(R) (2007)

"""
import numpy as np
import numpy.linalg as la
from scipy import interpolate
import mpmath

from DFTB.MolecularIntegrals.LebedevQuadrature import get_lebedev_grid, outerN, spherical_harmonics_it
from DFTB.MolecularIntegrals.lebedev_quad_points import LebedevGridPoints, Lebedev_Npts, Lebedev_Lmax
from DFTB.MolecularIntegrals.SphericalCoords import cartesian2spherical
# for testing
from DFTB.MolecularIntegrals.integrals import norm, coulomb_repulsion

from DFTB.AtomicData import bohr_to_angs, atom_names, slater_radii


def number_of_radial_points(Z):
    """
    select the number of radial grid points for the subintegral
    around an atom with atomic number Z
    """
    # Hydrogen atom receives an initial quota of 20 points
    Nr = 20
    # Each shell receives an additional 5 points
    if Z >= 2:
        Nr += 5
    if Z >= 11:
        Nr += 5
    if Z >= 19:
        Nr += 5
    if Z >= 37:
        Nr += 5
    if Z >= 55:
        Nr += 5
    if Z >= 87:
        Nr += 5
        
    return Nr

def select_angular_grid(lebedev_order):
    # find closest Lebedev grid of requested order
    n_lebedev = abs(np.array(Lebedev_Lmax) - lebedev_order).argmin()
    Lmax = Lebedev_Lmax[n_lebedev]
    if lebedev_order != Lmax:
        print "No grid for order %s, using grid which integrates up to Lmax = %s exactly instead." \
            % (lebedev_order, Lmax)
    # angular grid
    th,ph,angular_weights = np.array(LebedevGridPoints[Lebedev_Npts[n_lebedev]]).transpose()
    
    return Lmax, (th,ph,angular_weights)

def print_grid_summary(atomlist,
                       lebedev_order=23, radial_grid_factor=1):
    """
    print a table with number of radial and angular grid points
    for each atom. The total number of grid points in all overlapping
    grids is also shown.

    Parameters
    ----------
    atomlist     :  list of tuples (Zat,(x,y,z)) with
                    molecular geometry

    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    grid_sizes         : list of tuples (Nr, Nang)
                         giving the number of radial and angular grid points for
                         each atom
    
    """
    print " "
    print "   Size of grids"
    print " "
    print " Atom      #    radial points    angular points     radial x angular"
    print " -------------------------------------------------------------------"
    Ntot = 0
    grid_sizes = []
    for i,(Z,pos) in enumerate(atomlist):
        # radial grid
        Nr = number_of_radial_points(Z)
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        # angular grid
        Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
        Nang = th.size
        Ntot += Nr*Nang
        
        grid_sizes.append( (Nr,Nang) )
        
        print " %2s      %3.1d      %7.1d           %7.1d           %10.1d" % (atom_names[Z-1], i+1, Nr, Nang, Nr*Nang)

    print " -------------------------------------------------------------------"
    print " Total                                                %10.1d" % Ntot
    print ""
    
    return grid_sizes
    
def multicenter_integration(f, atomic_coordinates, atomic_numbers, lebedev_order=23, radial_grid_factor=1):
    """
    compute the integral

             / 
         I = | f(x,y,z) dV
             / 

    numerically on a multicenter spherical grid using Becke's scheme 
   
    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the function at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    I       : value of the integral
    """
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
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

    integral = 0.0
    # atom-centered subintegral
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        xr = np.cos(k/(Nr+1.0) * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
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
    
        # weight function
        wr = P[:,I]/Ptot

        # evaluate function on the grid
        fI = wr * f(x,y,z)

        sub_integral = np.sum( weights * fI )
        #print "I= %d    sub_integral= %s" % (I, sub_integral)
        integral += sub_integral

    return integral

def multicenter_poisson(f, atomic_coordinates, atomic_numbers, lebedev_order=23, radial_grid_factor=1):
    """
    solve the Poisson equation
      __2
      \/  V(r) = -4 pi f(r)

    numerically on a multicenter spherical grid

    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the charge distribution at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    V       : callable, V(x,y,z) evaluates the electrostatic potential generated
              by the charge distribution f.
    """    
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
        f = mu
        for ik in range(0, k):
            f = 1.5 * f -0.5 * f**3
        return 0.5*(1-f)
    
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

    radial_potentials = []
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
        # from variable transformation
        g = 2 * rm**3 * np.sqrt(((1+xr)/(1-xr)**3)**3)
        radial_weights *= g
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # cartesian coordinates of grid
        x = (outerN(r, sc) + atomic_coordinates[0,I])
        y = (outerN(r, ss) + atomic_coordinates[1,I])
        z = (outerN(r, c ) + atomic_coordinates[2,I])
        weights = outerN(radial_weights, 4.0*np.pi * angular_weights)
        #
        Npts = Nr*Nang
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # evaluate function on the grid
        fI = wr * f(x,y,z)

        # solve the Poisson equation
        # __2  (I)
        # \/  V    = - 4 pi f_I(r)
        #

        # total charge qI
        qI = np.sum( weights * fI)
        
        # expand charge distribution f_I(r) into spherical harmonics
        #
        #  f_I(x,y,z) = sum_l sum_{m=-l}^{l}  f_lm(r) Y_{l,m}(th,ph)
        #
        # and solve Poisson equation for each (l,m) component
        radial_potentials.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            fI_lm = 4.0*np.pi * np.sum(fI*wYlm, axis=-1)

            # boundary conditions for u(z)
            #  z=0 <-> x=1  <-> r=+infinity
            #  z=1 <-> x=-1 <-> r=0
            if (l,m) == (0,0):
                u0 = np.sqrt(4.0*np.pi) * qI
            else:
                u0 = 0.0
            u1 = 0.0
            #
            omx = 1-xr
            opx = 1+xr
            #   l*(l+1)
            # - ------- u = c0 u
            #     r^2
            c0 = -l*(l+1)/rm**2 * (omx/opx)**2
            # The transformation of partial derivatives from
            # r-coordinates to z-coordinates is accomplished as
            #  d^2 u      d u      d^2 u
            #  ----- = c1 --- + c2 -----
            #  d r^2      d z      d z^2
            c1 = 1.0/(4.0*np.pi*rm**2) * omx**(2.5) * (1+opx) / opx**(1.5)
            c2 = 1.0/(4.0*np.pi**2*rm**2) * omx**3 / opx
            # source term
            #   -4 pi r fI_lm  = -4 pi rm (1+x)/(1-x) fI_lm
            source = -4.0*np.pi * rm * opx / omx * fI_lm

            u_lm = solve_radial_dgl(c0,c1,c2,source, u0, u1)
            # r*V(r) = u(r)
            VI_lm = u_lm/r
            
            # z-array is sorted in assending order,
            #   0 < z[0] < z[1] < ... < z[-1] < 1
            # while associated r-array is sorted in descending order
            #   infinity > r[0] > r[1] > ... > r[-1] > 0

            spline_lm_real = interpolate.splrep(zr, VI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, VI_lm.imag, s=0)
            radial_potentials[-1][(l,m)] = spline_lm_real, spline_lm_imag

            # Up to (including) l = (Lmax-1)/2 the integration on
            # the Lebedev grid is exact. If higher spherical harmonics
            # were to be included, interpolated function can blow
            # up after several iterations.
            if m == -(Lmax-1)/2:
                break

    def electrostatic_potential_func(x,y,z):
        """
        function for evaluating the electrostatic potential V(x,y,z)
        """
        V = 0j*x
        # build the total potential from the solutions to the subproblems
        #  V = sum_n  V^(n)
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:
                
                spline_lm_real, spline_lm_imag = radial_potentials[I][(l,m)]
                # interpolate
                VI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)
                V += VI_lm*Ylm

                if m == -(Lmax-1)/2:
                    break

        return V.real
    
    return electrostatic_potential_func

    
def solve_radial_dgl(c0,c1,c2, source, u0, u1):
    """
    solve the differential equation on the interval [0,1]

                        du         d^2 u
     c0(z) u(z) + c1(z) -- + c2(z) -----  = source(z)
                        dz         dz^2

    for u(z) subject to the boundary conditions

    u(z=0) = u0
    u(z=1) = u1

    Central finite difference formula for 1st and 2nd derivatives are taken from Bickley [1].


    Parameters
    ----------
    c0,c1,c2, source  :   values of c0(z), c1(z), c2(z), source(z) on an equidistant
                          grid z_i = i/(N+1)  i=1,...,N
    u0, u1            :   scalar value, boundary conditions u(z=0), u(z=1)

    Returns
    -------
    u                 :   solution u(z) on the equidistant grid

    References
    ----------
    [1] W. Bickley, "Formulae for Numerical Differentiation",
        The Mathematical Gazette, vol. 25, no. 263, pp. 19-27 (1941)

    """
    # convert the differential equation into an algebraic equation
    #   A.u = b
    N = len(c0)
    h = 1.0/(N+1)    # separation between equidistant points

    # operators d/dz and d^2/dz^2
    D1 = np.zeros((N,N))
    D2 = np.zeros((N,N))
    # terms from boundary conditions
    b1 = np.zeros(N)
    b2 = np.zeros(N)
    # non-centered five-point formulae for i=0
    D1[0,0:4] = np.array([                    -20.0, +36.0, -12.0, 2.0])/(24.0*h)
    b1[0]     =           -6.0/(24.0*h) * u0
    D2[0,0:4] = np.array([                    -20.0,  +6.0,  +4.0, -1.0])/(12.0*h**2)
    b2[0]     =           +11.0/(12.0*h**2) * u0
    # non-centered six-point formulae for i=1
    D1[1,0:5] = np.array([                    -60.0, -40.0, +120.0, -30.0, +4.0])/(120.0*h)
    b1[1]     =            6.0/(120.0*h) * u0
    D2[1,0:5] = np.array([                    +80.0,-150.0,  +80.0, -5.0,  0.0])/(60.0*h**2)
    b2[1]     =           -5.0/(60.0*h**2) * u0
    # centered seven-point formulae for i=2
    D1[2,0:6] = np.array([                     +108.0, -540.0, 0.0, 540.0, -108.0, 12.0])/(720.0*h)
    b1[2]     =           -12.0/(720.0*h) * u0
    D2[2,0:6] = np.array([                      -54.0, +540.0, -980.0, +540.0, -54.0, +4.0])/(360.0*h**2)
    b2[2]     =             4.0/(360.0*h**2) * u0
    # centered seven-point formulae for i=3,...,N-4
    for i in range(3, N-3):
        D1[i,i-3:i+4] = np.array([-12.0, 108.0, -540.0, 0.0, +540.0, -108.0, +12.0])/(720.0*h)
        D2[i,i-3:i+4] = np.array([  4.0, -54.0, +540.0, -980.0, +540.0, -54.0, +4.0])/(360.0*h**2)
    # centered seven-point formulae for i=N-3
    D1[N-3,N-6:] = np.array([-12.0, +108.0, -540.0, 0.0, 540.0, -108.0              ])/(720.0*h)
    b1[N-3]      =                                                       +12.0/(720.0*h) * u1
    D2[N-3,N-6:] = np.array([+4.0, -54.0, 540.0, -980.0, +540.0, -54.0              ])/(360.0*h**2)
    b2[N-3]      =                                                        +4.0/(360.0*h**2) * u1
    # non-centered six-point formulae for i=N-2
    D1[N-2,N-5:] = np.array([-4.0, +30.0, -120.0, +40.0, +60.0                  ])/(120.0*h)
    b1[N-2]      =                                                -6.0/(120.0*h) * u1
    D2[N-2,N-5:] = np.array([0.0, -5.0, +80.0, -150.0, +80.0                    ])/(60.0*h**2)
    b2[N-2]      =                                                -5.0/(60.0*h**2) * u1
    # non-centered five-point formulae for i=N-1
    D1[N-1,N-4:] = np.array([-2.0, +12.0, -36.0, +20.0            ])/(24.0*h)
    b1[N-1]      =                                      +6.0/(24.0*h) * u1
    D2[N-1,N-4:] = np.array([-1.0, +4.0, +6.0, -20.0              ])/(12.0*h**2)
    b2[N-1]      =                                      +11.0/(12.0*h**2) * u1

    # The differential equation is transformed into the following linear equation
    #
    #  sum_j [ c0(i) delta_ij + c1(i) D1(i,j) + c2(i) D2(i,j) ] u(j) = source(i) - c1(i)*b1(i) - c2(i)*b2(i)
    #
    
    # build matrix A on the left hand side of the equation
    A = np.zeros((N,N))
    for i in range(0, N):
        A[i,i] = c0[i]
        A[i,:] += c1[i]*D1[i,:] + c2[i]*D2[i,:]

    # right hand side
    rhs = source - c1*b1 - c2*b2
    # solve matrix equation
    u = la.solve(A, rhs)
    
    return u


def multicenter_laplacian_OLD(f, atomic_coordinates, atomic_numbers, lebedev_order=23, radial_grid_factor=1):
    """
    compute the action of the Laplace operator on a function f
      __2
      \/  f(r)

    numerically on a multicenter spherical grid.

    The algorithm is very similar to but simpler than `multicenter_poisson()`: 
    The function f is decomposed using the weight functions wI(r)

         f(r) = sum_I fI(r)             with  fI(r) = wI(r) f(r)

    Each part fI(r) is expanded into spherical harmonics, fI = sum_lm fI_lm(r) Y_lm(th,ph),
    and the Laplace operator is applied in spherical harmonics,
         __2       d^2      L
         \/  = 1/r --- r - --- .
                   dr^2    r^2
    At the end the Laplacians of the parts are summed
         __2           __2
         \/  f = sum_I \/ fI(r)

    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the function at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    lapf               : callable, lapf(x,y,z) evaluates the Laplacian of f
    """    
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
        f = mu
        for ik in range(0, k):
            f = 1.5 * f -0.5 * f**3
        return 0.5*(1-f)
    
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

    radial_laplacians = []
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
        # from variable transformation
        g = 2 * rm**3 * np.sqrt(((1+xr)/(1-xr)**3)**3)
        radial_weights *= g
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # cartesian coordinates of grid
        x = (outerN(r, sc) + atomic_coordinates[0,I])
        y = (outerN(r, ss) + atomic_coordinates[1,I])
        z = (outerN(r, c ) + atomic_coordinates[2,I])
        weights = outerN(radial_weights, 4.0*np.pi * angular_weights)
        #
        Npts = Nr*Nang
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # evaluate function on the grid
        fI = wr * f(x,y,z)

        # compute Laplacian
        #    (I)    __2  (I)
        # lap    =  \/  f    
        #

        # expand function f_I(r) into spherical harmonics
        #
        #  f_I(x,y,z) = sum_l sum_{m=-l}^{l}  f_lm(r) Y_{l,m}(th,ph)
        #
        # and apply Laplacian to each (l,m) component
        radial_laplacians.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            fI_lm = 4.0*np.pi * np.sum(fI*wYlm, axis=-1)

            omx = 1-xr
            opx = 1+xr
            # coefficients for Laplacian operator after coordinate transformation
            c0 = -l*(l+1)/rm**2 * (omx/opx)**2
            c1 = 1.0/(4.0*np.pi*rm**2) * omx**(2.5) * (1.0+opx) / opx**(1.5)
            c2 = 1.0/(4.0*np.pi**2*rm**2) * omx**3 / opx

            #  
            # f_lm(r) = 1/r u_lm
            #
            uI_lm = r*fI_lm

            # apply Laplacian operator
            #            d^2     l(l+1)
            # L u_lm = ( ----  - ------ ) u_lm
            #            dr^2     r^2
            LuI_lm = radial_difop(c0,c1,c2, uI_lm)
            
            # __2  (I)              1    d^2     l(l+1)
            # \/  f    = sum_{l,m} --- ( ----  - ------ ) u_lm(r) Y_lm(th,ph)
            #                       r    dr^2     r^2
            #                         
            #          = sum_{l,m} lapI_lm(r) Y_lm(th,ph)
            #
            #lapI_lm = LuI_lm/r
            
            # z-array is sorted in assending order,
            #   0 < z[0] < z[1] < ... < z[-1] < 1
            # while associated r-array is sorted in descending order
            #   infinity > r[0] > r[1] > ... > r[-1] > 0

            # It seems better to spline LuI_lm instead of LuI_lm/r,
            # because the first function is smoother and has no
            # singularity at r=0.
            spline_lm_real = interpolate.splrep(zr, LuI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, LuI_lm.imag, s=0)
            radial_laplacians[-1][(l,m)] = spline_lm_real, spline_lm_imag

            if m == -(Lmax-1)/2:
                break

    def laplacian_func(x,y,z):
        """
        function for evaluating the Laplacian of f
             __2
             \/  f
        """
        lap = 0j*x
        # build the total Laplacian from the Laplacians of the parts
        #  __2           __2
        #  \/ f = sum_n  \/  f^(n)
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:
                
                spline_lm_real, spline_lm_imag = radial_laplacians[I][(l,m)]
                # interpolate
                LuI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                lapI_lm = LuI_lm/rI
                lap += lapI_lm*Ylm

                if m == -(Lmax-1)/2:
                    break

        return lap.real
    
    return laplacian_func

def multicenter_laplacian(f, atomic_coordinates, atomic_numbers,
                          cusps_separate=True,
                          lebedev_order=23, radial_grid_factor=1):
    """
    compute the action of the Laplace operator on a function f
      __2
      \/  f(r)

    numerically on a multicenter spherical grid.

    The algorithm is very similar to but simpler than `multicenter_poisson()`: 
    The function f is decomposed using the weight functions wI(r)

         f(r) = sum_I fI(r)             with  fI(r) = wI(r) f(r)

    Each part fI(r) is expanded into spherical harmonics, fI = sum_lm fI_lm(r) Y_lm(th,ph),
    and the Laplace operator is applied in spherical harmonics,
         __2       d^2      L
         \/  = 1/r --- r - --- .
                   dr^2    r^2
    At the end the Laplacians of the parts are summed
         __2           __2
         \/  f = sum_I \/ fI(r)

    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the function at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    cusps_separate     : if True, Laplacians of `f` around the atoms are calculated
                         separately to avoid numerical artifacts that arise from the cusps
                         if `f` represents a wavefunction
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    lapf               : callable, lapf(x,y,z) evaluates the Laplacian of f
    """
    Nat = atomic_coordinates.shape[1]
    atomic_names = [atom_names[Z-1] for Z in atomic_numbers]
    ###
    # Calculating the Laplacian is complicated by the fact that
    # electronic wavefunctions have cusps at the nuclear positions.
    # We split the wavefunction into a sum of spherically symmetric
    # wavefunctions g_i around each atom with cusps and the remainder.
    #             
    #   f = sum  g (|r-Ri|)  +  [ f  -  sum  g (|r-Rj|) ]
    #          i  i                j       j  j
    #
    # The Laplacian is computed separately for the g_i's
    # using the fact the for a spherically symmetric function
    #
    #  __2             d^2           
    #  \/  g (r) = 1/r ---- (r g (r)) 
    #       i          dr^2     i     
    # 
    # The remainder should be relatively smooth everywhere and pose
    # no problem to numerical differentiation.
    #
    
    # 1) First we calculate the spherical averages of f around each atom
    #    and spline the resulting functions g_i(z) and their Laplacian
    #    using z-coordinates.
    # list of splines g_i
    radial_avg = []
    # list of splines of the Laplacian of the spherical average 
    radial_lap_avg = []
    for I in  range(0, Nat):
        # center around which the average is performed
        Zat = atomic_numbers[I]
        pos = atomic_coordinates[:,I]
        atom = (Zat, pos)
        #
        gI_func = spherical_average_func(atom, f,
                                         lebedev_order=lebedev_order,
                                         radial_grid_factor=radial_grid_factor)

        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # variable transformation
        omx = 1-xr
        opx = 1+xr

        gI = gI_func(r)
        
        # coefficients for Laplacian operator after coordinate transformation
        c0 = 0.0*omx  # -l*(l+1)/rm**2 * (omx/opx)**2  is zero since l=0 
        c1 = 1.0/(4.0*np.pi*rm**2) * omx**(2.5) * (1+opx) / opx**(1.5)
        c2 = 1.0/(4.0*np.pi**2*rm**2) * omx**3 / opx

        # apply Laplacian operator to spherically symmetric function
        #         1  d^2 
        # L gI = --- ---- ( r gI(r) )
        #         r  dr^2

        # Is better to interpolate r*L(gI) instead of L(gI)
        # in order to avoid the singularity at r=0
        #  r L gI = d^2/dr^2 (r gI)
        rLgI = radial_difop(c0,c1,c2, r*gI)
        
        # interpolate gI(z) and LgI(z) on a grid
        spline_avg = interpolate.splrep(zr, gI, s=0)
        radial_avg.append(spline_avg)

        spline_lap_avg = interpolate.splrep(zr, rLgI, s=0)
        radial_lap_avg.append( spline_lap_avg )

    """
    ### DEBUG
    import matplotlib.pyplot as plt
    zr = np.linspace(0.0, 10.0, 1000)
    for I in range(0, Nat):
        plt.plot(zr, interpolate.splev(zr, radial_avg[I], der=0, ext=0),     label="around atom %d" % (I+1))
        plt.plot(zr, interpolate.splev(zr, radial_lap_avg[I], der=0, ext=0), label="Laplacian, around atom %d" % (I+1))
    plt.legend()
    plt.show()
    ###
    """

    ####
    
    # 2) Define weight functions for fuzzy Voronoi decomposition
    #
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
        f = mu
        for ik in range(0, k):
            f = 1.5 * f -0.5 * f**3
        return 0.5*(1-f)
    
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

    # 3) Evaluate function f on the multicenter grid and perform a spherical
    #    wave decomposition

    # 
    radial_laplacians = []
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
        # from variable transformation
        g = 2 * rm**3 * np.sqrt(((1+xr)/(1-xr)**3)**3)
        radial_weights *= g
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # cartesian coordinates of grid
        x = (outerN(r, sc) + atomic_coordinates[0,I])
        y = (outerN(r, ss) + atomic_coordinates[1,I])
        z = (outerN(r, c ) + atomic_coordinates[2,I])
        weights = outerN(radial_weights, 4.0*np.pi * angular_weights)
        #
        Npts = Nr*Nang
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # evaluate function on the grid
        fI = wr * f(x,y,z)

        if cusps_separate:
            ####
            # subtract spherical averages around each atom
            #
            #  f_i   = w (r) [ f(r) - sum  g (|r-Rj|) ]
            #           i                j  j
            for j in range(0, Nat):
                # rJ = |r-Rj]
                rJ = dist[:,:,j]
                # transform to z-coordinates
                rmJ = 0.5*slater_radii[atomic_names[j]]
                xrJ = (rJ-rmJ)/(rJ+rmJ)
                zrJ = np.arccos(xrJ) / np.pi
                
                gJ_spline = radial_avg[j]
                fI -= wr * interpolate.splev(zrJ, gJ_spline, der=0, ext=0)
            ####
        
        # compute Laplacian
        #    (I)    __2  (I)
        # lap    =  \/  f    
        #

        # expand function f_I(r) into spherical harmonics
        #
        #  f_I(x,y,z) = sum_l sum_{m=-l}^{l}  f_lm(r) Y_{l,m}(th,ph)
        #
        # and apply Laplacian to each (l,m) component
        radial_laplacians.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            fI_lm = 4.0*np.pi * np.sum(fI*wYlm, axis=-1)

            omx = 1-xr
            opx = 1+xr
            # coefficients for Laplacian operator after coordinate transformation
            c0 = -l*(l+1)/rm**2 * (omx/opx)**2
            c1 = 1.0/(4.0*np.pi*rm**2) * omx**(2.5) * (1+opx) / opx**(1.5)
            c2 = 1.0/(4.0*np.pi**2*rm**2) * omx**3 / opx

            #  
            # f_lm(r) = 1/r u_lm
            #
            uI_lm = r*fI_lm

            # apply Laplacian operator
            #            d^2     l(l+1)
            # L u_lm = ( ----  - ------ ) u_lm
            #            dr^2     r^2
            LuI_lm = radial_difop(c0,c1,c2, uI_lm)

            # __2  (I)              1    d^2     l(l+1)
            # \/  f    = sum_{l,m} --- ( ----  - ------ ) u_lm(r) Y_lm(th,ph)
            #                       r    dr^2     r^2
            #                         
            #          = sum_{l,m} lapI_lm(r) Y_lm(th,ph)
            #
            #lapI_lm = LuI_lm/r
            
            # z-array is sorted in assending order,
            #   0 < z[0] < z[1] < ... < z[-1] < 1
            # while associated r-array is sorted in descending order
            #   infinity > r[0] > r[1] > ... > r[-1] > 0

            # It seems better to spline LuI_lm instead of LuI_lm/r,
            # because the first function is smoother and has no
            # singularity at r=0.

            spline_lm_real = interpolate.splrep(zr, LuI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, LuI_lm.imag, s=0)
            radial_laplacians[-1][(l,m)] = spline_lm_real, spline_lm_imag

            if m == -(Lmax-1)/2:
                break

    #                                   __2
    # 4) Define function for evaluating \/ f
    def laplacian_func(x,y,z):
        """
        function for evaluating the Laplacian of f
             __2
             \/  f
        """
        lap = 0j*x
        # build the total Laplacian from the Laplacians of the parts
        #  __2           __2
        #  \/ f = sum_I  \/  f^(I)
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:
                
                spline_lm_real, spline_lm_imag = radial_laplacians[I][(l,m)]
                # interpolate
                LuI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)
                lapI_lm = LuI_lm/rI
                lap += lapI_lm*Ylm
                
                if cusps_separate:
                    #####
                    if l == 0 and m == 0:
                        # add back the Laplacian due to spherical
                        # averaged part g_i that has been subtracted previously:
                        #   __2    __2                           __2
                        #   \/ f = \/ [ f - sum_i g_i ]  + sum_i \/ g_i
                        #   
                        rLgI = interpolate.splev(zr, radial_lap_avg[I], der=0, ext=0)
                        lap += rLgI/rI
                    #####
                    
                if m == -(Lmax-1)/2:
                    break

        return lap.real
    
    return laplacian_func

def radial_difop(c0,c1,c2, u):
    """
    apply the differential operator L defined by the coefficients c0,c1 and c2
    to a function u(z) defined on a grid covering the range [0,1]

                                   du         d^2 u
       L u(z) = c0(z) u(z) + c1(z) -- + c2(z) -----
                                   dz         dz^2

    Central finite difference formula for 1st and 2nd derivatives are taken from Bickley [1].


    Parameters
    ----------
    c0,c1,c2          :   values of c0(z), c1(z), c2(z) on an equidistant
                          grid z_i = i/(N+1)  i=1,...,N
                          defining the operator L
    u                 :   values of u(z) on the grid

    Returns
    -------
    Lu                :   result of applying the operator L to u(z)

    References
    ----------
    [1] W. Bickley, "Formulae for Numerical Differentiation",
        The Mathematical Gazette, vol. 25, no. 263, pp. 19-27 (1941)

    """
    N = len(u)      # number of grid points
    h = 1.0/(N+1)    # separation between equidistant points
    
    # operators d/dz and d^2/dz^2
    D1 = np.zeros((N,N))
    D2 = np.zeros((N,N))
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
    # centered seven-point formulae for i=3,...,N-4
    for i in range(3, N-3):
        D1[i,i-3:i+4] = np.array([-12.0, +108.0, -540.0,    0.0, +540.0, -108.0, +12.0])/(720.0*h)
        D2[i,i-3:i+4] = np.array([ +4.0,  -54.0, +540.0, -980.0, +540.0,  -54.0,  +4.0])/(360.0*h**2)
    # non-centered seven-point formulae for i=N-3
    # D^1 u_{N-3}   ~   D^1 u_4
    D1[N-3,N-7:] = np.array([+12.0, -96.0, +360.0, -960.0, +420.0, +288.0, -24.0 ])/(720.0*h)
    # D^2 u_{N-3}   ~   D^2 u_4
    D2[N-3,N-7:] = np.array([ +4.0, -24.0,  +30.0, +400.0, -840.0, +456.0, -26.0 ])/(360.0*h**2)
    # non-centered seven-point formulae for i=N-2
    # D^1 u_{N-2}   ~   D^1 u_5
    D1[N-2,N-7:] = np.array([ -24.0, +180.0, -600.0, +1200.0, -1800.0, +924.0, +120.0])/(720.0*h)
    # D^2 u_{N-2}   ~   D^2 u_5
    D2[N-2,N-7:] = np.array([ -26.0, +186.0, -570.0,  +940.0,  -510.0, -294.0, +274.0])/(360.0*h**2)
    # non-centered seven-point formulae for i=N-1
    # D^1 u_{N-1}   ~   D^1 u_6
    D1[N-1,N-7:] = np.array([ +120.0,  -864.0, +2700.0,  -4800.0,  +5400.0, -4320.0, +1764.0])/(720.0*h)
    # D^2 u_{N-1}   ~   D^2 u_6
    D2[N-1,N-7:] = np.array([ +274.0, -1944.0, +5940.0, -10160.0, +10530.0, -6264.0, +1624.0])/(360.0*h**2)

    # finite difference formula converts differential operators into matrices
    Lu = c0*u + c1*np.dot(D1,u) + c2*np.dot(D2,u)

    return Lu


def multicenter_laplacian_SMOOTH(f, atomic_coordinates, atomic_numbers, lebedev_order=23, radial_grid_factor=1):
    """
    compute the action of the Laplace operator on a function f
      __2
      \/  f(r)

    numerically on a multicenter spherical grid.

    The algorithm is very similar to but simpler than `multicenter_poisson()`: 
    The function f is decomposed using the weight functions wI(r)

         f(r) = sum_I fI(r)             with  fI(r) = wI(r) f(r)

    Each part fI(r) is expanded into spherical harmonics, fI = sum_lm fI_lm(r) Y_lm(th,ph),
    and the Laplace operator is applied in spherical harmonics,
         __2       d^2      L
         \/  = 1/r --- r - --- .
                   dr^2    r^2
    At the end the Laplacians of the parts are summed
         __2           __2
         \/  f = sum_I \/ fI(r)

    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the function at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    lapf               : callable, lapf(x,y,z) evaluates the Laplacian of f
    """    
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
        f = mu
        for ik in range(0, k):
            f = 1.5 * f -0.5 * f**3
        return 0.5*(1-f)
    
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

    radial_laplacians = []
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
#        # The end point z=+1 is included!
#        zr = k/(Nr+0.0)
        zr = k/(Nr+1.0)        
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(zr * np.pi)**2
        # from variable transformation
        g = 2 * rm**3 * np.sqrt(((1+xr)/(1-xr)**3)**3)
        radial_weights *= g
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # cartesian coordinates of grid
        x = (outerN(r, sc) + atomic_coordinates[0,I])
        y = (outerN(r, ss) + atomic_coordinates[1,I])
        z = (outerN(r, c ) + atomic_coordinates[2,I])
        weights = outerN(radial_weights, 4.0*np.pi * angular_weights)
        #
        Npts = Nr*Nang
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # evaluate function on the grid
        fI = wr * f(x,y,z)

        # compute Laplacian
        #    (I)    __2  (I)
        # lap    =  \/  f    
        #

        # expand function f_I(r) into spherical harmonics
        #
        #  f_I(x,y,z) = sum_l sum_{m=-l}^{l}  f_lm(r) Y_{l,m}(th,ph)
        #
        # and apply Laplacian to each (l,m) component
        radial_laplacians.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            fI_lm = 4.0*np.pi * np.sum(fI*wYlm, axis=-1)
            # f_lm(r) = 1/r u_lm
            uI_lm = r*fI_lm
            #
            # __2  (I)              1    d^2     l(l+1)
            # \/  f    = sum_{l,m} --- ( ----  - ------ ) u_lm(r) Y_lm(th,ph)
            #                       r    dr^2     r^2
            #                         
            #          = sum_{l,m} lapI_lm(r) Y_lm(th,ph)
            #

            # apply Laplacian operator
            #
            #                       d^2     l(l+1)
            #   1/r  L u(r) = 1/r ( ----  - ------ ) u(r)
            #                       dr^2     r^2
            # using splines
            #
            # The radial grid points have to be provided in increasing order.            
            # z-array is sorted in assending order,
            #   0 < z[0] < z[1] < ... < z[-1] <= 1
            # while associated r-array is sorted in descending order
            #   infinity > r[0] > r[1] > ... > r[-1] >= 0

            omx = 1-xr
            opx = 1+xr
            #   l*(l+1)
            # - ------- u = c0 u
            #     r^2
            c0 = -l*(l+1)/rm**2 * (omx/opx)**2
            # The transformation of partial derivatives from
            # r-coordinates to z-coordinates is accomplished as
            #  d^2 u      d u      d^2 u
            #  ----- = c1 --- + c2 -----
            #  d r^2      d z      d z^2
            c1 = 1.0/(4.0*np.pi*rm**2) * omx**(2.5) * (1+opx) / opx**(1.5)
            c2 = 1.0/(4.0*np.pi**2*rm**2) * omx**3 / opx

            LuI_lm = radial_difop(c0, c1, c2, uI_lm)

            eps = 1.0e-3
            lapI_lm_real_func = spline_radial_wavefunction(r[r > eps][::-1], LuI_lm[r > eps][::-1].real)
            lapI_lm_imag_func = spline_radial_wavefunction(r[r > eps][::-1], LuI_lm[r > eps][::-1].imag)
            radial_laplacians[-1][(l,m)] = (lapI_lm_real_func, lapI_lm_imag_func)
            """
            lapI_lm = LuI_lm/r

            ###
            import matplotlib.pyplot as plt
            plt.cla()
            plt.clf()
            plt.plot(zr, lapI_lm.real, label="lapI_lm")
            plt.legend()
            plt.show()
            
#            lapI_lm = radial_laplacian_spline(r[::-1], uI_lm[::-1], l)[::-1]
                
            spline_lm_real = interpolate.splrep(zr, lapI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, lapI_lm.imag, s=0)

            radial_laplacians[-1][(l,m)] = spline_lm_real, spline_lm_imag
            """
            """
            if l == 0 and m == 0:
                analyze_cusp(r[::-1], uI_lm[::-1], LuI_lm[::-1])
            
            # __2  (I)              1    d^2     l(l+1)
            # \/  f    = sum_{l,m} --- ( ----  - ------ ) u_lm(r) Y_lm(th,ph)
            #                       r    dr^2     r^2
            #                         
            #          = sum_{l,m} 1/r LuI_lm(r) Y_lm(th,ph)
            #
            
            # z-array is sorted in assending order,
            #   0 < z[0] < z[1] < ... < z[-1] <= 1
            # while associated r-array is sorted in descending order
            #   infinity > r[0] > r[1] > ... > r[-1] >= 0

            # It seems better to spline LuI_lm instead of LuI_lm/r,
            # because the first function is smoother and has no
            # singularity at r=0.

            spline_lm_real = interpolate.splrep(zr, LuI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, LuI_lm.imag, s=0)

            radial_laplacians[-1][(l,m)] = spline_lm_real, spline_lm_imag
            """

            if m == -(Lmax-1)/2:
                break

    def laplacian_func(x,y,z):
        """
        function for evaluating the Laplacian of f
             __2
             \/  f
        """
        lap = 0j*x
        # build the total Laplacian from the Laplacians of the parts
        #  __2           __2
        #  \/ f = sum_n  \/  f^(n)
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:

                """
                spline_lm_real, spline_lm_imag = radial_laplacians[I][(l,m)]

                # interpolate
                lapI_lm =       interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                         + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)
                """
                lapI_lm_real_func, lapI_lm_imag_func = radial_laplacians[I][(l,m)]
                lapI_lm = lapI_lm_real_func(rI) + 1.0j*lapI_lm_imag_func(rI)
                
                lap += lapI_lm*Ylm

                if m == -(Lmax-1)/2:
                    break

        return lap.real
    
    return laplacian_func


def radial_laplacian_spline(r, u, l):
    """
     apply Laplacian operator

                   1    d^2     l(l+1)
     1/r L u(r) = --- ( ----  - ------ ) u(r)
                   r    dr^2     r^2

    The function u(r) is first splined, then the second derivative
    is computed from the spline.

    u(r)/r can only take a finite value for r --> 0, if the Taylor
    series of u(r) around r=0 has the form

                        1+i
      u(r) = sum    C  r 
                i=0  i

    Parameters
    ----------
    r           : 1d numpy array with grid points in decreasing order,
                  i.e. r[0] < r[1] < ... r[-1]
    u           : 1d numpy array with function values of u(r) at the
                  grid points
    l           : integer >= 0, angular momentum quantum number

    Returns
    -------
    lap         : 1d numpy array with value of `1/r L u(r)` at the grid points
    """
    # spline u(r)
    spl_re = interpolate.splrep(r, u.real, s=0, k=5)
    spl_im = interpolate.splrep(r, u.imag, s=0, k=5)
    # compute d^2u/dr^2 using spline
    d2udr2 =        interpolate.splev(r, spl_re, der=2, ext=0) \
             + 1.0j*interpolate.splev(r, spl_im, der=2, ext=0)
    # compute d^3u/dr^3(r=0) using spline
    d3udr3_0 =  interpolate.splev([0.0], spl_re, der=3, ext=0) \
         + 1.0j*interpolate.splev([0.0], spl_im, der=3, ext=0)
    
    #Lu = d2udr2 - l*(l+1)/r**2 * u
    if l == 0:
        lap = 0.0j*r
        lap[r>0] = d2udr2[r>0] / r[r>0]
        # r=0 limit
        lap[r==0] = d3udr3_0
    else:
        lap = 0.0j*r
        lap[r>0] = ( d2udr2[r>0] - l*(l+1)/r[r>0]**2 * u[r>0] ) / r[r>0]

    #### DEBUG
    import matplotlib.pyplot as plt
    plt.cla()
    plt.clf()
    
    d3udr3 =  interpolate.splev(r, spl_re, der=3, ext=0) \
              + 1.0j*interpolate.splev(r, spl_im, der=3, ext=0)

    
    plt.plot(r, lap.real, label="l=%d" % l)
    plt.plot(r, lap.imag, label="l=%d" % l)
    plt.plot(r, d2udr2.real, label="d2udr2")
    plt.plot(r, d3udr3.real, label="d3udr3")
    plt.plot([0.0], [d3udr3_0.real], "o")
    plt.legend()
    plt.show()
    ####
        
    return lap

def analyze_cusp(r, u, Lu):
    """
    The cusp of the wavefunction at the nucleus can provide information
    about the nuclear charge and the local energy of the wavefunction.
    At the cusp the Schroedinger equation

             __2       Z
        -1/2 \/ phi = --- phi + E phi        (1)
                       r

    has to be fulfilled. If we know the Laplacian and the wavefunction
    around r=0, we can find Z and E.

    Parameters
    ----------
    r           : 1d numpy array with grid points in decreasing order,
                  i.e. r[0] < r[1] < ... r[-1]
    u           : 1d numpy array with function values of radial wavefunction
                  u(r)=r*phi(r) at the grid points
    Lu          : radial part of Laplacian
                       __2
                       \/  phi = Lu/r

    Returns
    -------
    Z           : float, nuclear charge deduced from cusp of wavefunction
    E           : float, local energy deduced from cusp
    """
    #                      __2
    # Using phi = r*u  and \/ phi = L.u/r  equation (1) turns into
    #              Z
    #   L.u = -2 (--- + E) u         (2)
    #              r

    # Radii r < rs are considered small enough so that eqn. (1) is
    # valid.
    rs = 0.05
    # fit 3rd order polynomial to u(r) for small r
    #                  i
    #  u(r) sum    C  r              (no constant term)
    #          i=1  i
    # `polyfit` returns the coefficient vector in the wrong order.
    C  = np.polyfit(r[r < rs],   u[r < rs], 3)[::-1].real
    print "C = %s" % C
    # In the Taylor expansion of u(r) there must be no constant
    # term, for if there were a constant term, i.e.
    #   u(r) = u0 + C1 r + ...
    # phi(r) = u/r would not be finite at the origin.
#    assert abs( C[0] ) < 1.0e-5
    # fit 3rd order polynomial to L.u
    #                  i
    #  u(r) sum    D  r
    #          i=0  i
    D  = np.polyfit(r[r < rs],  Lu[r < rs], 3)[::-1].real
    print "D = %s" % D
    # Substituting these expansions for u and L.u into eqn. (2)
    # leads to the
    #
    #    Z = -1/2 D0/C1
    #    E =  1/(2 C1) ( C2/C1 D0 - D1 )
    #

    # nuclear charge
    Z = -0.5 * D[0]/C[1]
    # local energy
    E = 1.0/(2*C[1]) * ( C[2]/C[1] * D[0] - D[1] )

    print "nuclear charge  Z= %e" % Z
    print "local energy    E= %e" % E
    
    return Z, E


def multicenter_residual(f, potential_ee, e, 
                         atomic_coordinates, atomic_numbers,
                         lebedev_order=23, radial_grid_factor=1):
    """
    compute the action of the residual operator R=(H-e) on a function f
                               __2
      (H-E) f(x,y,z) =   -1/2  \/  f(x,y,z)  +  (V(x,y,z) - e) f(x,y,z)

    numerically on a multicenter spherical grid. The potential V is assumed to
    have singularities at the atomic nuclei which are the centers of the grid:

                       Z(I)
         V = sum   - --------  +   V      +  V
                I    |r-R(I)|       coul      xc

    Only V_coul+V_xc need to be provided, V_nuc is implied in the definition of the grid.
    The action of the Laplacian and the nuclear potential energy are computed 
    together, so as to avoid numerical errors at the nuclear positions where
    two large numbers of opposite signs (T phi ---> +oo, Vnuc phi ---> -oo) 
    would have to be added. 
    The algorithm is very similar to but simpler than `multicenter_poisson()`: 
    The function f is decomposed using the weight functions wI(r)

         f(r) = sum_I fI(r)             with  fI(r) = wI(r) f(r)

    Each part fI(r) is expanded into spherical harmonics, fI = sum_lm fI_lm(r) Y_lm(th,ph),
    and the kinetic operator and nuclear attraction at center I are applied in spherical harmonics,

             __2   Z(I)     1         d^2          L
       -1/2  \/  - ---- =  --- ( -1/2 --- r + 1/2 --- - Z ) .
                     r      r         dr^2         r

    At the end the parts are summed and the rest of the effective potential, Vcoul+Vxc, and the
    nuclear attraction to the other centers J!=I is added:
              __2                             __2   Z(I)                 Z(J)
       (-1/2  \/  + V_nuc) f = sum_I ( { -1/2 \/  - ---- } fI - sum      ---- fI )
                                                     rI            J!=I   rJ

                                + (Vcoul + Vxc) f
    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the function at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    potential_ee       : callable potential_ee(x,y,z) evaluates the effective potential 
                         WITHOUT the nuclear potential, only Vcoul + Vxc. The nuclear potential
                         is constructed from the atomic coordinates and numbers.
    e                  : float, energy
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    resf               : callable, resf(x,y,z) evaluates the residual (H-e) f
    """    
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
        f = mu
        for ik in range(0, k):
            f = 1.5 * f -0.5 * f**3
        return 0.5*(1-f)
    
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
    radial_hamiltonians = []
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
        # from variable transformation
        g = 2 * rm**3 * np.sqrt(((1+xr)/(1-xr)**3)**3)
        radial_weights *= g
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # cartesian coordinates of grid
        x = (outerN(r, sc) + atomic_coordinates[0,I])
        y = (outerN(r, ss) + atomic_coordinates[1,I])
        z = (outerN(r, c ) + atomic_coordinates[2,I])
        weights = outerN(radial_weights, 4.0*np.pi * angular_weights)
        #
        Npts = Nr*Nang
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # evaluate function on the grid
        fI = wr * f(x,y,z)

        # compute local atomic Hamiltonian
        #  (I)          __2   Z(I)    (I)
        # H    = { -1/2 \/  - ---- } f    
        #                      r
        #
        
        # expand function f_I(r) into spherical harmonics
        #
        #  f_I(x,y,z) = sum_l sum_{m=-l}^{l}  f_lm(r) Y_{l,m}(th,ph)
        #
        # and apply atomic Hamiltonian (kinetic + nuclear attraction)
        # to each (l,m) component
        radial_functions.append( {} )
        radial_hamiltonians.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            fI_lm = 4.0*np.pi * np.sum(fI*wYlm, axis=-1)

            spline_lm_real = interpolate.splrep(zr, fI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, fI_lm.imag, s=0)
            radial_functions[-1][(l,m)] = spline_lm_real, spline_lm_imag
            
            omx = 1-xr
            opx = 1+xr
            # coefficients for Laplacian operator after coordinate transformation
            c0 = -l*(l+1)/rm**2 * (omx/opx)**2
            c1 = 1.0/(4.0*np.pi*rm**2) * omx**(2.5) * (1.0+opx) / opx**(1.5)
            c2 = 1.0/(4.0*np.pi**2*rm**2) * omx**3 / opx

            #  
            # f_lm(r) = 1/r u_lm
            #
            uI_lm = r*fI_lm

            # apply Laplacian operator
            #            d^2     l(l+1)
            # L u_lm = ( ----  - ------ ) u_lm
            #            dr^2     r^2
            LuI_lm = radial_difop(c0,c1,c2, uI_lm)
            
            # __2  (I)              1    d^2     l(l+1)
            # \/  f    = sum_{l,m} --- ( ----  - ------ ) u_lm(r) Y_lm(th,ph)
            #                       r    dr^2     r^2
            #                         
            #          = sum_{l,m} 1/r LuI_lm(r) Y_lm(th,ph)
            #

            # locally the Hamiltonian is dominated by the kinetic and nuclear
            # energy of center I
            #      __2  (I)    Z(I)  (I)
            # -1/2 \/  f    -  ---- f    = sum_{l,m} 1/r  HuI_lm(r)  Y_lm(th,ph)
            #                   r             
            HuI_lm = -0.5 * LuI_lm - atomic_numbers[I] * fI_lm
            
            # z-array is sorted in assending order,
            #   0 < z[0] < z[1] < ... < z[-1] < 1
            # while associated r-array is sorted in descending order
            #   infinity > r[0] > r[1] > ... > r[-1] > 0

            # It seems better to spline HuI_lm instead of HuI_lm/r,
            # because the first function is smoother and has no
            # singularity at r=0.
            spline_lm_real = interpolate.splrep(zr, HuI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, HuI_lm.imag, s=0)
            radial_hamiltonians[-1][(l,m)] = spline_lm_real, spline_lm_imag
            
            if m == -(Lmax-1)/2:
                break

    def residual_func(x,y,z):
        """
        function for evaluating the residual (H-e) f
                 __2
           -1/2  \/  f  + (Vnuc + Vcoul + Vxc - e) f
        """
        res = 0j*x
        # build the total residual from parts
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:

                #              __2   Z(I)
                # sum    -1/2 (\/  - ---- ) fI(x,y,z)
                #    I                rI
                spline_lm_real, spline_lm_imag = radial_hamiltonians[I][(l,m)]
                HuI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                res += HuI_lm/rI * Ylm

                #                   Z(J)
                # sum   sum       - ----  fI(x,y,z)
                #    I     J!=I      rJ
                #
                spline_lm_real, spline_lm_imag = radial_functions[I][(l,m)]
                # interpolate
                fI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                # nuclear potential without atom I
                Vnuc = 0*x
                for J in range(0, Nat):
                    if I == J:
                        continue
                    xJ = x - atomic_coordinates[0,J]
                    yJ = y - atomic_coordinates[1,J]
                    zJ = z - atomic_coordinates[2,J]
                    rJ = np.sqrt(xJ*xJ+yJ*yJ+zJ*zJ)
                    
                    Vj = -atomic_numbers[J]/rJ
                    Vnuc += Vj

                res += Vnuc * fI_lm * Ylm
                    
                if m == -(Lmax-1)/2:
                    break


        #  (Vcoul + Vxc - e) f(x,y,z)
        res += (potential_ee(x,y,z) - e) * f(x,y,z)
        
        return res.real
    
    return residual_func


def multicenter_inhomogeneous_schroedinger(potential, source, energy,
                                           atomic_coordinates, atomic_numbers, lebedev_order=23, radial_grid_factor=1):
    """
    solve the inhomogeneous Schroedinger equation
           __2
     -1/2  \/  psi(r) + (V - E) psi(r) = source(r)

    numerically on a multicenter spherical grid (see Ref. [3])

    Parameters
    ----------
    potential          : callable, potential(x,y,z) should evaluate the potential V at the
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    source             : callable, source(x,y,z) should evaluate the right hand side of
                         Schroedinger's equation at the grid points
    energy             : float, energy E
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    psi                : callable, psi(x,y,z) evaluates the wavefunction which is the solution
                         to the inhomogeneous Schroedinger equation
    """    
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
        f = mu
        for ik in range(0, k):
            f = 1.5 * f -0.5 * f**3
        return 0.5*(1-f)

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

    radial_wavefunctions = []
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points if requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
        # from variable transformation
        g = 2 * rm**3 * np.sqrt(((1+xr)/(1-xr)**3)**3)
        radial_weights *= g
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # cartesian coordinates of grid
        x = (outerN(r, sc) + atomic_coordinates[0,I])
        y = (outerN(r, ss) + atomic_coordinates[1,I])
        z = (outerN(r, c ) + atomic_coordinates[2,I])
        weights = outerN(radial_weights, 4.0*np.pi * angular_weights)
        #
        Npts = Nr*Nang
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # evaluate source term on the grid
        sI = wr * source(x,y,z)

        # solve the inhomogeneous Schroedinger equation around each center I
        #
        #      __2    (I)     (sph)            (I)         (I)
        # -1/2 \/  psi    + (V     (r) - E) psi    = source 
        #
        
        # In order to decouple different (l,m) channels the potential
        # is spherically averaged around atom I
        #
        #    (sph)              /                             (I)
        #   V     (r) = 1/(4pi) | wI V dOmega = 1/sqrt(4 pi) V   (r)
        #                       /                             l=0,m=0
        # differential for solid angle 
        dOmega = 4*np.pi * outerN(np.ones(Nr), angular_weights)
        # 
        vI = wr * potential(x,y,z)
        # After spherical averaging only the l=0,m=0 component survives
        v_sph = 1.0/(4*np.pi) * np.sum(vI*dOmega, axis=-1)
        
        # expand source term s_I(r) into spherical harmonics
        #
        #  s_I(x,y,z) = sum_l sum_{m=-l}^{l}  s_lm(r) Y_{l,m}(th,ph)
        #
        # and solve Schroedinger equation for each (l,m) component
        radial_wavefunctions.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            sI_lm = 4.0*np.pi * np.sum(sI*wYlm, axis=-1)

            # boundary conditions for u(z)
            #  z=0 <-> x=1  <-> r=+infinity
            #  z=1 <-> x=-1 <-> r=0
            #      (I)
            # 1/r u    (r) -------> 0        exponentially fast
            #      l,m      r->oo
            u0 = 0.0
            u1 = 0.0
            # variable transformation
            omx = 1-xr
            opx = 1+xr
            #     l*(l+1)        (sph)
            # - ( ------- + 2 [ V     (r) - energy ] ) u = c0 u
            #      r^2
            c0 = -l*(l+1.0)/rm**2 * (omx/opx)**2 - 2*(v_sph - energy)
            # The transformation of partial derivatives from
            # r-coordinates to z-coordinates is accomplished as
            #  d^2 u      d u      d^2 u
            #  ----- = c1 --- + c2 -----
            #  d r^2      d z      d z^2
            c1 = 1.0/(4.0*np.pi*rm**2) * omx**(2.5) * (1+opx) / opx**(1.5)
            c2 = 1.0/(4.0*np.pi**2*rm**2) * omx**3 / opx
            # source term on right hand side
            #   -2 r sI_lm  = -2 rm (1+x)/(1-x) sI_lm
            rhs = -2.0 * rm * opx / omx * sI_lm
            
            #               d^2 u     l*(l+1)        (sph)
            # Now solve     -----  - (------- + 2 [ V     (r) - energy ] ) u = -2 r s
            #               d r^2       r^2
            #
            # which in z-coordinates turns into
            #               d u      d^2 u
            #     c0 u + c1 --- + c2 ----- = rhs
            #               d z      d z^2 
            u_lm = solve_radial_dgl(c0,c1,c2,rhs, u0, u1)
            # r*psi(r) = u(r)
            psiI_lm = u_lm/r
            
            # z-array is sorted in assending order,
            #   0 < z[0] < z[1] < ... < z[-1] < 1
            # while associated r-array is sorted in descending order
            #   infinity > r[0] > r[1] > ... > r[-1] > 0

            spline_lm_real = interpolate.splrep(zr, psiI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, psiI_lm.imag, s=0)
            radial_wavefunctions[-1][(l,m)] = spline_lm_real, spline_lm_imag

            if m == -(Lmax-1)/2:
                break

    def wavefunction(x,y,z):
        """
        function for evaluating the solution psi(x,y,z)
        """
        psi = 0j*x
        # build the total wavefunction from the solutions to the subproblems
        #  psi = sum_n  psi^(n)
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:
                
                spline_lm_real, spline_lm_imag = radial_wavefunctions[I][(l,m)]
                # interpolate
                psiI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)
                psi += psiI_lm*Ylm

                if m == -(Lmax-1)/2:
                    break

        return psi.real
    
    return wavefunction

def multicenter_continuum_schroedinger(potential, energy, charge, 
                                       atomic_coordinates, atomic_numbers,
                                       lebedev_order=23, radial_grid_factor=1):
    """
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CONCEPTUAL MISTAKE:
                              (I)
         If   (T + V  - E) psi    = 0       for all cells I
                    I
    
         adding the equations for all I leads to
                            (I)             (I)
            (T - E) sum  psi    + sum  V  psi   = 0
                       I             I  I

         which is not equivalent to
                                      (J)
            (T + sum  V  - E) sum  psi    = 0         
                    I  I         J  

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    solve the Schroedinger equation
           __2
     -1/2  \/  psi(r) + (V - E) psi(r) = 0                for E > 0

    numerically on a multicenter spherical grid

    Parameters
    ----------
    potential          : callable, potential(x,y,z) should evaluate the potential V at the
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    energy             : float, energy E
    charge             : asymptotically the potential takes the form 
                                    r-> oo
                                V --------->  -charge/r

    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    psi                : callable, psi(x,y,z) evaluates the wavefunction 
    """
    # energy and wavevector
    E = energy
    k = np.sqrt(2*E)
    
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, kk=3):
        f = mu
        for ik in range(0, kk):
            f = 1.5 * f -0.5 * f**3
        return 0.5*(1-f)

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
            
    radial_wavefunctions = []
    phase_shifts = []
    radial_grid = "linear" #"exponential" #"linear"
    for I in  range(0, Nat):
        # radial grid
        Nr = 500001
        rmin = 0.0
        rmax = 20.0
        if radial_grid == "linear":
            # The first point, where r=0, is excluded to avoid
            # the Coulomb singularity
            r = np.linspace(rmin, rmax, Nr)
        elif radial_grid == "exponential":
            #
            # Variable transformation
            #    r = exp(rho)     <=>   rho = log(r)
            #
            rmin = 1.0e-4
            rho = np.linspace(np.log(rmin), np.log(rmax), Nr)
            r = np.exp(rho)
        else:
            raise ValueError("`radial_grid` should be 'linear' or 'exponential'")
            
        # cartesian coordinates of grid
        x = (outerN(r, sc) + atomic_coordinates[0,I])
        y = (outerN(r, ss) + atomic_coordinates[1,I])
        z = (outerN(r, c ) + atomic_coordinates[2,I])
        #
        Npts = Nr*Nang
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # solve the continuum Schroedinger equation around each center I
        #
        #      __2    (I)     (sph)            (I)    
        # -1/2 \/  psi    + (V     (r) - E) psi    = 0
        #
        
        # In order to decouple different (l,m) channels the potential
        # is spherically averaged around atom I
        #
        #    (sph)              /                             (I)
        #   V     (r) = 1/(4pi) | wI V dOmega = 1/sqrt(4 pi) V   (r)
        #                       /                             l=0,m=0
        # differential for solid angle 
        dOmega = 4*np.pi * outerN(np.ones(Nr), angular_weights)
        # 
        vI = wr * potential(x,y,z)
        # After spherical averaging only the l=0,m=0 component survives
        v_sph = 1.0/(4*np.pi) * np.sum(vI*dOmega, axis=-1)

        """
        ### DEBUG
        print "atom I=%d" % I
        import matplotlib.pyplot as plt
        # 
        plt.plot(r, v_sph, ls="-", label=r"$V^{(sph)}_{%d}$" % I)
        plt.plot(r, -charge/r, ls="-.", label=r"$-Q/r$")
        plt.plot(r, -1/r, ls="--", label=r"$-1/r$")
        # show positions of grid points
        plt.plot(r, 0.0*r+1, "x", label="grid points")
        plt.legend()
        plt.show()
        ###
        """
        
        # Now we solve the Schroedinger equation for each l component
        radial_wavefunctions.append( {} )
        phase_shifts.append( {} )
        for l in range(0, Lmax/2+1):

            if radial_grid == "linear":
                #
                #              d^2 u     l*(l+1)        (sph)
                # We solve     -----  = (------- + 2 [ V     (r) - energy ] ) u 
                #              d r^2       r^2
                #
                #                     = F(r) u(r)
                #
                # using Numerov's method
                #
                #
                # The radial wavefunction is propagated outward starting
                # from r=0. We need the values of the wavefunction at
                # the first two grid points in order to initiate the propagation.
                # Around r=0 the differential equation for u simplifies to
                #
                #  d^2 u   l*(l+1)                                       l+1
                #  ----- = ------ u      with the solution   u(r) = C * r
                #  d r^2     r^2
                #
                # The unknown constant C is determined from matching the wavefunction
                # to a Coulomb wave at large r.

                if l == 0:
                    F = 2*(v_sph - energy)

                    # initial values at first two points
                    u0 = r[0]**(l+1)
                    u1 = r[1]**(l+1)

                    # propagate radial wavefunction
                    u_l = integrate_numerov(r, F, u0, u1)
                else:
                    # l > 0
                    # We have to skip the first point, where r=0
                    # so as not to divide by zero.
                    F = l*(l+1)/r[1:]**2 + 2*(v_sph[1:] - energy)

                    # initial values
                    u0 = r[0]**(l+1)  # == 0.0
                    u1 = r[1]**(l+1)
                    u2 = r[2]**(l+1)

                    u_l = 0*r

                    u_l[0]  = u0
                    # start the propagation at the second point
                    u_l[1:] = integrate_numerov(r[1:], F, u1, u2)
                    # 
                    
            elif radial_grid == "exponential":
                #
                # Using rho=log(r) as radial coordinate and solving for
                #
                #   v(rho) = r^{-1/2} u(r)
                #
                # leads to the transformed radial equation:
                #
                #   d^2 v
                #   ------  = [ 1/4 + l*(l+1) - 2 exp(2 rho) (E - V) ] v
                #   drho^2
                #
                #           = F(rho) v

                # using exp(2 rho) = r**2
                F = 0.25 + l*(l+1) - 2.0 * r**2 * (energy - v_sph)
                
                # initial values at first points
                #   v(rho) = exp(-rho/2) u(r(rho)) = r^{-1/2} u(r)
                #
                v0 = r[0]**(l+0.5)
                v1 = r[1]**(l+0.5)

                # propagate radial wavefunction
                v_l = integrate_numerov(rho, F, v0, v1)

                # transform back to r
                u_l = np.sqrt(r) * v_l
                
            # To determine the phase shift and normalization factor the
            # radial wavefunction has to be matched to the asymptotic solution,
            # which is a regular of irregular Coulomb function

            # The matching is done at the last two grid points
            r1 = r[-2]
            r2 = r[-1]
            # where the radial wavefunction has the following values
            ul1 = u_l[-2]
            ul2 = u_l[-1]
            # regular and irregular Coulomb functions at r1 
            solF1 = float(mpmath.coulombf(l, -charge/k, r1*k))/r1
            solG1 = float(mpmath.coulombg(l, -charge/k, r1*k))/r1
            # and r2
            solF2 = float(mpmath.coulombf(l, -charge/k, r2*k))/r2
            solG2 = float(mpmath.coulombg(l, -charge/k, r2*k))/r2
            
            Det = solF1*solG2 - solF2*solG1
            aa = (solG2 * ul1/r1 - solG1 * ul2/r2)/Det
            bb = (-solF2* ul1/r1 + solF1 * ul2/r2)/Det
            # phase shift
            delta = np.arctan2(-bb,aa)
            # scale factor
            scale_fac = 1.0/np.sqrt(aa*aa + bb*bb)

            # normalize u
            u_l *= scale_fac

            # spline the radial wavefunction u(r) = r*psi(r)
            #
            # There are no values for the wavefunction beyond r=rmax, however, at large
            # distances the wavefunction should agree with a Coulomb wave. Outside
            # the interpolation range u_l is substituted by ug_asymptotic, which is a shifted
            # Coulomb wave matched to u_l at r=rmax.

            # For r < rmin, the solution should be u(r) = scale_fac * r^(l+1)
            #
            
            Rl_splined = spline_radial_wavefunction(r, u_l,
                                                    ug_asymptotic=coulombf_ufunc(E, charge, l, delta))

            ### DEBUG
            if l == 0:
                print "scale_fac= %s" % scale_fac
                print u_l[1]/r[1]
                print Rl_splined(np.array([r[0], r[1], r[2]]))
                import matplotlib.pyplot as plt
                plt.cla()
                plt.clf()
                plt.plot(r[:20], Rl_splined(r[:20]))
                rnew = np.linspace(r[0], r[20], 1000)
                plt.plot(rnew, Rl_splined(rnew), ls="-.")
                plt.show()
            ###
            
            radial_wavefunctions[-1][l] = Rl_splined
            
            # save phase shifts for partial wave l of atom I
            phase_shifts[-1][l] = delta
            
    def wavefunction(x,y,z, l=0, m=0):
        """
        function for evaluating the continuum solution psi(x,y,z) which 
        has asymptotic angular momentum (l,m)
        """
        psi = 0j*x
        # build the total wavefunction from the solutions to the subproblems
        #  psi = sum_n  psi^(n)
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            for Ylm,ll,mm in sph_it:

                if ll == l and mm == m:
                    # radial wavefunction with angular momentum l
                    # for potential V_I = w_I V
                    Rl_splined = radial_wavefunctions[I][l]
                    # interpolate
                    psiI_l = Rl_splined(rI)
                    psi += psiI_l*Ylm

                    break
                
                if m == -(Lmax-1)/2:
                    break

        return psi.real
    
    return wavefunction

def integrate_numerov(r, f, u0, u1):
    """
    Integrate the second order linear differential equation

       d^2 u
       ----- = f(r) u
       dr^2 

    numerically on an equidistant grid r using Numerov's method.
    The propagation is initiated by providing the values of u
    at the first two grid points.

    Parameters
    ----------
    r           :   1d numpy array, equidistant grid
    f           :   1d numpy array, f(r) evaluated at the grid points
    u0, u1      :   floats, initial values u0 = u(r[0]) and u1 = u(r[1])

    Returns
    -------
    u           :   1d numpy array, solution u(r) at the grid points

    """
    h = np.ediff1d(r, to_end=r[-1]-r[-2])
    w = np.zeros(r.shape)
    w[0] = (1.0 - 1.0/12.0*pow(h[0],2) * f[0])*u0
    w[1] = (1.0 - 1.0/12.0*pow(h[1],2) * f[1])*u1    
    # propagate solution w
    for i,ri in enumerate(r):
        if i < 2:
            # first two values given as initial conditions
            continue
        w[i] = 2.0 * (1.0 + 5.0/12.0*pow(h[i],2)*f[i-1])/(1.0 - 1.0/12.0*pow(h[i],2)*f[i-1]) * w[i-1] - w[i-2]
    # transform back u -> w
    u = 1.0/(1.0 - pow(h,2)/12.0 * f) * w
    # plug in initial values at first two points
    u[0] = u0
    u[1] = u1
    
    return u

def spline_radial_wavefunction(rg,ug,
                               ug_asymptotic=lambda r: 0*r):
    """
    Create a spline of the radial part R(r) from a 
    a discretized version u of u(r)=r*R(r) on the radial grid r.

    Parameters:
    -----------
    rg: numpy array with r-points in increasing order, 
        rg[0] < rg[1] < ... < rg[-1]
    ug: numpy array with values of u on r

    Optional:
    ---------
    ug_asymptotic :   function u(r) for large r > rg.max(), 
                      outside the interpolation range this asymptotic form of the wavefunction
                      is used, for bound wavefunctions the function should give 0, while 
                      for scattering wavefunctions the oscillating solution with the correct phase shift 
                      should be returned.

    Returns:
    --------
    a callable function R(r) = u(r)/r
    """
    rmin, rmax = min(rg[rg > 0.0]), max(rg)

    R = ug[rg > rmin]/rg[rg > rmin]
    # spline Rg
    tck = interpolate.splrep(rg[rg > rmin],R, s=0)

    # On the interval [0,rmin] R is approximated by a polynomial
    #   p(r) = a + 1/2 b r^2 + 1/24 c r^4
    # The coefficients are chosen such that
    #   R''(0) = 0
    # and that R is continuous and continuously differentiable
    # at r=rmin

    # evaluate R and its derivatives at r=rmin
    # R(rmin)
    R0 = interpolate.splev([rmin], tck, der=0)
    # R'(rmin)
    R1 = interpolate.splev([rmin], tck, der=1)
    # R''(rmin)
    R2 = interpolate.splev([rmin], tck, der=2)

    # From
    #   p(rmin)   = R(rmin)
    #   p'(rmin)  = R'(rmin)
    #   p''(rmin) = R''(rmin)
    # follows that
    
    c = 3.0/rmin**3 * (rmin*R2 - R1)
    b = R2 - 0.5*c*rmin**2
    a = R0 - 0.5*b*rmin**2 - 1.0/24.0 * c * rmin**4
    
    def R_poly(r):
        return a + 1.0/2.0 * b * r**2 + 1.0/24.0 * c * r**4
    
    def R_func(r):
        # The input array `r` can have any shape.
        #
        # keep track of original shape
        shape = r.shape
        # flatten arrays 
        r = r.flatten()
        Rg = np.zeros(r.shape)

        # below interpolation range
        Rg[r < rmin] = R_poly(r[r < rmin])
        # inside interpolation range
        Rg[(rmin <= r) & (r < rmax)] = interpolate.splev(r[(rmin <= r) & (r < rmax)], tck, der=0)
        # above interpolation range
        Rg[rmax <= r] = ug_asymptotic(r[rmax <= r]) / r[rmax <= r]

        # At the end the array is cast to the original shape again
        Rg = np.reshape(Rg, shape)
        
        return Rg

    return R_func

def monomial_ufunc(scale_fac, l):
    """
    create function     
          scale_fac * r**l
    """
    def u(r):
        return scale_fac * r**l

    return u
    
def coulombf_ufunc(E, charge, l, delta):
    """
    create regular Coulomb function that can operate on arrays
    """
    k = np.sqrt(2*E)
    
    def coulombf(r):
        kr = k*r
        # radial part R_{E,l}(r) of wavefunction
        rhos = kr.flatten()
        fs = []
        for rho in rhos:
            fi = mpmath.coulombf(l, -charge/k, rho+delta)
            fs.append(complex(fi).real)
        cf = np.reshape(np.array(fs), kr.shape)
        return cf

    return coulombf



def multicenter_spherical_remainder(potential, 
                                       atomic_coordinates, atomic_numbers,
                                       lebedev_order=23, radial_grid_factor=1):
    """
    When solving the inhomogeneous Schroedinger equation on a multicenter
    grid, the potential is replaced by a spherical average around each
    atom in order to decouple different angular momenta.

    The wavefunction is thus a not a solution for the full potential

      V  =  sum  V                  V (x,y,z) = w (x,y,z) V(x,y,z)
               I  I                  I           I
    
    but for the simpler potential

       (sph)          (sph)                /
      V      =  sum  V      = sum  1/(4pi) | V (r-R ) dOmega
                   I  I          I         /  I    I

    To assess the quality of this approximation, we need to know the remaining
    non-spherical potential

       (rem)        (sph)
      V      = V - V

    This function splits the potential into the two parts, V^(sph) and V^(rem).

    Parameters
    ----------
    potential          : callable, potential(x,y,z) should evaluate the potential V at the
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    spherical          : callable, spherical(x,y,z) evaluates V^(sph)
    remainder          : callable, remainder(x,y,z) evaluates the difference between the
                         full potential V and the simple potential V^(sph)
    """    
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
        f = mu
        for ik in range(0, k):
            f = 1.5 * f -0.5 * f**3
        return 0.5*(1-f)

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

    radial_spherical = []
    radial_remainder = []
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
        # from variable transformation
        g = 2 * rm**3 * np.sqrt(((1+xr)/(1-xr)**3)**3)
        radial_weights *= g
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # cartesian coordinates of grid
        x = (outerN(r, sc) + atomic_coordinates[0,I])
        y = (outerN(r, ss) + atomic_coordinates[1,I])
        z = (outerN(r, c ) + atomic_coordinates[2,I])
        weights = outerN(radial_weights, 4.0*np.pi * angular_weights)
        #
        Npts = Nr*Nang
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # evaluate full potential on the grid
        vI = wr * potential(x,y,z)
        
        # expand potential v_I(r) into spherical harmonics
        #
        #  v_I(x,y,z) = sum_l sum_{m=-l}^{l}  v^(I)_lm(r) Y_{l,m}(th,ph)
        #
        # and separate the (l=0,m=0) component
        radial_spherical.append( {} )
        radial_remainder.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            vI_lm = 4.0*np.pi * np.sum(vI*wYlm, axis=-1)
            
            spline_lm_real = interpolate.splrep(zr, vI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, vI_lm.imag, s=0)

            if l == 0 and m == 0:
                radial_spherical[-1][(l,m)] = spline_lm_real, spline_lm_imag
            else:
                radial_remainder[-1][(l,m)] = spline_lm_real, spline_lm_imag

            if m == -(Lmax-1)/2:
                break

    def spherical(x,y,z):
        """
        function for evaluating the V^(sph)
        """
        vsph = 0j*x
        # build the total potential from the individual parts
        #  V^(sph) = sum_I  V^(sph)_I
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:
                if l == 0 and m == 0:
                    # only (0,0) component in V^(sph)
                    spline_lm_real, spline_lm_imag = radial_spherical[I][(l,m)]
                    # interpolate
                    vsphI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                               + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)
                    vsph += vsphI_lm*Ylm

                    break

        return vsph.real
            
    def remainder(x,y,z):
        """
        function for evaluating the remainder V - V^(sph)
        """
        vrem = 0j*x
        # build the total remainder from the individual parts
        #  V^(rem) = sum_I  V^(rem)_I
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:
                if l == 0 and m == 0:
                    # no (0,0) component in remainder
                    continue
                
                spline_lm_real, spline_lm_imag = radial_remainder[I][(l,m)]
                # interpolate
                vremI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)
                vrem += vremI_lm*Ylm

                if m == -(Lmax-1)/2:
                    break

        return vrem.real
    
    return spherical, remainder


def multicenter_operation(fs, op, 
                          atomic_coordinates, atomic_numbers,
                          lebedev_order=23, radial_grid_factor=1):
    """
    given a list of functions fs = [f1,f2,...,fn] perform the n-ary operation

         h = op(f1,f2,...,fn)

    to create a new function h. Examples for binary operations are
    addition

        op = lambda fs: a*fs[0]+b*fs[1]

    or multiplication

        op = lambda fs: fs[0]*fs[1]

    The new function h is defined by radial splines on the multicenter grid.
    
    At first it appears easier to simplify define a new function h as

        def h(x,y,z):
            return op(f1(x,y,z),f2(x,y,z),...)

    However, every evaluation of h requires the evaluation of all f1,f2,...,fn. 
    If we would use this approach to build up more complicated functions, the effort
    for evaluating the combined functions would grow with the number of elementary
    functions used to define them. By interpolating h over the grid, the effort
    for evaluating h remains the same as that for each fi. 


    Parameters
    ----------
    fs                :  list of callables, fs[i](x,y,z) should evaluate the i-th function
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    op                 : operator, callable or lambda function, takes a list of numpy grids
                         as a single argument, e.g. op([f1(x,y,z),f2(x,y,z)]) for a binary operator.
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    h                  : callable, h(x,y,z) evaluates the function h(x,y,z) = op([f1(x,y,z),f2(x,y,z),...])
    """    
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
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

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
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
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # evaluate h = op(f1,f2,...) on the grid
        
        # First we need to evaluate the functions fi on the grid
        fs_xyz = []
        for f in fs:
            fs_xyz.append( f(x,y,z) )
        # Then we pass the values of the functions as arguments to the operator
        hI = wr * op(fs_xyz)

        radial_functions.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            hI_lm = 4.0*np.pi * np.sum(hI*wYlm, axis=-1)

            spline_lm_real = interpolate.splrep(zr, hI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, hI_lm.imag, s=0)
            radial_functions[-1][(l,m)] = spline_lm_real, spline_lm_imag

            if m == -(Lmax-1)/2:
                break

    def h_func(x,y,z):
        """
        function for evaluating h = op(f,g)
        """
        h = 0j*x
        # sum over centers
        #  h = sum_I  h^(I)
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            sph_it = spherical_harmonics_it(thI,phI)

            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:
                
                spline_lm_real, spline_lm_imag = radial_functions[I][(l,m)]
                # interpolate
                hI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)
                h += hI_lm*Ylm

                if m == -(Lmax-1)/2:
                    break

        return h.real
    
    return h_func

def multicenter_gradient2(f,
                          atomic_coordinates, atomic_numbers,
                          cusps_separate=True,
                          lebedev_order=23, radial_grid_factor=1):
    """
    given a scalar function f(x,y,z) create a function for evaluating
    the dot product of the gradient of f with itself:
                        __     __
        sigma(x,y,z) = (\/ f).(\/ f)

    In the context of DFT, where f is the electron density, this quantity
    is needed for GGA functionals.

    Parameters
    ----------
    f                  : callable, f(x,y,z) should evaluate the function at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    cusps_separate     : if True, gradients of `f` around the atoms are calculated
                         separately to avoid numerical artifacts that arise from the cusps
                         if `f` represents the electron density.
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    sigma              : callable, sigma(x,y,z) evaluates the square of the gradient of f,
                         i.e. (grad f)^2
    """
    Nat = atomic_coordinates.shape[1]
    atomic_names = [atom_names[Z-1] for Z in atomic_numbers]
    ####
    # The computation of sigma = (grad rho)^2 is difficult
    # probably because at the nuclei the density has cusps, which
    # makes numerical differentiation very difficult. Around an
    # atom with nuclear charge Z_i and position R_i the electron density
    # behaves as
    #
    #   rho_i(r) = N exp(-2*Z_i*r)
    #
    # At the atom positions the density has a cusp and is approximately
    # spherically symmetric. We can split the total density into a sum
    # of spherically symmetric densities with cusps and a relatively
    # smooth remainder drho(r):
    #
    #   rho(r) = sum  rho (r) + drho(r)
    #               i    i
    #
    # The gradient (grad rho) is then calculated separately for drho
    # on a multicenter grid and for the atomic densities exploiting the
    # spherical symmetry to avoid the numerical problems at r=Ri.
    # We do not know the atomic density rho_i, therefore we replace this
    # term by a spherical average of the total density around atom i:
    #                                /
    #   rho  --->    g (r) = 1/(4pi) | rho(r+Ri) dOmega
    #      i          i              /
    # So that the total density becomes
    #
    #   rho(r) = sum   g (|r-Ri|)  +  [ rho(r) - sum  g (|r-Rj|) ]
    #               i   i                           j  j
    # To compute the gradient we first subtract the spherical average
    # of rho around each atom and add the contribution to the gradient
    # from the subtracted part at the end.
    #   __          __                                        __
    #   \/ rho(r) = \/ [ rho(r) - sum  g (|r-Rj|) ]   +  sum  \/ g (|r-Ri|)
    #                                j  j                   i     i

    # 1) First we calculate the spherical averages of f around each atom
    #    and spline the resulting functions g_i(z) and their derivatives
    #    g_i'(z) = d/dz g_i(z) using z-coordinates.
    # list of splines g_i
    radial_avg = []
    # list of splines of the derivative g_i' of the spherical average
    radial_deriv_avg = []
    for I in  range(0, Nat):
        # center around which the average is performed
        Zat = atomic_numbers[I]
        pos = atomic_coordinates[:,I]
        atom = (Zat, pos)
        #
        gI_func = spherical_average_func(atom, f,
                                         lebedev_order=lebedev_order,
                                         radial_grid_factor=radial_grid_factor)

        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # variable transformation
        omx = 1-xr
        opx = 1+xr

        # d f   d f
        # --- = --- (-1/(2*pi*rm) (1-x)^(3/2) (1+x)(-1/2))
        # d r   d z
        #          d f
        #     = c1 ---
        #          d z

        c0 = 0.0*omx
        c1 = -1.0/(2.0*np.pi*rm) * omx**(1.5) * opx**(-0.5)
        c2 = 0.0*omx

        gI = gI_func(r)
        gI_der = radial_difop(c0, c1, c2, gI)
        
        # interpolate gI(z) and gI'(z) on a grid
        spline_avg = interpolate.splrep(zr, gI, s=0)
        radial_avg.append(spline_avg)

        spline_der_avg = interpolate.splrep(zr, gI_der, s=0)
        radial_deriv_avg.append( spline_der_avg )

    """
    ### DEBUG
    import matplotlib.pyplot as plt
    r = np.linspace(0.0, 10.0, 1000)
    for I in range(0, Nat):
        plt.plot(r, interpolate.splev(r, radial_avg[I], der=0, ext=0), label="around atom %d" % (I+1))
    plt.legend()
    plt.show()
    ###
    """
    ####

    # 2) Define weight functions for fuzzy Voronoi decomposition
    #
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
        ff = mu
        for ik in range(0, k):
            ff = 1.5 * ff -0.5 * ff**3
        return 0.5*(1-ff)
    
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

    # 3) Evaluate function f on the multicenter grid and perform a spherical
    #    wave decomposition 
    
    # These lists contain for each center I a dictionary with
    # splines of the angular components f^I_lm(r) and (df/dr)^I_lm.
    radial_functions_f = []
    radial_functions_dfdr = []
    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
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
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # evaluate f on the grid
        fI = wr * f(x,y,z)

        if cusps_separate:
            ####
            # subtract spherical averages around each atom
            #
            #  f_i   = w (r) [ f(r) - sum  g (|r-Rj|) ]
            #           i                j  j
            for j in range(0, Nat):
                # rJ = |r-Rj]
                rJ = dist[:,:,j]
                # transform to z-coordinates
                rmJ = 0.5*slater_radii[atomic_names[j]]
                xrJ = (rJ-rmJ)/(rJ+rmJ)
                zrJ = np.arccos(xrJ) / np.pi
                
                gJ_spline = radial_avg[j]
                fI -= wr * interpolate.splev(zrJ, gJ_spline, der=0, ext=0)
            ####

        # Perform spherical wave decomposition of fI
        #     (I)              m=+l    (I)
        #    f     = sum    sum       f  (r)  Y (th,ph)
        #               l=0    m=-l     l,m    l,m
        # and compute the derivatives
        #    d(fI_{l,m})/dr
        # The radial functions fI_{l,m} and dfI_{l,m}/dr are splined
        # in z-coordinates for later use.

        # List of spherical wave components f^(I)_{l,m},
        # the spline for f^(I)_{l,m} is stored in   radial_functions_f[I][(l,m)]
        radial_functions_f.append( {} )
        # and the spline for df/dr^(I)_{l,m} is stored in   radial_functions_dfdr[I][(l,m)]
        radial_functions_dfdr.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            # project fI onto spherical harmonic Y^*_{l,m}
            fI_lm = 4.0*np.pi * np.sum(fI*wYlm, axis=-1)

            # spline angular components of f
            spline_lm_real = interpolate.splrep(zr, fI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, fI_lm.imag, s=0)
            radial_functions_f[-1][(l,m)] = spline_lm_real, spline_lm_imag
            
            # variable transformation
            omx = 1-xr
            opx = 1+xr

            # d f   d f
            # --- = --- (-1/(2*pi*rm) (1-x)^(3/2) (1+x)(-1/2))
            # d r   d z
            #          d f
            #     = c1 ---
            #          d z

            c0 = 0.0*omx
            c1 = -1.0/(2.0*np.pi*rm) * omx**(1.5) * opx**(-0.5)
            c2 = 0.0*omx
            
            dfIdr_lm = radial_difop(c0, c1, c2, fI_lm)

            # spline angular components of radial gradient df/dr
            spline_lm_real = interpolate.splrep(zr, dfIdr_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, dfIdr_lm.imag, s=0)
            radial_functions_dfdr[-1][(l,m)] = spline_lm_real, spline_lm_imag
            
            if m == -(Lmax-1)/2:
                break

    # 4) Define function for evalulating (grad f)^2
    def sigma_func(x,y,z):
        """
        function for evaluating sigma = (grad f).(grad f)

        In spherical coordinates the gradient is

        grad f = e  g  + e   g   + e   g 
                  r  r    th  th    ph  ph

        where e_r, e_th and e_ph are unit vectors and 

            g  = df/dr   ;   g  = 1/r df/dth   ;   g  = 1/(r sin(th)) df/dph
             r                th                    ph

        For resummation over the centers we need to convert the gradients to cartesian
        coordinates so that all unit vectors point into the same direction.
        """
        grad_x = 0.0j*x
        grad_y = 0.0j*x
        grad_z = 0.0j*x
        # sum over centers
        #   grad_x = sum_I grad_x^(I)
        #   etc.
        # 
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            # unit vectors in cartesian cordinates
            stI = np.sin(thI)
            ctI = np.cos(thI)
            spI = np.sin(phI)
            cpI = np.cos(phI)
            #
            unit_r  = [stI*cpI, stI*spI, ctI]
            unit_th = [ctI*cpI, ctI*spI, -stI]
            unit_ph = [-spI   , cpI    , 0*stI]

            cosI = ctI
            sinI = stI
            s2 = sinI**2

            #
            sph_it = spherical_harmonics_it(thI,phI)
            # inverse transformation rI -> zr
            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:
                
                # interpolate fI_{l,m}(x,y,z)
                spline_lm_real, spline_lm_imag = radial_functions_f[I][(l,m)]
                fI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                # interpolate dfI_{l,m}/dr(x,y,z)
                spline_lm_real, spline_lm_imag = radial_functions_dfdr[I][(l,m)]
                dfIdr_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                # compute g_r = df/dr(x,y,z) = sum_{l,m} dfIdr_{l,m}(r) Y_{l,m}(th,ph)
                grad_r = dfIdr_lm * Ylm


                if cusps_separate:
                    #####
                    if l == 0 and m == 0:
                        # add back the radial gradient due to spherical
                        # averaged part g_i that has been subtracted previously:
                        #   __     __                            __
                        #   \/ f = \/ [ f - sum_i g_i ]  + sum_i \/ g_i
                        #   
                        grad_gI_r = interpolate.splev(zr, radial_deriv_avg[I], der=0, ext=0)
                        grad_r += grad_gI_r

                    #####

                # compute g_th = 1/r df/dth = sum_{l,m} 1/r f_{l,m}(r) d/dth Y_{l,m}(th,ph)
                #  = sum_{l,m} [ 1/r f_{l,m} m cot(th)
                #                + Heaviside(-l+1 <= m) 1/r f_{l,m-1} sqrt((l-m+1)(l+m)) exp(-i ph) ] Y_{l,m}
                #
                # The formula for d(Y_lm)/dth is taken from http://functions.wolfram.com/05.10.20.0001.01
                #
                # To avoid dividing by 0, we first multiply all terms by sin(th) and at the end
                # divide those elements of the array where sin(th) != 0 by sin(th) again.
                # I am not sure if removing the singularity at th=0,pi in this way is correct.
                #
                grad_th = 1.0/rI * fI_lm * m * cosI * Ylm
                if -l+1 <= m:
                    # interpolate fI_{l,m-1}
                    spline_lm_real, spline_lm_imag = radial_functions_f[I][(l,m-1)]
                    fI_lmm1 = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                            + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                    grad_th += 1.0/rI * fI_lmm1 * np.sqrt((l-m+1)*(l+m)) * np.exp(-1.0j*phI) * sinI * Ylm
                # Those elements where division by sin(th) would lead to NaN are set to zero. 
                grad_th[s2 > 0.0] /= sinI[s2 > 0.0]
                grad_th[s2 <= 0.0] = 0.0
                    
                # compute g_ph = 1/(r*sin(th)) df/dph = sum_{l,m} (i m) 1/(r*sin(th)) f_{l,m} Y_{l,m}
                grad_ph = 1.0/rI * fI_lm * (1.0j*m) * Ylm
                # Those elements where division by sin(th) would lead to NaN are set to zero. 
                grad_ph[s2 > 0.0] /= sinI[s2 > 0.0]
                grad_ph[s2 <= 0.0] = 0.0
                
                # We cannot add the gradients belonging to different spherical coordinate systems
                # since the unit vectors point in different directions. 

                grad_x += unit_r[0] * grad_r + unit_th[0] * grad_th + unit_ph[0] * grad_ph
                grad_y += unit_r[1] * grad_r + unit_th[1] * grad_th + unit_ph[1] * grad_ph
                grad_z += unit_r[2] * grad_r + unit_th[2] * grad_th + unit_ph[2] * grad_ph

                if m == -(Lmax-1)/2:
                    break
                
        # Now we compute the square of the gradient vector. 
        sigma = grad_x**2 + grad_y**2 + grad_z**2
        
        return sigma.real
    
    return sigma_func


def multicenter_gradient_product(f1, f2, 
                                 atomic_coordinates, atomic_numbers,
                                 cusps_separate=True,
                                 lebedev_order=23, radial_grid_factor=1):
    """
    given two scalar functions f1(x,y,z) and f2(x,y,z) create a function 
    for evaluating the dot product of the gradients:
                          __      __
        sigma12(x,y,z) = (\/ f1).(\/ f2)

    This function is needed in the finite-volume variational method (L. Rouzo and Raseev (1984)),
    where we have to integrate sigma12 over a finite volume.

    Parameters
    ----------
    f1,f2              : callable, f1(x,y,z) and f2(x,y,z) should evaluate the functions at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]
    atomic_coordinates : numpy array with shape (3,Nat), atomic_coordinates[:,i] is the 
                         cartesian position of atom i
    atomic_numbers     : numpy array with shape (Nat)
    
    Optional
    --------
    cusps_separate     : if True, gradients of `f1` and `f2` around the atoms are calculated
                         separately to avoid numerical artifacts that arise from the cusps
                         of wavefunctions at the nuclei.
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor

    Returns
    -------
    sigma12            : callable, sigma12(x,y,z) evaluates the product of the gradients,
                         i.e. (grad f1).(grad f2)
    """
    Nat = atomic_coordinates.shape[1]
    atomic_names = [atom_names[Z-1] for Z in atomic_numbers]
    ####
    # This function is very similar to `multicenter_gradient2`, but much more difficult
    # to understand because lots of code is duplicate, so see also the comments there.
    # 
    
    # 1) First we calculate the spherical averages of f1 and f2 around each atom
    #    and spline the resulting functions g_i(z) and their derivatives
    #    g_i'(z) = d/dz g_i(z) using z-coordinates.
    # list of splines g_i
    radial_avg_1 = []         #  for f1
    radial_avg_2 = []         #  for f2
    # list of splines of the derivative g_i' of the spherical average
    radial_deriv_avg_1 = []   #  for f1
    radial_deriv_avg_2 = []   #  for f2
    for I in  range(0, Nat):
        # center around which the average is performed
        Zat = atomic_numbers[I]
        pos = atomic_coordinates[:,I]
        atom = (Zat, pos)
        # for f1
        gI_func_1 = spherical_average_func(atom, f1,
                                           lebedev_order=lebedev_order,
                                           radial_grid_factor=radial_grid_factor)
        # same for f2
        gI_func_2 = spherical_average_func(atom, f2,
                                           lebedev_order=lebedev_order,
                                           radial_grid_factor=radial_grid_factor)
        
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # radial grid points on interval [0,infinity]
        r = rm * (1+xr)/(1-xr)

        # variable transformation
        omx = 1-xr
        opx = 1+xr

        # d f   d f
        # --- = --- (-1/(2*pi*rm) (1-x)^(3/2) (1+x)(-1/2))
        # d r   d z
        #          d f
        #     = c1 ---
        #          d z

        c0 = 0.0*omx
        c1 = -1.0/(2.0*np.pi*rm) * omx**(1.5) * opx**(-0.5)
        c2 = 0.0*omx

        # for f1
        gI_1 = gI_func_1(r)
        gI_der_1 = radial_difop(c0, c1, c2, gI_1)
        # ... and for f2
        gI_2 = gI_func_2(r)
        gI_der_2 = radial_difop(c0, c1, c2, gI_2)

        # interpolate gI(z) and gI'(z) on a grid
        # for f1
        spline_avg_1 = interpolate.splrep(zr, gI_1, s=0)
        radial_avg_1.append(spline_avg_1)

        spline_der_avg_1 = interpolate.splrep(zr, gI_der_1, s=0)
        radial_deriv_avg_1.append( spline_der_avg_1 )
        # ... and for f2
        spline_avg_2 = interpolate.splrep(zr, gI_2, s=0)
        radial_avg_2.append(spline_avg_2)

        spline_der_avg_2 = interpolate.splrep(zr, gI_der_2, s=0)
        radial_deriv_avg_2.append( spline_der_avg_2 )
        
    ####

    # 2) Define weight functions for fuzzy Voronoi decomposition
    #
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)
    # for nuclear weight functions
    def s(mu, k=3):
        ff = mu
        for ik in range(0, k):
            ff = 1.5 * ff -0.5 * ff**3
        return 0.5*(1-ff)
    
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

    # 3) Evaluate functions f1 and f2 on the multicenter grid and perform a spherical
    #    wave decomposition 
    
    # These lists contain for each center I a dictionary with
    # splines of the angular components f^I_lm(r) and (df/dr)^I_lm.
    #   for f1
    radial_functions_f1 = []
    radial_functions_df1dr = []
    # ... and for f2
    radial_functions_f2 = []
    radial_functions_df2dr = []

    for I in  range(0, Nat):
        # radial grid
        Nr = number_of_radial_points(atomic_numbers[I])
        # increase number of grid points is requested
        Nr *= radial_grid_factor
        rm = 0.5*slater_radii[atomic_names[I]]

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
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
        # distance between grid points and atom i
        dist = np.zeros((Nr,Nang, Nat))
        for i in range(0, Nat):
            dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

        # evaluate f1 and f2 on the grid
        f1I = wr * f1(x,y,z)
        f2I = wr * f2(x,y,z)

        if cusps_separate:
            ####
            # subtract spherical averages around each atom
            #
            #  f_i   = w (r) [ f(r) - sum  g (|r-Rj|) ]
            #           i                j  j
            for j in range(0, Nat):
                # rJ = |r-Rj]
                rJ = dist[:,:,j]
                # transform to z-coordinates
                rmJ = 0.5*slater_radii[atomic_names[j]]
                xrJ = (rJ-rmJ)/(rJ+rmJ)
                zrJ = np.arccos(xrJ) / np.pi

                #  for f1
                gJ_spline_1 = radial_avg_1[j]
                f1I -= wr * interpolate.splev(zrJ, gJ_spline_1, der=0, ext=0)
                #  for f2
                gJ_spline_2 = radial_avg_2[j]
                f2I -= wr * interpolate.splev(zrJ, gJ_spline_2, der=0, ext=0)

            ####

        # Perform spherical wave decomposition of fI
        #     (I)              m=+l    (I)
        #    f     = sum    sum       f  (r)  Y (th,ph)
        #               l=0    m=-l     l,m    l,m
        # and compute the derivatives
        #    d(fI_{l,m})/dr
        # The radial functions fI_{l,m} and dfI_{l,m}/dr are splined
        # in z-coordinates for later use.

        # List of spherical wave components f^(I)_{l,m},
        # the spline for f^(I)_{l,m} is stored in   radial_functions_f[I][(l,m)]
        radial_functions_f1.append( {} )
        radial_functions_f2.append( {} )
        # and the spline for df/dr^(I)_{l,m} is stored in   radial_functions_dfdr[I][(l,m)]
        radial_functions_df1dr.append( {} )
        radial_functions_df2dr.append( {} )
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,l,m in sph_it:
            wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
            # project f1I onto spherical harmonic Y^*_{l,m}
            f1I_lm = 4.0*np.pi * np.sum(f1I*wYlm, axis=-1)
            f2I_lm = 4.0*np.pi * np.sum(f2I*wYlm, axis=-1)

            # spline angular components of f1
            spline_lm_real = interpolate.splrep(zr, f1I_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, f1I_lm.imag, s=0)
            radial_functions_f1[-1][(l,m)] = spline_lm_real, spline_lm_imag

            # spline angular components of f2
            spline_lm_real = interpolate.splrep(zr, f2I_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, f2I_lm.imag, s=0)
            radial_functions_f2[-1][(l,m)] = spline_lm_real, spline_lm_imag

            
            # variable transformation
            omx = 1-xr
            opx = 1+xr

            # d f   d f
            # --- = --- (-1/(2*pi*rm) (1-x)^(3/2) (1+x)(-1/2))
            # d r   d z
            #          d f
            #     = c1 ---
            #          d z

            c0 = 0.0*omx
            c1 = -1.0/(2.0*np.pi*rm) * omx**(1.5) * opx**(-0.5)
            c2 = 0.0*omx
            
            df1Idr_lm = radial_difop(c0, c1, c2, f1I_lm)
            df2Idr_lm = radial_difop(c0, c1, c2, f2I_lm)

            # spline angular components of radial gradient df1/dr
            spline_lm_real = interpolate.splrep(zr, df1Idr_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, df1Idr_lm.imag, s=0)
            radial_functions_df1dr[-1][(l,m)] = spline_lm_real, spline_lm_imag

            # spline angular components of radial gradient df2/dr
            spline_lm_real = interpolate.splrep(zr, df2Idr_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, df2Idr_lm.imag, s=0)
            radial_functions_df2dr[-1][(l,m)] = spline_lm_real, spline_lm_imag

            
            if m == -(Lmax-1)/2:
                break

    # 4) Define function for evalulating (grad f1).(grad f2)
    def sigma12_func(x,y,z):
        """
        function for evaluating sigma = (grad f1).(grad f2)

        In spherical coordinates the gradient is

        grad f = e  g  + e   g   + e   g 
                  r  r    th  th    ph  ph

        where e_r, e_th and e_ph are unit vectors and 

            g  = df/dr   ;   g  = 1/r df/dth   ;   g  = 1/(r sin(th)) df/dph
             r                th                    ph

        For resummation over the centers we need to convert the gradients to cartesian
        coordinates so that all unit vectors point into the same direction.
        """
        grad1_x = 0.0j*x
        grad1_y = 0.0j*x
        grad1_z = 0.0j*x

        grad2_x = 0.0j*x
        grad2_y = 0.0j*x
        grad2_z = 0.0j*x

        # sum over centers
        #   grad_x = sum_I grad_x^(I)
        #   etc.
        # 
        for I in range(0, Nat):
            xI = x - atomic_coordinates[0,I]
            yI = y - atomic_coordinates[1,I]
            zI = z - atomic_coordinates[2,I]
            # spherical coordinates
            rI,thI,phI = cartesian2spherical((xI,yI,zI))
            #
            # unit vectors in cartesian cordinates
            stI = np.sin(thI)
            ctI = np.cos(thI)
            spI = np.sin(phI)
            cpI = np.cos(phI)
            #
            unit_r  = [stI*cpI, stI*spI, ctI]
            unit_th = [ctI*cpI, ctI*spI, -stI]
            unit_ph = [-spI   , cpI    , 0*stI]

            cosI = ctI
            sinI = stI
            s2 = sinI**2

            #
            sph_it = spherical_harmonics_it(thI,phI)
            # inverse transformation rI -> zr
            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:

                ############ for f1 ######################################################
                # interpolate fI_{l,m}(x,y,z)
                spline_lm_real, spline_lm_imag = radial_functions_f1[I][(l,m)]
                fI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                # interpolate dfI_{l,m}/dr(x,y,z)
                spline_lm_real, spline_lm_imag = radial_functions_df1dr[I][(l,m)]
                dfIdr_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                # compute g_r = df/dr(x,y,z) = sum_{l,m} dfIdr_{l,m}(r) Y_{l,m}(th,ph)
                grad_r = dfIdr_lm * Ylm

                if cusps_separate:
                    #####
                    if l == 0 and m == 0:
                        # add back the radial gradient due to spherical
                        # averaged part g_i that has been subtracted previously:
                        #   __     __                            __
                        #   \/ f = \/ [ f - sum_i g_i ]  + sum_i \/ g_i
                        #   
                        grad_gI_r = interpolate.splev(zr, radial_deriv_avg_1[I], der=0, ext=0)
                        grad_r += grad_gI_r

                    #####

                # compute g_th = 1/r df/dth = sum_{l,m} 1/r f_{l,m}(r) d/dth Y_{l,m}(th,ph)
                #  = sum_{l,m} [ 1/r f_{l,m} m cot(th)
                #                + Heaviside(-l+1 <= m) 1/r f_{l,m-1} sqrt((l-m+1)(l+m)) exp(-i ph) ] Y_{l,m}
                #
                # The formula for d(Y_lm)/dth is taken from http://functions.wolfram.com/05.10.20.0001.01
                #
                # To avoid dividing by 0, we first multiply all terms by sin(th) and at the end
                # divide those elements of the array where sin(th) != 0 by sin(th) again.
                # I am not sure if removing the singularity at th=0,pi in this way is correct.
                #
                grad_th = 1.0/rI * fI_lm * m * cosI * Ylm
                if -l+1 <= m:
                    # interpolate fI_{l,m-1}
                    spline_lm_real, spline_lm_imag = radial_functions_f1[I][(l,m-1)]
                    fI_lmm1 = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                            + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                    grad_th += 1.0/rI * fI_lmm1 * np.sqrt((l-m+1)*(l+m)) * np.exp(-1.0j*phI) * sinI * Ylm
                # Those elements where division by sin(th) would lead to NaN are set to zero. 
                grad_th[s2 > 0.0] /= sinI[s2 > 0.0]
                grad_th[s2 <= 0.0] = 0.0
                    
                # compute g_ph = 1/(r*sin(th)) df/dph = sum_{l,m} (i m) 1/(r*sin(th)) f_{l,m} Y_{l,m}
                grad_ph = 1.0/rI * fI_lm * (1.0j*m) * Ylm
                # Those elements where division by sin(th) would lead to NaN are set to zero. 
                grad_ph[s2 > 0.0] /= sinI[s2 > 0.0]
                grad_ph[s2 <= 0.0] = 0.0
                
                # We cannot add the gradients belonging to different spherical coordinate systems
                # since the unit vectors point in different directions. 

                grad1_x += unit_r[0] * grad_r + unit_th[0] * grad_th + unit_ph[0] * grad_ph
                grad1_y += unit_r[1] * grad_r + unit_th[1] * grad_th + unit_ph[1] * grad_ph
                grad1_z += unit_r[2] * grad_r + unit_th[2] * grad_th + unit_ph[2] * grad_ph

                ############ for f2 ######################################################
                # This is a direct copy of the above code with f1 replaced by f2
                ##########################################################################
                # interpolate fI_{l,m}(x,y,z)
                spline_lm_real, spline_lm_imag = radial_functions_f2[I][(l,m)]
                fI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                # interpolate dfI_{l,m}/dr(x,y,z)
                spline_lm_real, spline_lm_imag = radial_functions_df2dr[I][(l,m)]
                dfIdr_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                # compute g_r = df/dr(x,y,z) = sum_{l,m} dfIdr_{l,m}(r) Y_{l,m}(th,ph)
                grad_r = dfIdr_lm * Ylm


                if cusps_separate:
                    #####
                    if l == 0 and m == 0:
                        # add back the radial gradient due to spherical
                        # averaged part g_i that has been subtracted previously:
                        #   __     __                            __
                        #   \/ f = \/ [ f - sum_i g_i ]  + sum_i \/ g_i
                        #   
                        grad_gI_r = interpolate.splev(zr, radial_deriv_avg_2[I], der=0, ext=0)
                        grad_r += grad_gI_r

                    #####

                # compute g_th = 1/r df/dth = sum_{l,m} 1/r f_{l,m}(r) d/dth Y_{l,m}(th,ph)
                #  = sum_{l,m} [ 1/r f_{l,m} m cot(th)
                #                + Heaviside(-l+1 <= m) 1/r f_{l,m-1} sqrt((l-m+1)(l+m)) exp(-i ph) ] Y_{l,m}
                #
                # The formula for d(Y_lm)/dth is taken from http://functions.wolfram.com/05.10.20.0001.01
                #
                # To avoid dividing by 0, we first multiply all terms by sin(th) and at the end
                # divide those elements of the array where sin(th) != 0 by sin(th) again.
                # I am not sure if removing the singularity at th=0,pi in this way is correct.
                #
                grad_th = 1.0/rI * fI_lm * m * cosI * Ylm
                if -l+1 <= m:
                    # interpolate fI_{l,m-1}
                    spline_lm_real, spline_lm_imag = radial_functions_f2[I][(l,m-1)]
                    fI_lmm1 = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                            + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)

                    grad_th += 1.0/rI * fI_lmm1 * np.sqrt((l-m+1)*(l+m)) * np.exp(-1.0j*phI) * sinI * Ylm
                # Those elements where division by sin(th) would lead to NaN are set to zero. 
                grad_th[s2 > 0.0] /= sinI[s2 > 0.0]
                grad_th[s2 <= 0.0] = 0.0
                    
                # compute g_ph = 1/(r*sin(th)) df/dph = sum_{l,m} (i m) 1/(r*sin(th)) f_{l,m} Y_{l,m}
                grad_ph = 1.0/rI * fI_lm * (1.0j*m) * Ylm
                # Those elements where division by sin(th) would lead to NaN are set to zero. 
                grad_ph[s2 > 0.0] /= sinI[s2 > 0.0]
                grad_ph[s2 <= 0.0] = 0.0
                
                # We cannot add the gradients belonging to different spherical coordinate systems
                # since the unit vectors point in different directions. 

                grad2_x += unit_r[0] * grad_r + unit_th[0] * grad_th + unit_ph[0] * grad_ph
                grad2_y += unit_r[1] * grad_r + unit_th[1] * grad_th + unit_ph[1] * grad_ph
                grad2_z += unit_r[2] * grad_r + unit_th[2] * grad_th + unit_ph[2] * grad_ph

                
                if m == -(Lmax-1)/2:
                    break

        # Now we compute the product of the gradient vectors
        sigma12 = grad1_x*grad2_x+ grad1_y*grad2_y + grad1_z*grad2_z
        
        return sigma12.real
    
    return sigma12_func


def spherical_average_func(atom, f,
                           lebedev_order=23, radial_grid_factor=1):
    """
    create a function avg(r) that evaluates the spherical average
    of f(x,y,z) around the center defined by atom (Zat,(X,Y,Z)).

                          /             /
        avg(r) = 1/(4 pi) | sin(th) dth | dph  f(r-Rc)
                          /             /

    Parameters
    ----------
    atom              :  tuple (Zat,(xc,yc,zc)) where Zat is the atomic
                         number used to select the integration grid
                         and Rc=(xc,yc,zc) are the coordinates of the center
                         which is taken as the origin
    f                 :  callable f(x,y,z)

    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor


    Returns
    -------
    avg                : callable avg(r) that evaluates the spherical
                         average of f around the atom
    """
    # origin
    Zat,(xc,yc,zc) = atom
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)

    # radial grid
    Nr = number_of_radial_points(Zat)
    # increase number of grid points is requested
    Nr *= radial_grid_factor
    rm = 0.5*slater_radii[atom_names[Zat-1]]

    k = np.array(range(1,Nr+1))
    # grid points on interval [-1,1]
    zr = k/(Nr+1.0)
    xr = np.cos(zr * np.pi)
    # radial grid points on interval [0,infinity]
    r = rm * (1+xr)/(1-xr)

    # cartesian coordinates of grid
    x = (outerN(r, sc) + xc)
    y = (outerN(r, ss) + yc)
    z = (outerN(r, c ) + zc)
    #
    Npts = Nr*Nang

    # evaluate function on the grid
    fI = f(x,y,z)
    
    sph_it = spherical_harmonics_it(th,ph)
    for Ylm,l,m in sph_it:
        wYlm = outerN(np.ones(Nr), angular_weights*Ylm.conjugate())
        # We only need the l=0,m=0 term
        if l == 0 and m == 0:
            fI_00 = 4.0*np.pi * np.sum(fI*wYlm, axis=-1)
            break

    spline_00_real = interpolate.splrep(zr, fI_00.real, s=0)
    spline_00_imag = interpolate.splrep(zr, fI_00.imag, s=0)
        
    def avg_func(r):
        xr = (r-rm)/(r+rm)
        zr = np.arccos(xr) / np.pi

        avg00 =        interpolate.splev(zr, spline_00_real, der=0, ext=0) \
               + 1.0j*interpolate.splev(zr, spline_00_imag, der=0, ext=0)
        # avg00 is the projection of f onto Y_00 = 1/sqrt(4*pi)
        # To average over the solid angle, another factor of
        # 1/sqrt(4*pi) is missing.
        avg = avg00 * np.sqrt(1.0/(4.0*np.pi))
        
        return avg.real

    return avg_func


def radial_component_func(atom, f, l, m,  
                          lebedev_order=23, radial_grid_factor=1):
    """
    create a function f_{l,m}(r) that evaluates the projection 
    of f(x,y,z) onto the real spherical harmonic Y_{l,m}(th,ph).
    The origin of the coordinate system and the integration grid
    is defined by atom (Zat,(X,Y,Z)).

                  / /
        f  (r) =  | | Y   (th,ph)  f(r,th,ph)  sin(th) dth dph  
         l,m      / /  l,m           

    Parameters
    ----------
    atom              :  tuple (Zat,(xc,yc,zc)) where Zat is the atomic
                         number used to select the integration grid
                         and xc,yc,zc are the coordinates of the center
                         which is taken as the origin
    f                 :  callable f(x,y,z)
    l,m               :  integers, -l <= m <= l, angular quantum numbers

    Optional
    --------
    lebedev_order      : order Lmax of the Lebedev grid
    radial_grid_factor : the number of radial grid points is increased by this factor


    Returns
    -------
    f_lm               : callable f_lm(r) that evaluates the radial part of the
                         l,m-component of f(x,y,z)
    """
    Zat,(xc,yc,zc) = atom
    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)

    # radial grid
    Nr = number_of_radial_points(Zat)
    # increase number of grid points is requested
    Nr *= radial_grid_factor
    rm = 0.5*slater_radii[atom_names[Zat-1]]

    k = np.array(range(1,Nr+1))
    # grid points on interval [-1,1]
    zr = k/(Nr+1.0)
    xr = np.cos(zr * np.pi)
    # radial grid points on interval [0,infinity]
    r = rm * (1+xr)/(1-xr)

    # cartesian coordinates of grid
    x = (outerN(r, sc) + xc)
    y = (outerN(r, ss) + yc)
    z = (outerN(r, c ) + zc)
    #
    Npts = Nr*Nang

    # evaluate function on the grid
    fI = f(x,y,z)
    
    sph_it = spherical_harmonics_it(th,ph)
    for Ylm,ll,mm in sph_it:
        # We only need the l=0,m=0 term
        if ll == l and mm == m:
            # real spherical harmonics
            if m < 0:
                Ylm_real = -np.sqrt(2.0) * Ylm.imag
            elif m > 0:
                Ylm_real =  np.sqrt(2.0) * (-1)**m * Ylm.real
            else:
                # m == 0
                Ylm_real = Ylm.real

            wYlm = outerN(np.ones(Nr), angular_weights*Ylm_real)
                
            fI_lm = 4.0*np.pi * np.sum(fI*wYlm, axis=-1)
            break

    spline_lm = interpolate.splrep(zr, fI_lm, s=0)
        
    def f_lm_func(r):
        xr = (r-rm)/(r+rm)
        zr = np.arccos(xr) / np.pi

        f_lm = interpolate.splev(zr, spline_lm, der=0, ext=0)
        
        return f_lm

    return f_lm_func


def atomlist2arrays(atomlist):
    """
    convert geometry specification to numpy arrays

    Parameters
    ==========
    atomlist        :   list of tuples (Z,[x,y,z]) for each atom,
                        molecular geometry

    Returns
    =======
    atomic_numbers     : numpy array with atomic numbers
    atomic_coordinates : numpy array of shape (3,Nat) 
                         with coordinates
    """
    Nat = len(atomlist)
    atomic_numbers = np.zeros(Nat, dtype=int)
    atomic_coordinates = np.zeros((3,Nat))
    for i in range(0, Nat):
        Z,pos = atomlist[i]
        atomic_numbers[i] = Z
        atomic_coordinates[:,i] = pos
    return atomic_numbers, atomic_coordinates


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
        xr = np.cos(k/(Nr+1.0) * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
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

def join_grids(points, weights, volumes):
    """
    combine the multicenter grids into a single grid so that we get
    a quadrature rule for integration
         /
         | f(x,y,z) dV = sum  w  f(x ,y ,z ) 
         /                  i  i    i  i  i

    Parameters
    ----------
    points, weights, volumes:  return values of `multicenter_grids`

    Returns
    -------
    x,y,z     :  1d numpy arrays with cartesian coordinates of grid points
    w         :  1d numpy array with weights
    """
    # weights of quadrature rule
    w = []
    # sampling points of quadrature rule
    x,y,z = [],[],[]
    # xI,yI,zI : grid points of spherical grid around atom I
    # wI : weights of spherical grid around atom I
    # vI : volume elements of spherical grid around atom I
    for (xI,yI,zI), wI, dVI in zip(points, weights, volumes):
        x += [xI]
        y += [yI]
        z += [zI]
        # The weights are the product of the weight function
        # (from fuzzy Voronoi decomposition of space) and the volume element.
        w += [wI*dVI]
    # join arrays
    w = np.hstack(w)
    x = np.hstack(x)
    y = np.hstack(y)
    z = np.hstack(z)

    return x,y,z, w
        
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

        k = np.array(range(1,Nr+1))
        # grid points on interval [-1,1]
        zr = k/(Nr+1.0)
        xr = np.cos(zr * np.pi)
        # weights
        radial_weights = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
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
                dist[:,:,i] = np.sqrt(  (x - atomic_coordinates[0,i])**2   \
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

            spline_lm_real = interpolate.splrep(zr, fI_lm.real, s=0)
            spline_lm_imag = interpolate.splrep(zr, fI_lm.imag, s=0)
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

            rm = 0.5*slater_radii[atomic_names[I]]
            xr = (rI-rm)/(rI+rm)
            zr = np.arccos(xr) / np.pi

            for Ylm,l,m in sph_it:
                
                spline_lm_real, spline_lm_imag = radial_functions[I][(l,m)]
                # interpolate
                fI_lm = interpolate.splev(zr, spline_lm_real, der=0, ext=0) \
                        + 1.0j*interpolate.splev(zr, spline_lm_imag, der=0, ext=0)
                
                f += fI_lm*Ylm

                if m == -(Lmax-1)/2:
                    break

        return f.real
    
    return interpolation_func


###############################################################################################
#
# TESTING
#
################################################################################################

from DFTB.MolecularIntegrals.fchkfile import G09ResultsDFT
from DFTB.MolecularIntegrals import gtobasis
from DFTB.MolecularIntegrals import integrals

def radial_integration(f, rm, Nr):
    """
    Gauss-Chebyshev integration of the second kind for radial integral
          /oo    2
      I = |  dr r  f(r)
          /0
    """
    k = np.array(range(1,Nr+1))
    # grid points on interval [-1,1]
    x = np.cos(k/(Nr+1.0) * np.pi)
    # weights
    w = np.pi/(Nr+1.0) * np.sin(k/(Nr+1.0) * np.pi)**2
    g = 2 * np.sqrt(((1+x)/(1-x)**3)**3)
    w *= g
    # radial grid points on interval [0,infinity]
    r = rm * (1+x)/(1-x)
    #
    I = np.sum( rm**3 * w * f(r) )
    return I

def test_radial_integration():
    def f(r):
        return np.exp(-r**2)

    print "Nr        I"
    for Nr in range(10,100):
        print "%d        %e" % (Nr, radial_integration(f, 1.0, Nr))
    print "exact       %e" % (np.sqrt(np.pi)/4.0)


def electron_density_function(orbs_alpha, orbs_beta, nelec_alpha, nelec_beta, basis):
    
    # build density matrix
    occ_orbs_alpha = orbs_alpha[:,:nelec_alpha]
    occ_orbs_beta  = orbs_beta[:,:nelec_beta]
    Pa = np.dot(occ_orbs_alpha, occ_orbs_alpha.transpose())
    Pb = np.dot(occ_orbs_beta , occ_orbs_beta.transpose())
    P = Pa + Pb

    # overlap matrix
    S = integrals.basis_overlap(basis)
    electronic_charge = np.sum(S*P)
    
    def density(x,y,z):
        rho = 0.0*x
        for i in range(0, basis.nbfs):
            aoI = gtobasis.wavefunction(basis.exponents[i], basis.powers[:,i], basis.centers[:,i],
                                     x,y,z)
            for j in range(0, basis.nbfs):
                aoJ = gtobasis.wavefunction(basis.exponents[j], basis.powers[:,j], basis.centers[:,j],
                                            x,y,z)
                rho += P[i,j]*aoI*aoJ
        return rho

    return electronic_charge, density

def test_charge_from_numerical_integration(res):
    """
    compare the total electronic charge with the charge obtained by numerically
    integrating the electron density
    """
    exact_charge, density = electron_density_function(res.orbs_alpha, res.orbs_beta,
                                        res.nelec_alpha, res.nelec_beta, res.basis)

    print ""
    print "Lebedev order    radial grid factor         integrated charge"
    print "-------------------------------------------------------------"

    for radial_grid_factor in [1,2,3,4]:
        for lebedev_order in [11, 17, 23]:
            integrated_charge = multicenter_integration(density, res.coordinates, res.atomic_numbers,
                                                        radial_grid_factor=radial_grid_factor,
                                                        lebedev_order=lebedev_order)
            print "    %2.d                  %2.d                    %e" % (lebedev_order, radial_grid_factor, integrated_charge)

    print "        exact                                 %e" % exact_charge
    
def test_coulomb_integral_from_numerical_integration(res):
    """
    compute the classical Coulomb energy

              / / rho(1) rho(2)                    /
      I = 1/2 | | ------------- d^3r1 d^3r2  = 1/2 | V(r) rho(r) d^3r
              / /    r_12                          /

    for the total electron density rho using radial and angular grids of increasing
    resolution.

    Parameters
    ----------
    res         : instance of G09ResultsDFT with basis functions and orbital coefficients
                  
    """
    exact_charge, density = electron_density_function(res.orbs_alpha, res.orbs_beta,
                                        res.nelec_alpha, res.nelec_beta, res.basis)

    print ""
    print "Lebedev order    radial grid factor       classical Coulomb energy   "
    print "---------------------------------------------------------------------"

    for radial_grid_factor in [1,2,3,4]:
        for lebedev_order in [11, 17, 23]:
            # solve for the electrostatic potential generated by the charge density rho(r)
            elec_potential_func = multicenter_poisson(density, res.coordinates, res.atomic_numbers,
                                                      radial_grid_factor=radial_grid_factor,
                                                      lebedev_order=lebedev_order)
            def elec_energy_func(x,y,z):
                return 0.5 * density(x,y,z) * elec_potential_func(x,y,z)
            
            en_elec = multicenter_integration(elec_energy_func, res.coordinates, res.atomic_numbers,
                                              radial_grid_factor=radial_grid_factor,
                                              lebedev_order=lebedev_order)
            print "    %2.d                  %2.d                    %e" % (lebedev_order, radial_grid_factor, en_elec)
    #
    # compute classical Coulomb energy from basis functions using analytical formulae
    # for electron repulsion integrals
    #
    basis = res.basis
    norms = np.zeros(basis.nbfs)
    for i in range(0, basis.nbfs):
        norms[i] = norm(basis.exponents[i], basis.powers[:,i])
    
    # build density matrix
    occ_orbs_alpha = res.orbs_alpha[:,:res.nelec_alpha]
    occ_orbs_beta  = res.orbs_beta[:,:res.nelec_beta]
    Pa = np.dot(occ_orbs_alpha, occ_orbs_alpha.transpose())
    Pb = np.dot(occ_orbs_beta , occ_orbs_beta.transpose())
    P = Pa + Pb

    # compute
    #
    #  I = 1/2 sum_{a,b,c,d} P_{a,b} (ab|cd) P_{c,d}
    #
    I = 0.0
    n = res.basis.nbfs
    for a in range(0, n):
        for b in range(0, n):
            for c in range(0, n):
                for d in range(0, n):
                    abcd = coulomb_repulsion(basis.centers[:,a], norms[a], basis.powers[:,a], basis.exponents[a],
                                             basis.centers[:,b], norms[b], basis.powers[:,b], basis.exponents[b],
                                             basis.centers[:,c], norms[c], basis.powers[:,c], basis.exponents[c],
                                             basis.centers[:,d], norms[d], basis.powers[:,d], basis.exponents[d])
                    I += 0.5 * P[a,b] * abcd * P[c,d]

    print "          exact                               %e" % I

    
def test_h2_coulomb_exchange_integrals():
    """
    compute the Coulomb and exchange integrals

         (AA|BB), (AB|AB)

    for the hydrogen molecule H2 for a bond length of 1.40 bohr. The orbitals A and B
    are 1s hydrogen orbitals centered on either proton. The numbers from table I in 
    ref. [2] are reproduced. 
    """
    def hydrogen_1s(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        psi = 2.0/np.sqrt(4.0*np.pi) * np.exp(-r)
        return psi
    # distance between the protons (in bohr)
    dist = 1.40
    #
    def rhoAB(x,y,z):
        # density of 1s orbital on nuclei A and B
        psiA = hydrogen_1s(x,y,z)
        psiB = hydrogen_1s(x,y,z-dist)
        return psiA*psiB
    def rhoAA(x,y,z):
        # density of 1s orbital on nucleus A
        psiA = hydrogen_1s(x,y,z)
        return psiA**2
    def rhoBB(x,y,z):
        # density of 1s orbital on nucleus B
        psiB = hydrogen_1s(x,y,z-dist)
        return psiB**2

    rfac=2
    Lmax=23
    
    atomic_coordinates = np.zeros((3,2))
    atomic_coordinates[2,1] = dist
    atomic_numbers = np.array([1,1])

    normA = multicenter_integration(rhoAA, atomic_coordinates, atomic_numbers,
                                    radial_grid_factor=rfac, lebedev_order=Lmax)
    print "normalization (a|a)= %s" % normA

    print "hydrogen molecule H2"
    print ""
    print "         Integration Mesh                                                          "
    print "Lebedev order    radial grid factor       Coulomb (aa|bb)      Exchange (ab|ab)    "
    print "-----------------------------------------------------------------------------------"

    for radial_grid_factor in [1,2,3,4]:
        for lebedev_order in [11, 17, 23]:
            # Coulomb integral
            Vbb = multicenter_poisson(rhoBB, atomic_coordinates, atomic_numbers,
                                      radial_grid_factor=radial_grid_factor,
                                      lebedev_order=lebedev_order)
            def Iaabb_integrand(x,y,z):
                return rhoAA(x,y,z) * Vbb(x,y,z)
            Iaabb = multicenter_integration(Iaabb_integrand, atomic_coordinates, atomic_numbers,
                                            radial_grid_factor=radial_grid_factor,
                                            lebedev_order=lebedev_order)

            # exchange integral
            Vab = multicenter_poisson(rhoAB, atomic_coordinates, atomic_numbers,
                                      radial_grid_factor=radial_grid_factor,
                                      lebedev_order=lebedev_order)
            def Iabab_integrand(x,y,z):
                return rhoAB(x,y,z) * Vab(x,y,z)
            Iabab = multicenter_integration(Iabab_integrand, atomic_coordinates, atomic_numbers,
                                            radial_grid_factor=radial_grid_factor,
                                            lebedev_order=lebedev_order)

            print "    %2.d                  %2.d                  %.8f         %.8f" % (lebedev_order, radial_grid_factor, Iaabb, Iabab)
    print "         exact                              %.8f         %.8f" % (0.50352093, 0.32329114)

    
def test_hydrogen_1s_electrostatic_potential():
    """ 
    plot the electrostatic potential generated by an electron in a hydrogen 1s orbital

       psi(r) = 2/sqrt(4 pi) exp(-r)

    The exact result is

              1 - exp(-2 r) (1+r)
       V(r) = -------------------
                     r
    """
    def rho_1s(x,y,z):
        """electron density of 1s hydrogen electron"""
        r = np.sqrt(x*x+y*y+z*z)
        rho = 1.0/np.pi * np.exp(-2*r)
        return rho
    
    atomic_coordinates = np.zeros((3,1))
    atomic_numbers = np.array([1])

    rfac=2
    Lmax=23
    
    # Coulomb integral
    V1s = multicenter_poisson(rho_1s, atomic_coordinates, atomic_numbers,
                              radial_grid_factor=rfac, lebedev_order=Lmax)

    import matplotlib.pyplot as plt
    plt.xlabel("r / bohr", fontsize=17)
    plt.ylabel("electrostatic potential $V_{1s}(r)$ / a.u.", fontsize=17)

    r = np.linspace(1.0e-10, 10.0, 200)
    
    plt.plot(r, V1s(r,0*r,0*r), "x", label="numerical solution")
    plt.plot(r, (1.0-np.exp(-2*r)*(1+r))/r, lw=2, ls="-.", label="exact solution")

    plt.legend()
    plt.show()

def test_laplacian_hydrogen_1s():
    """
    compute the Laplacian for the 1s hydrogen wavefunction

       psi(r) = 2/sqrt(4 pi) exp(-r)

    The exact result is

       __2                             2
       \/  psi(r) = 2/sqrt(4 pi) (1 - --- ) exp(-r) 
                                       r
    """
    def psi_1s(x,y,z):
        """wavefunction of 1s hydrogen electron"""
        r = np.sqrt(x*x+y*y+z*z)
        psi = 1.0/np.sqrt(np.pi) * np.exp(-r)
        return psi
    
    atomic_coordinates = np.zeros((3,1))
    atomic_numbers = np.array([1])

    rfac=2
    Lmax=23
    
    # Laplacian
    lap_1s = multicenter_laplacian(psi_1s, atomic_coordinates, atomic_numbers,
                                   radial_grid_factor=rfac, lebedev_order=Lmax)

    import matplotlib.pyplot as plt
    plt.xlabel("r / bohr", fontsize=17)
    plt.ylabel(r"Laplacian $\nabla^2 \psi$ / a.u.", fontsize=17)

    r = np.linspace(1.0e-5, 10.0, 100000)
    
    plt.plot(r, lap_1s(r,0*r,0*r), label="numerical solution")
    plt.plot(r, 1/np.sqrt(np.pi) * (1.0-2.0/r)*np.exp(-r), lw=2, ls="-.", label="exact solution")

    plt.legend()
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
#    atomlist = [(1, (0,0,0))]


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

    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    # set radial and angular resolution of grid
    rfac=10
    Lmax=41
    
    # comupte Laplacian on the grid
    lap_sigma = multicenter_laplacian(psi_sigma, atomic_coordinates, atomic_numbers,
                                      cusps_separate=False,
                                      radial_grid_factor=rfac, lebedev_order=Lmax)
    lap_sigma_cusps_separate = multicenter_laplacian(psi_sigma, atomic_coordinates, atomic_numbers,
                                      cusps_separate=True,
                                      radial_grid_factor=rfac, lebedev_order=Lmax)

    # plot Laplacian along the H-H bond (the z-axis)
    import matplotlib.pyplot as plt
    plt.xlabel("z / bohr", fontsize=17)

    r = np.linspace(-10.0, 10.0, 100000)
    x, y, z = 0*r, 0*r, r

    plt.plot(r, psi_sigma(x,y,z), label=r"wavefunction $\psi_{\sigma}$")
    plt.plot(r, lap_sigma(x,y,z), label=r"$\nabla^2 \psi$, numerical solution")
    plt.plot(r, lap_sigma_cusps_separate(x,y,z), ls="-.", label=r"$\nabla^2 \psi$, numerical solution (cusps separate)")
    plt.plot(r, lap_sigma_exact(x,y,z), ls="--", label=r"$\nabla^2 \psi$, analytical solution")

    plt.legend()
    plt.show()

    
    
def test_multicenter_operations():
    """
    test addition of functions on the grid
    """
    atomlist = [(1, (0.0, 0.0, 0.0))]
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    
    def f(x,y,z):
        return np.sin(x+y+z)

    def g(x,y,z):
        return np.exp(-0.5 * (x**2 + y**2 + z**2))

    a = 0.3
    b = 0.7
    
    def h_test(x,y,z):
        return a*f(x,y,z) + b*g(x,y,z)

    rfac=3
    Lmax=23
    # h() = a*f() + b*g()
    h = multicenter_operation([f, g], lambda fs: a*fs[0]+b*fs[1], 
                              atomic_coordinates, atomic_numbers,
                              radial_grid_factor=rfac,
                              lebedev_order=Lmax)

    # check that h and h_test agree
    Npts = 100
    x = np.random.rand(Npts)
    y = np.random.rand(Npts)
    z = np.random.rand(Npts)

    err = la.norm( h(x,y,z) - h_test(x,y,z) )
    print "error |h - h_test|= %e" % err
    assert err < 1.0e-3

def test1_multicenter_gradient2():

    # Define multicenter grid, any grid centered around the
    # origin should work
    atomlist = [(6, (0,-0.25,0.01)),
                (6, (0,+0.25,0.01))]
    #atomlist = [(1, (0,0,0))]
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # resolution of grid
    rfac=10
    Lmax=41

    # We need to know the square of the gradient of this function
    def f(x,y,z):
        return x*y*np.sin(z)*np.exp(-(x**2+y**2+z**2))
    
    #        __   2
    # exact (\/ f)   computed in cartesian coordinates
    def sigma_test(x,y,z):
        expmr2 = np.exp(-(x**2+y**2+z**2))

        sinz = np.sin(z)
        cosz = np.cos(z)
        dfdx = (1 - 2*x**2)*y*sinz*expmr2
        dfdy = (1 - 2*y**2)*x*sinz*expmr2
        dfdz = (cosz - 2*z*sinz)*x*y*expmr2

        sgma = dfdx**2+dfdy**2+dfdz**2
        return sgma

    # compute square of gradient on multicenter grid
    sigma = multicenter_gradient2(f,
                                  atomic_coordinates, atomic_numbers,
                                  lebedev_order=Lmax, radial_grid_factor=rfac)
    # compute product of gradient with itself, should give the same result
    sigma12 = multicenter_gradient_product(f, f, 
                                           atomic_coordinates, atomic_numbers,
                                           lebedev_order=Lmax, radial_grid_factor=rfac)

    ###
    # compare visually
    import matplotlib.pyplot as plt
    Npts = 500
    r = np.linspace(-3,3,Npts)
    zero = np.zeros(Npts)

    plt.xlabel("r")
    plt.ylabel(r"$\sigma = \vert \nabla f \vert^2$")
    
    # cut along diagonal
    l, = plt.plot(r, sigma(r,r,r),        label="$\sigma(r,r,r)$ (approx.)")
    plt.plot(r, sigma_test(r,r,r),   "o", label="$\sigma(r,r,r)$ (exact)", mfc='none', color=l.get_color())
    plt.plot(r, sigma12(r,r,r),      "x", label="$\sigma12(r,r,r)$ (approx.)", mfc='none', color=l.get_color())

    """
    # cut along x-axis
    l, = plt.plot(r, sigma(r,0,0)+1,      label="$\sigma(r,0,0)+1$ (approx.)")
    plt.plot(r, sigma_test(r,0,0)+1, "o", label="$\sigma(r,0,0)+1$ (exact)", mfc='none', color=l.get_color())    
    # cut along z-axis
    l, = plt.plot(r, sigma(0,r,0)+2,      label="$\sigma(0,r,0)+2$ (approx.)")
    plt.plot(r, sigma_test(0,r,0)+2, "o", label="$\sigma(0,r,0)+2$ (exact)", mfc='none', color=l.get_color())
    # cut along z-axis
    l, = plt.plot(r, sigma(0,0,r)+3,      label="$\sigma(0,0,r)+3$ (approx.)")
    plt.plot(r, sigma_test(0,0,r)+3, "o", label="$\sigma(0,0,r)+3$ (exact)", mfc='none', color=l.get_color())
    """
    
    plt.legend()
    plt.show()
    ###
    # compare analytical and approximate sigma
    Npts = 100
    x = np.random.rand(Npts)
    y = np.random.rand(Npts)
    z = np.random.rand(Npts)

    s = sigma(x,y,z)
    s12 = sigma12(x,y,z)
    s_test = sigma_test(x,y,z)

    print "square of gradient (multicenter grid)"
    print s
    print "product of gradient with itself (multicenter grid)"
    print s12
    print "square of gradient (analytical)"
    print s_test
    print "difference s-s_test"
    print s-s_test
    print "difference s12-s_test"
    print s12-s_test
    
    err = la.norm(s-s_test)
    err12 = la.norm(s12-s_test)
    print "|(grad f)^2 (approx.) - (grad f)^2 (exact)|= %e" % err
    assert err < 5.0e-3
    print "|(grad f).(grad f) (approx.) - (grad f)^2 (exact)|= %e" % err
    assert err12 < 5.0e-3

def test2_multicenter_gradient2():
    # Define multicenter grid, any grid centered around the
    # origin should work
    atomlist = [(6, (0,-0.25,0.01)),
                (6, (0,+0.25,0.01))]
    #atomlist = [(1, (0,0,0))]
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # resolution of grid
    rfac=10
    Lmax=41

    alphas = [0.3, 0.5]
    # We need to know the square of the gradient of this function
    def f(x,y,z):
        val = 0.0*x
        for I,(Z,pos) in enumerate(atomlist):
            xI = x-pos[0]
            yI = y-pos[1]
            zI = z-pos[2]
            rI = np.sqrt(xI**2+yI**2+zI**2)

            val += np.exp(-alphas[I]*rI)

        return val
    
    #        __   2
    # exact (\/ f)   computed in cartesian coordinates
    def sigma_test(x,y,z):
        dfdx = 0.0*x
        dfdy = 0.0*x
        dfdz = 0.0*x
        
        for I,(ZatI,posI) in enumerate(atomlist):
            xI = x-posI[0]
            yI = y-posI[1]
            zI = z-posI[2]
            rI = np.sqrt(xI**2+yI**2+zI**2)
            expI = np.exp(-alphas[I]*rI) * (-alphas[I])
            dfdx += expI * xI/rI
            dfdy += expI * yI/rI
            dfdz += expI * zI/rI

        sgma = dfdx**2+dfdy**2+dfdz**2
        return sgma

    # compute square of gradient on multicenter grid
    sigma = multicenter_gradient2(f,
                                  atomic_coordinates, atomic_numbers,
                                  lebedev_order=Lmax, radial_grid_factor=rfac)
    ###
    # compare visually
    import matplotlib.pyplot as plt
    Npts = 500
    r = np.linspace(-3,3,Npts)
    zero = np.zeros(Npts)

    plt.xlabel("r")
    plt.ylabel(r"$\sigma = \vert \nabla f \vert^2$")
    
    # cut along diagonal
    l, = plt.plot(r, sigma(r,r,r),        label="$\sigma(r,r,r)$ (approx.)")
    plt.plot(r, sigma_test(r,r,r),   "o", label="$\sigma(r,r,r)$ (exact)", mfc='none', color=l.get_color())
    """
    # cut along x-axis
    l, = plt.plot(r, sigma(r,0,0)+1,      label="$\sigma(r,0,0)+1$ (approx.)")
    plt.plot(r, sigma_test(r,0,0)+1, "o", label="$\sigma(r,0,0)+1$ (exact)", mfc='none', color=l.get_color())    
    # cut along z-axis
    l, = plt.plot(r, sigma(0,r,0)+2,      label="$\sigma(0,r,0)+2$ (approx.)")
    plt.plot(r, sigma_test(0,r,0)+2, "o", label="$\sigma(0,r,0)+2$ (exact)", mfc='none', color=l.get_color())
    # cut along z-axis
    l, = plt.plot(r, sigma(0,0,r)+3,      label="$\sigma(0,0,r)+3$ (approx.)")
    plt.plot(r, sigma_test(0,0,r)+3, "o", label="$\sigma(0,0,r)+3$ (exact)", mfc='none', color=l.get_color())
    """
    
    plt.legend()
    plt.show()
    ###
    # compare analytical and approximate sigma
    Npts = 100
    x = np.random.rand(Npts)
    y = np.random.rand(Npts)
    z = np.random.rand(Npts)

    s = sigma(x,y,z)
    s_test = sigma_test(x,y,z)

    print "square of gradient (multicenter grid)"
    print s
    print "square of gradient (analytical)"
    print s_test
    print "difference s-s_test"
    print s-s_test
    
    err = la.norm(s-s_test)
    print "|(grad f)^2 (approx.) - (grad f)^2 (exact)|= %e" % err
    assert err < 5.0e-3


def test_multicenter_gradient_product():
    # Define multicenter grid, any grid centered around the
    # origin should work
    atomlist = [(6, (0,-0.25,0.01)),
                (6, (0,+0.25,0.01))]
    #atomlist = [(1, (0,0,0))]
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # resolution of grid
    rfac=10
    Lmax=41

    alphas1 = [0.3, 0.5]
    alphas2 = [0.4, 0.7]
    # We need to know the square of the gradient of these functions
    def f1(x,y,z):
        val = 0.0*x
        for I,(Z,pos) in enumerate(atomlist):
            xI = x-pos[0]
            yI = y-pos[1]
            zI = z-pos[2]
            rI = np.sqrt(xI**2+yI**2+zI**2)

            val += np.exp(-alphas1[I]*rI)

        return val

    def f2(x,y,z):
        val = 0.0*x
        for I,(Z,pos) in enumerate(atomlist):
            xI = x-pos[0]
            yI = y-pos[1]
            zI = z-pos[2]
            rI = np.sqrt(xI**2+yI**2+zI**2)

            val += np.exp(-alphas2[I]*rI)

        return val
    
    #        __      __
    # exact (\/ f1).(\/ f2)   computed in cartesian coordinates
    def sigma12_test(x,y,z):
        # cartesian gradient of f1
        df1dx = 0.0*x
        df1dy = 0.0*x
        df1dz = 0.0*x
        # cartesian gradient of f2 
        df2dx = 0.0*x
        df2dy = 0.0*x
        df2dz = 0.0*x

        for I,(ZatI,posI) in enumerate(atomlist):
            xI = x-posI[0]
            yI = y-posI[1]
            zI = z-posI[2]
            rI = np.sqrt(xI**2+yI**2+zI**2)
            # for f1
            exp1I = np.exp(-alphas1[I]*rI) * (-alphas1[I])
            df1dx += exp1I * xI/rI
            df1dy += exp1I * yI/rI
            df1dz += exp1I * zI/rI
            # for f2
            exp2I = np.exp(-alphas2[I]*rI) * (-alphas2[I])
            df2dx += exp2I * xI/rI
            df2dy += exp2I * yI/rI
            df2dz += exp2I * zI/rI

        sgma12 = df1dx*df2dx + df1dy*df2dy + df1dz*df2dz
        return sgma12

    # compute square of gradient on multicenter grid
    sigma12 = multicenter_gradient_product(f1, f2, 
                                  atomic_coordinates, atomic_numbers,
                                  lebedev_order=Lmax, radial_grid_factor=rfac)
    ###
    # compare visually
    import matplotlib.pyplot as plt
    Npts = 500
    r = np.linspace(-3,3,Npts)
    zero = np.zeros(Npts)

    plt.xlabel("r")
    plt.ylabel(r"$\sigma12 = (\nabla f_1)(\nabla f_2)$")
    
    # cut along diagonal
    l, = plt.plot(r, sigma12(r,r,r),        label="$\sigma12(r,r,r)$ (approx.)")
    plt.plot(r, sigma12_test(r,r,r),   "o", label="$\sigma12(r,r,r)$ (exact)", mfc='none', color=l.get_color())
    
    plt.legend()
    plt.show()
    ###
    # compare analytical and approximate sigma
    Npts = 100
    x = np.random.rand(Npts)
    y = np.random.rand(Npts)
    z = np.random.rand(Npts)

    s12 = sigma12(x,y,z)
    s12_test = sigma12_test(x,y,z)

    print "product of gradients (multicenter grid)"
    print s12
    print "product of gradients (analytical)"
    print s12_test
    print "difference s12-s12_test"
    print s12-s12_test
    
    err12 = la.norm(s12-s12_test)
    print "|(grad f1).(grad f2) (approx.) - (grad f1).(grad f2) (exact)|= %e" % err12
    assert err12 < 1.0e-2

    
def test_inhomogeneous_schroedinger():
    """
    check that we can solve
             __2
       - 1/2 \/ phi(r) + (V(r) - E) phi(r) = S(r)

    For a given phi(r) it is easy to calculate the right-hand side S(r).
    Then we use Becke's scheme for solving the equation for phi(r)
    given S(r) and check that we recover the original phi(r). 
    """
    E = 0.0
    # We use a spherically symmetric potential, so that we get an exact
    # solution even with the spherical averaging in the solver.
    def potential(x,y,z):
        r = np.sqrt(x**2+y**2+z**2)
        return -1.0/r

    # RHS computed for phi(r) = exp(-r^2)
    def source(x,y,z):
        r = np.sqrt(x**2+y**2+z**2)
        S = (3-2*r**2)*np.exp(-r**2) + (potential(x,y,z) - E)*np.exp(-r**2)
        return S

    # expected solution 
    def phi_test(x,y,z):
        return np.exp(-r**2)

    # define the multicenter grid
    atomlist = [(1, (0.0, 0.0, 0.0))]
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    rfac=5
    Lmax=23

    # solve for phi
    phi = multicenter_inhomogeneous_schroedinger(potential, source,  E,
                                                 atomic_coordinates, atomic_numbers,
                                                 radial_grid_factor=rfac,
                                                 lebedev_order=Lmax)

    # compare with expected solution
    import matplotlib.pyplot as plt
    Npts = 100
    zero = np.zeros(Npts)
    r = np.linspace(0.001, 5.0, Npts)

    plt.xlabel(" r / bohr ")
    plt.ylabel(r"wavefunction $\phi$")
    l, = plt.plot(r, phi(r,zero,zero),                             label=r"$\phi(r,0,0)$ (approx.)")
    plt.plot(r, phi_test(r,zero,zero), "o", color=l.get_color(),   label=r"$\phi(r,0,0)$ (exact)")
    l, = plt.plot(r, phi(zero,r,zero)+1,                           label=r"$\phi(0,r,0)+1$ (approx.)")
    plt.plot(r, phi_test(zero,r,zero)+1, "x", color=l.get_color(), label=r"$\phi(0,r,0)+1$ (exact)")
    l, = plt.plot(r, phi(zero,zero,r)+2,                           label=r"$\phi(0,0,r)+2$ (approx.)")
    plt.plot(r, phi_test(zero,zero,r)+2, "s", color=l.get_color(), label=r"$\phi(0,0,r)+2$ (exact)")
    
    plt.legend()
    plt.show()

def test_continuum_schroedinger():
    """
    check that we can solve
             __2
       - 1/2 \/ phi(r) + (V(r) - E) phi(r) = 0

    for the Coulomb potential

        V = -1/r

    """
    # Quantum numbers of continuum wavefunction
    # energy
    E = 0.2
    # asymptotic angular momentum
    l = 1
    m = 0    
    
    # Coulomb potential
    charge = +1
    def potential(x,y,z):
        r = np.sqrt(x**2+y**2+z**2)
        return -charge/r

    # define the grid
    atomlist = [(1, (0.0, 0.0, 0.0))]
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    from DFTB.MolecularIntegrals import settings
    settings.radial_grid_factor = 40
    settings.lebedev_order = 5
    
    # solve for phi
    phi = multicenter_continuum_schroedinger(potential, E, charge, 
                                             atomic_coordinates, atomic_numbers,
                                             radial_grid_factor=settings.radial_grid_factor, 
                                             lebedev_order=settings.lebedev_order)
    # __2
    # \/ phi
    lap_phi = multicenter_laplacian(phi,
                                    atomic_coordinates, atomic_numbers,
                                    radial_grid_factor=settings.radial_grid_factor, 
                                    lebedev_order=settings.lebedev_order)
    
    # exact solution
    from DFTB.MolecularIntegrals.BasissetFreeDFT import regular_coulomb_func, residual_func
    
    phi_exact = regular_coulomb_func(E, charge, l, m, 0.0)

    # residuals
    residual = residual_func(atomlist, lambda x,y,z: phi(x,y,z,l=l,m=m), potential, E)
    residual_exact = residual_func(atomlist, phi_exact, potential, E)
    
    # compare with expected solution
    import matplotlib.pyplot as plt
    Npts = 5000
    r = np.linspace(0.001, 50.0, Npts)

    #
    x = 0*r
    y = 0*r
    z = r
    
    plt.xlabel(" z / bohr ")
    # wavefunctions
    plt.plot(r, phi_exact(x,y,z), label=r"$\phi$ (exact)")
    plt.plot(r, phi(x,y,z,l=l,m=m), ls="-.", label=r"$\phi$ (approx.)")
    # residuals
    plt.plot(r, residual_exact(x,y,z), ls="-",  label=r"$(H-E)\phi$ (exact)")
    plt.plot(r, residual(x,y,z),       ls="-.", label=r"$(H-E)\phi$")

    # kinetic energy T phi
    plt.plot(r, -0.5 * lap_phi(x,y,z), ls="--", label=r"$-\frac{1}{2} \nabla^2 \phi$")
    # potential energy (V-E)phi
    plt.plot(r, (potential(x,y,z)-E) * phi(x,y,z), ls="--", label=r"$(V-E) \phi$")
    # residual (T+V-E)phi
    plt.plot(r, -0.5 * lap_phi(x,y,z) + (potential(x,y,z)-E) * phi(x,y,z), ls="--", label=r"$(T+V-E) \phi$")

    
    plt.legend()
    plt.show()


def test_spherical_average():
    # resolution of grid
    rfac=10
    Lmax=41
    
    def f(x,y,z):
        return x**2+y**2+z**2 + x + y + z
    def avg_exact(r):
        return r**2
    
    atom = (1, (0.0, 0.0, 0.0))
    
    avg = spherical_average_func(atom, f,
                                 lebedev_order=Lmax, radial_grid_factor=rfac)

    # compare avg() and avg_exact() for some points
    r = np.linspace(0.0, 20.0, 1000)
    err = la.norm(avg(r) - avg_exact(r))
    print "|avg(numerical)-avg(analytical)|= %e" % err
    
    import matplotlib.pyplot as plt
    plt.plot(r, avg(r),       label="numerical spherical average")
    plt.plot(r, avg_exact(r), label="analytical spherical average")
    plt.legend()
    
    plt.show()

def test_radial_component_func():
    # resolution of grid
    rfac=10
    Lmax=41
    
    def f(x,y,z):
        return x**2+y**2+z**2 + x + y + z
    def avg_exact(r):
        return r**2
    
    atom = (1, (0.0, 0.0, 0.0))

    # projecting onto Y_{0,0} differs from averaging over the solid angle
    # by a factor 1/sqrt(4 pi)
    f_00 = radial_component_func(atom, f, 0,0, 
                                lebedev_order=Lmax, radial_grid_factor=rfac)

    # compare avg() and avg_exact() for some points
    r = np.linspace(0.0, 20.0, 1000)
    err = la.norm(f_00(r)/np.sqrt(4.0*np.pi) - avg_exact(r))
    print "|avg(numerical)-avg(analytical)|= %e" % err
    
    import matplotlib.pyplot as plt
    plt.plot(r, f_00(r)/np.sqrt(4.0*np.pi), label="numerical spherical average")
    plt.plot(r, avg_exact(r), label="analytical spherical average")
    plt.legend()
    
    plt.show()


def plot_water_grid():
    # experimental geometry of water
    #  r(OH) = 0.958 Ang, angle(H-O-H) = 104.4776 degrees
    atomlist = [
        (8, (0.000000000000000,  0.000000000000000, -0.222540557483415)),
        (1, (0.000000000000000, +1.431214118579765,  0.886071388908105)),
        (1, (0.000000000000000, -1.431214118579765,  0.886071388908105))]

    grid_points, grid_weights, grid_volumes = multicenter_grids(atomlist, kmax=20,
                                                                lebedev_order=41, radial_grid_factor=3)

    import matplotlib.pyplot as plt
    # show cut along x=0
    eps = 0.2
    for i,((xi,yi,zi),wi) in enumerate(zip(grid_points, grid_weights)):
        # select points close to the yz-plane
        y2d = yi[(-eps <= xi) & (xi <= eps)]
        z2d = zi[(-eps <= xi) & (xi <= eps)]
        w2d = wi[(-eps <= xi) & (xi <= eps)]

        # show position of atom, on which the grid is centered
        Zat, pos = atomlist[i]
        plt.plot([pos[1]],[pos[2]], "o")
        plt.text(pos[1],pos[2], atom_names[Zat-1].upper(),
                 fontsize=20, 
                 horizontalalignment='center', verticalalignment='center')
        # show grid points, the size of the marke is proportional to the weight
        plt.scatter(y2d,z2d,s=w2d*100, alpha=0.5)

    plt.xlim((-2.5, 2.5))
    plt.ylim((-2.0, 3.0))

    plt.gca().axis("off")
    plt.show()
    
def plot_1d_weight_functions(kmax=12):
    """
    divide the line into two fuzzy cells and plot the weight functions
    """
    # for nuclear weight functions
    def s(mu, k=kmax):
        f = mu
        for ik in range(0, k):
            f = 1.5 * f -0.5 * f**3
        return 0.5*(1-f)

    # positions of atom i and j
    Ri = 3.0
    Rj = 6.0
    # 
    Rij = abs(Rj-Ri)
    r = np.linspace(0.0, 10.0, 1000)
    ri = abs(r-Ri)
    rj = abs(r-Rj)
    mu_ij = (ri-rj)/Rij
    mu_ji = -mu_ij
    
    # P_i = prod_{j != i} s(mu_ij)
    Pi = s(mu_ij)
    Pj = s(mu_ji)

    # w_i = Pi/(sum_j Pj)
    wi = Pi/(Pi+Pj)
    wj = Pj/(Pi+Pj)

    import matplotlib.pyplot as plt
    # atom i
    plt.plot([Ri], [0.5], "o", color="blue")
    plt.plot(r, wi, lw=2, color="blue")
    plt.text(Ri, 0.9, r"$w_i(r)$", color="blue", fontsize=18, horizontalalignment="center")
    plt.text(Ri, 0.4, "atom $i$", color="blue", fontsize=18, horizontalalignment="center")
    # atom j
    plt.plot([Rj], [0.5], "o", color="red")
    plt.plot(r, wj, lw=2, color="red")
    plt.text(Rj, 0.9, r"$w_j(r)$", color="red", fontsize=18, horizontalalignment="center")
    plt.text(Rj, 0.4, "atom $j$", color="red", fontsize=18, horizontalalignment="center")

    # sum of weight functions adds to 1
    plt.plot(r, wi+wj, lw=0.5, color="black")

    # hide x-axis
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels(["0", "", "", "", "", "1"])
    plt.xlabel("r", fontsize=18)
    plt.xlim((2.0, 8.0))
    plt.ylim((0.0, 1.05))
    plt.show()

def test_multicenter_residual():
    from DFTB.MolecularIntegrals.BasissetFreeDFT import effective_potential_func, residual_func

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

    # overlap between AOs
    from DFTB.MolecularIntegrals import Ints1e
    Sao = Ints1e.overlap(atomlist,
                         lambda x,y,z: psi_1s(x,y,z-R/2.0),
                         lambda x,y,z: psi_1s(x,y,z+R/2.0))
    
    def psi_sigma(x,y,z):
        """normalized LCAO wavefunction of H-H"""
        # sigma orbital is a linear combination of two
        # 1s orbitals of hydrogen
        psi_unnorm = psi_1s(x,y,z-R/2.0) + psi_1s(x,y,z+R/2.0)
        psi = psi_unnorm / np.sqrt(2.0 + 2.0*Sao)
        return psi

    def lap_sigma_exact(x,y,z):
        """                       __2
        analytical expression for \/  psi_sigma
        """
        lap_unnorm = lap_1s(x,y,z-R/2.0) + lap_1s(x,y,z+R/2.0)
        lap = lap_unnorm / np.sqrt(2.0 + 2.0*Sao)
        return lap
        
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    # set radial and angular resolution of grid
    rfac=10
    Lmax=41

    # electron density
    def rho(x,y,z):
        return 2.0*psi_sigma(x,y,z)**2

    # check normalization (integral of density should equal 2.0)
    nelec = Ints1e.integral(atomlist, rho)
    print "number of electrons  int rho(r) dr = %e" % nelec
    
    # effective potential (without nuclear attraction)
    # The keyword nuclear=False removes the nuclear potential from Veff
    potential_ee = effective_potential_func(atomlist, rho, None,
                                            nelec=2,
                                            nuclear=False)
    # effective potential (with nuclear attraction)
    potential = effective_potential_func(atomlist, rho, None,
                                         nelec=2,
                                         nuclear=True)

    # nuclear attraction potential
    def potential_nuc(x,y,z):
        # nuclear attraction
        Vnuc = 0.0*x
        for Zi,posi in atomlist:
            r = np.sqrt((x-posi[0])**2 + (y-posi[1])**2 + (z-posi[2])**2)
            Vnuc -= Zi/r
        return Vnuc
    
    def res_sigma_exact(x,y,z):
        #
        #  (T + V - en) phi
        #
        kin = -0.5 * lap_sigma_exact(x,y,z)
        # electron-electron repulsion + exchange-correlation
        Vee = potential_ee(x,y,z)
        Vnuc = potential_nuc(x,y,z)
        # (V-en)\psi
        pot = (Vnuc + Vee - en) * psi_sigma(x,y,z)
        # (T+V-en)\psi
        res = kin + pot
        return res

    # expectation value of energy
    from DFTB.MolecularIntegrals.BasissetFreeDFT import energy
    en = energy(atomlist, psi_sigma, psi_sigma, potential)
    print "energy expectation value  <psi|H|psi> = %e" % en
    
    # compute residual on the grid treating kinetic energy and nuclear
    # attraction together to minimize numerical errors
    res_sigma_together = multicenter_residual(psi_sigma,
                                     potential_ee, en,
                                     atomic_coordinates, atomic_numbers,
                                     radial_grid_factor=rfac, lebedev_order=Lmax)

    # compute residual treating kinetic energy and nuclear attraction
    # separately
    res_sigma_separate = residual_func(atomlist, psi_sigma,
                                       potential, en)

    # Quantify the deviation between the analytical residual and the
    # the numerical ones
    #  ||(H-E)\psi (analytic) - (H-E)\psi (numeric)||^2
    err2_together = multicenter_integration(lambda x,y,z:
            abs(res_sigma_together(x,y,z) - res_sigma_exact(x,y,z))**2,
                                         atomic_coordinates, atomic_numbers,
                                         radial_grid_factor=rfac, lebedev_order=Lmax)
    err2_separate = multicenter_integration(lambda x,y,z:
            abs(res_sigma_separate(x,y,z) - res_sigma_exact(x,y,z))**2,
                                         atomic_coordinates, atomic_numbers,
                                         radial_grid_factor=rfac, lebedev_order=Lmax)
    print "error of residual (separate) = %e" % np.sqrt(err2_together)
    print "error of residual (together) = %e" % np.sqrt(err2_separate)
    
    # plot residual along the H-H bond (the z-axis)
    import matplotlib.pyplot as plt
    plt.xlabel("z / bohr", fontsize=17)

    r = np.linspace(-10.0, 10.0, 100000)
    x, y, z = 0*r, 0*r, r

    plt.plot(r, psi_sigma(x,y,z), label=r"wavefunction $\psi_{\sigma}$")
    plt.plot(r, potential_nuc(x,y,z), label=r"$V_{nuc}$")
    plt.plot(r, potential_ee(x,y,z), label=r"$V_{coul} + V_{xc}$")
    plt.plot(r, res_sigma_together(x,y,z),          label=r"$(H-E) \psi$, numerical solution (new)")
    plt.plot(r, res_sigma_separate(x,y,z), ls="--", label=r"$(H-E) \psi$, numerical solution (old)")
    plt.plot(r, res_sigma_exact(x,y,z), ls="-.", label=r"$(H-E) \psi$, analytic solution")

    plt.legend()
    plt.show()


def test_join_grids():
    """
    check quadrature rule obtained by combining multicenter grids
    into a single grid
    """
    # set radial and angular resolution of grid
    rfac=10
    Lmax=21
    
    # bond length in bohr 
    R = 2.0
    # 
    atomlist = [(1, (0,0,-R/2.0)),
                (1, (0,0,+R/2.0))]
    # spherical grids around each atom
    points, weights, volumes = multicenter_grids(atomlist,
                                                 radial_grid_factor=rfac,
                                                 lebedev_order=Lmax)
    # combine multicenter grids into a single grid
    x,y,z, w = join_grids(points, weights, volumes)

    sigma = 0.5
    def gauss(r):
        return 1.0/(np.sqrt(2*np.pi)*sigma) * np.exp(-0.5*(r/sigma)**2)

    def f(x,y,z):
        return gauss(x)*gauss(y)*gauss(z)

    # The 3d Gaussian function should be normalized
    gnorm = np.sum(w * f(x,y,z))
    
    print "norm of Gaussian = %e" % gnorm
    
    
if __name__ == "__main__":
    import sys
    import os.path

    test_join_grids()
    
    if len(sys.argv) < 2:
        print "Usage: %s  <formatted checkpoint file>" % os.path.basename(sys.argv[0])
        print "  test Becke's multicenter integration scheme for electron density"
        exit(-1)

    # Formatted checkpoint files for water at the HF/SV level can be found in
    # test/h2o_hf_sv.fchk. The classical Coulomb energies should
    # be compared with table II in ref. [2].
        
    filename = sys.argv[1]
    res = G09ResultsDFT(filename)

    test_charge_from_numerical_integration(res)
    test_coulomb_integral_from_numerical_integration(res)
    test_h2_coulomb_exchange_integrals()
    test_hydrogen_1s_electrostatic_potential()

    test_laplacian_hydrogen_1s()
    test_multicenter_laplacian()
    
    test_multicenter_operations()

    test_inhomogeneous_schroedinger()
    test_continuum_schroedinger()

    test1_multicenter_gradient2()
    test2_multicenter_gradient2()
    test_multicenter_gradient_product()
    
    test_spherical_average()
    
    test_radial_component_func()

    test_multicenter_residual()
    
    plot_water_grid()
    plot_1d_weight_functions()

