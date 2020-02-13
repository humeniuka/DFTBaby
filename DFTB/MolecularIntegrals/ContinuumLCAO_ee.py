#!/usr/bin/env python
"""
Molecular continuum orbitals are approximated as a linear combination
of atomic continuum orbitals with energy E

                     atoms     oo      +li
    phi (x,y,z) = sum       sum     sum       c         u (x,y,z)        
                     i         li=0    mi=-li  i,li,mi   i,li,mi

The u's are solutions of the atomic Kohn-Sham equation:

  (T + V  - E) u        = 0                  V_i : effective potential of atom i
        i       i,li,mi


The LCAO coefficients c's are determined by minimizing the residual

                       2                                     
           ||(H-E)phi||   =  <phi|(H-E)(H-E)|phi>  

                                    *                  
                          = sum    c  c   <u |(H-E)|u >
                               i,j  i  j    i        j

                                    *
                          = sum    c  c  R
                               i,j  i  j  i,j

subject to the constraint that the c's are orthonormal. Since the u's satisfy the
atomic Kohn-Sham equations, the residual matrix can be obtained from integrals
                       +           
   R    = [ (H-E)|u > ]  (H-E)|u > 
    i,j            i            j   
                                          +
        = [ (T + V - E)|u > + (V-V )|u > ]  [ (T + V  - E)|u > + (V-V )|u > ]
                  i      i        i   i             j       j        j   j

       = <u |(V-V )(V-V )|u >
           i     i     j   j

"""

from DFTB.MolecularIntegrals import settings
from DFTB.MolecularIntegrals.Ints1e import integral, electronic_dipole, overlap
from DFTB.MolecularIntegrals.MulticenterIntegration import select_angular_grid, outerN
from DFTB.MolecularIntegrals.BasissetFreeDFT import BasissetFreeDFT, density_func, effective_potential_func, residual_func, residual_ee_func, orbital_transformation, energy_correction, inhomogeneous_schroedinger, add_two_functions, laplacian_func, regular_coulomb_func, nuclear_potential_func
from DFTB.SlaterKoster import XCFunctionals

from DFTB.Scattering.SlakoScattering import AtomicScatteringBasisSet
from DFTB.SlaterKoster.SKIntegrals import spline_effective_potential
from DFTB.BasisSets import import_pseudo_atom

from DFTB.Analyse import Cube
from DFTB import Parameters

import numpy as np
import scipy.linalg as sla

class AtomicPotential:
    def __init__(self, atom, center):
        self.center = center
        self.atom = atom
        r = atom.r
        self.veff = spline_effective_potential(atom.r[r > 0], atom.effective_potential[r > 0], atom.r0)
        # smallest radius > 0.0 for which we have a value of Veff
        imin = np.argmin(r[r > 0])
        self.rmin = r[r > 0][imin]
        self.vmin = atom.effective_potential[r > 0][imin]
        
        """
        ### DEBUG
        import matplotlib.pyplot as plt
        plt.cla()
        plt.clf()
        r = np.linspace(0.0001, 5.0, 10000)
        plt.plot(r, self.veff(r), label=r"radial $V_{eff}$")
        plt.legend()
        plt.show()
        ###
        """
        
    def __call__(self, x,y,z):
        return self.amp(x,y,z)
    
    def amp(self, x,y,z):
        x0,y0,z0 = self.center
        r = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
        v = self.veff(r)
        # join the curves Veff(r) and -Z/r continuously at r=rmin
        v[r < self.rmin] = -self.atom.Z * (1.0/r[r < self.rmin] - 1.0/self.rmin) + self.vmin
        
        return v

    
class AtomicPotentialSet:
    def __init__(self, atomlist):
        atom_data_dic = {}
        atomtypes = list(set([Zi for (Zi,posi) in atomlist]))
        atomtypes.sort()
        for Zi in atomtypes:
            confined_atom, free_atom = import_pseudo_atom(Zi)
            atom_data_dic[Zi] = free_atom

        # list of atomic potentials
        self.pots = []
        for (Zi,posi) in atomlist:
            pot = AtomicPotential(atom_data_dic[Zi], posi)
            self.pots.append(pot)            

def continuum_overlap(bfs, E, rlarge=50000.0):
    """
    Continuum orbitals cannot be normalized, but a quantity equivalent to the
    overlap matrix may be defined by integrating the product of two continuum
    wavefunctions over one wavelength at a large distance R away from the 
    origin

               /R+lambda /2pi  /pi              2    *
        S   =  |  dr     | dph | sin(th) dth   r  phi (r,th,ph)   phi (r,th,ph)
         ij    /R        /     /0                    i               j
    
    where lambda = 2pi/k = 2 pi / sqrt(2*E) is the wavelength

    Parameters
    ----------
    bfs        :  list of callables, continuum basis functions with energy E
    E          :  float, energy E in Hartree
    
    Optional
    --------
    rlarge     :  float, radius R in bohr where the integration starts
    """
    # number of basis functions
    Nbfs = len(bfs)
    # 'overlap' matrix
    S = np.zeros((Nbfs,Nbfs))

    # angular grid
    Lmax, (th,ph,angular_weights) = select_angular_grid(settings.lebedev_order)
    Nang = len(th)
    sc = np.sin(th)*np.cos(ph)
    ss = np.sin(th)*np.sin(ph)
    c  = np.cos(th)

    # wavelength
    k = np.sqrt(2*E)
    lamb = 2.0*np.pi/k
    print "k = %s" % k
    
    # radial grid on the interval [rlarge,rlarge+lambda]
    Nr = 1000
    r = np.linspace(rlarge,rlarge+lamb, Nr)
    # r^2 dr
    radial_weights = r**2 * np.ediff1d(r, to_end=r[-1]-r[-2])

    # cartesian coordinates of grid
    x = outerN(r, sc).flatten()
    y = outerN(r, ss).flatten()
    z = outerN(r, c ).flatten()
    weights = outerN(radial_weights, 4.0*np.pi * angular_weights).flatten()
    #
    Npts = Nr*Nang

    for i in range(0, Nbfs):
        phi_i = bfs[i]
        for j in range(0, Nbfs):
            phi_j = bfs[j]
            # If the continuum orbitals are normalized such that the
            # radial functions have the form
            #          r->oo
            #  phi(r) --------> 1/r sin(kr + delta_l)
            #
            # then the diagonal elements of the overlap matrix
            # should be
            #  S[i,i] = pi/k
            S[i,j] = np.sum( phi_i(x,y,z) * phi_j(x,y,z) * weights )

    print "overlap (unnormalized)"
    print S
    
    # For printing we normalize the overlap matrix
    # so that the diagonal elements are equal to 1.
    print "overlap (normalized)"
    print (S * k / np.pi)
            
    return S
            
def residual2_matrix(atomlist, potential, ps, bs):
    """
    compute matrix elements of (H-E)^2
                2
       <u |(H-E) |u > = <u |(V-V )(V-V )|u >
         i         j      i     i     j   j

    Parameters
    ----------
    atomlist  :  list of tuples (Zat,[x,y,z]) with atomic numbers
                 and positions
    potential :  callable, potential(x,y,z) evaluates the effective
                 molecular Kohn-Sham potential V
    ps        :  instance of `AtomicPotentialSet` with the effective
                 atomic Kohn-Sham potentials V_i
    bs        :  instance of `AtomicScatteringBasisSet` with the atomic
                 continuum orbitals |ui> (all for the same energy)
    """
    # number of basis functions
    nb = len(bs.bfs)
    R = np.zeros((nb,nb))
    
    for i in range(0, nb):
        # i-th basis function
        ui = bs.bfs[i]
        # effective potential of the atom to which the basis function belongs
        I = ui.atom_index
        poti = ps.pots[I]
        for j in range(i, nb):
            print "computing R[%d,%d] = " % (i,j), 
            uj = bs.bfs[j]
            J = uj.atom_index
            potj = ps.pots[J]


            ### DEBUG
            debug = 0
            if debug > 0:
                import matplotlib.pyplot as plt
                r = np.linspace(-5.0, 5.0, 100000)
                x = 0.0*r
                y = 0.0*r
                z = r
                plt.plot(r, potential(x,y,z), label=r"$V$")
                plt.plot(r, poti(x,y,z), label=r"$V_{%d}$" % (I+1))
                plt.plot(r, potj(x,y,z), ls="--", label=r"$V_{%d}$" % (J+1))            
                plt.plot(r, (potential(x,y,z) - poti(x,y,z))*(potential(x,y,z) - potj(x,y,z)), ls="-.", label=r"$(V-V_{%d})(V-V_{%d})$" % (I+1,J+1))
                VmVi = potential(x,y,z)-poti(x,y,z)
                xi,yi,zi = poti.center
                ri = np.sqrt((x-xi)**2+(y-yi)**2+(z-zi)**2)
                VmVi[ri < 2*poti.rmin] = 0.0
            
                VmVj = potential(x,y,z)-potj(x,y,z)
                xj,yj,zj = potj.center
                rj = np.sqrt((x-xj)**2+(y-yj)**2+(z-zj)**2)
                VmVi[rj < 2*potj.rmin] = 0.0

                plt.plot(r, VmVi*VmVj, ls="--", label=r"$(V-V_{%d})(V-V_{%d})$ (outer region only)" % (I+1,J+1))
                
                plt.ylim((-100.0, +100.0))
                plt.legend()
                plt.show()
            ###

            
            def integrand(x,y,z):
                V = potential(x,y,z)
                Vi = poti(x,y,z)
                Vj = potj(x,y,z)
                # Very close to each nucleus, the molecular Kohn-Sham potential
                # should be dominated by the nuclear attraction potential -Z/r.
                # In this region the atomic and molecular Kohn-Sham potentials should
                # cancel exactly, V-Vi = 0. However, subtracting two large numbers that
                # tend to infinity, will never give exactly 0 due to numerical errors.
                # Therefore a circle of radius 2*rmin is excised around atom i, where
                #  V-Vi is explicitly set to zero.
                VmVi = V-Vi
                xi,yi,zi = poti.center
                ri = np.sqrt((x-xi)**2+(y-yi)**2+(z-zi)**2)
                VmVi[ri < 2*poti.rmin] = 0.0

                VmVj = V-Vj
                xj,yj,zj = potj.center
                rj = np.sqrt((x-xj)**2+(y-yj)**2+(z-zj)**2)
                VmVi[rj < 2*potj.rmin] = 0.0

                # ui (V-Vi) (V-Vj) uj
                return ui(x,y,z)*VmVi*VmVj*uj(x,y,z)
#                return ui(x,y,z)*(V-Vi)*(V-Vj)*uj(x,y,z)

            R[i,j] = integral(atomlist, integrand)
            R[j,i] = R[i,j]

            print R[i,j]

    return R


def variational_mixture_continuum(atomlist, phi, delta_phi, potential_ee, E,
                                  rsmall=0.0, rlarge=10.0):
    """
    Since the orbital correction is only approximate due to the 
    spherical averaging, for the improved wavefunction we make
    the ansatz
              
          psi(r)  = a phi(r) + b delta_phi(r)

    where the coefficients a and b are optimized variationally.
    
    a,b are the coefficients belonging to the smallest in magnitude
    eigenvalue of the generalized eigenvalue equation

      R v = S v

    or 

      ( Raa Rab ) (a)     ( Saa Sab ) ( a )
      (         ) ( ) = E (         ) (   )
      ( Rba Rbb ) (b)     ( Sba Sbb ) ( b )

    where Raa = <phi|(H-E)^2|phi>, Rab = <phi|(H-E)^2|delta_phi>, etc.
    are the matrix elements of the residual operator (H-E)^2 and Sab
    is the overlap between the two wavefunctions integrated over one
    one wavelength at a large distance.

    Parameters
    ==========
    atomlist      :  list of tuples (Zat,[x,y,z]) with atomic numbers
                     and positions
    phi           :  callable, phi(x,y,z) evaluates the initial orbital
    delta_phi     :  callable, delta_phi(x,y,z) evaluates the orbital correction
    potential_ee  :  callable, potential_ee(x,y,z) evaluates the effective
                     molecular Kohn-Sham potential WITHOUT the nuclear attraction,
                     Vee = Vcoul + Vxc
    E             :  float, energy of continuum orbital, E=1/2 k^2

    Optional
    ========
    rsmall     :  regions closer than `rsmall` to any nucleus are excluded
                  from the integration
    rlarge     :  regions further away from the origin than `rlarge` are also
                  excluded, the overlap matrix is compute on the interval
                  [rlarge, rlarge+2*pi/k]

    Returns
    =======
    a,b        :  floats, mixing coefficients
    """
    # (H-E)phi
    residual_a = residual_ee_func(atomlist, phi, potential_ee, E)    
    # (H-E)delta_phi
    residual_b = residual_ee_func(atomlist, delta_phi, potential_ee, E)

    residuals = [residual_a, residual_b]
    n = len(residuals)
    
    R = np.zeros((n,n))
    
    for a in range(0, n):
        for b in range(0, n):
            def integrand(x,y,z):
                integ = residuals[a](x,y,z) * residuals[b](x,y,z)
                # Around the nuclei, it is very difficult to accuratly
                # calculate the residual (H-E)phi, because both the nuclear attraction (negative)
                # and the kinetic energy (positive) are very large, but cancel mostly,
                # so that we get very large numerical errors. To avoid this, we excise a circular
                # region of radius rsmall around each nucleus and set the integrand
                # to zero in this region, so as to exclude regions with large errors from
                # the calculation.
                for i,(Zi,posi) in enumerate(atomlist):
                    # distance of point (x,y,z) from nucleus i
                    xi,yi,zi = posi
                    ri = np.sqrt((x-xi)**2+(y-yi)**2+(z-zi)**2)
                    integ[ri < rsmall] = 0.0
                # Since density of the grid decreases very rapidly far away from the
                # molecule, we cannot expect the residual to be accurate at large
                # distances. Therefore we exclude all points outside the radius `rlarge`
                # from the integraion
                r = np.sqrt(x**2+y**2+z**2)
                integ[rlarge < r] = 0.0
                    
                return integ

            # compute the integrals R_ab = <a|(H-E)(H-E)|b> where |a> and |b>
            # are phi or delta_phi
            print "integration R[%d,%d]" % (a,b)
            R[a,b] = integral(atomlist, integrand)

    print "matrix elements of R = (H-E)^2"
    print R

    # 
    S = continuum_overlap([phi, delta_phi], E, rlarge=rlarge)
    print "continuum overlap matrix"
    print S
    
    # solve generalized eigenvalue problem
    eigvals, eigvecs = sla.eigh(R, S)
    print "eigenvalues"
    print eigvals
    print "eigenvectors"
    print eigvecs
    
    # select eigenvector with smallest (in magnitude) eigenvalue
    imin = np.argmin(abs(eigvals))

    a,b = eigvecs[:,imin]

    # Eigenvectors come out with arbitrary signs,
    # the phase of phi should always be positive.
    if a < 0.0:
        a *= -1.0
        b *= -1.0
    
    return a,b


def improve_continuum_orbital(atomlist, phi0, potential_ee, E,
                              max_iter=100, thresh=1.0e-6):
    """
    improves an initial guess phi0 for a continuum orbital with energy
    E by repeatedly computing orbital corrections 

         phi_{i+1} = a_i * phi_{i} + b_i * dphi_{i}

    The orbital correction delta_phi is obtained as the solution of
    the inhomogeneous Schroedinger equation

        (H-E) dphi  = -(H-E) phi
                  i             i

    The Schroedinger equation is solved approximately by replacing the
    effective potential with a spherical average around each atom. The
    best linear combination of phi and the correction is chosen. The 
    coefficients a_i and b_i minimize the expectation value of the residual
                  2
         R = (H-E)
                                     2
       a ,b  = argmin  <phi   |(H-E)  |phi   >
        i  i               i+1            i+1

    Parameters
    ==========
    atomlist      :  list of tuples (Zat,[x,y,z]) with atomic numbers
                     and positions
    phi0          :  callable, phi0(x,y,z) evaluates the initial orbital
    potential_ee  :  callable, potential_ee(x,y,z) evaluates the effective
                     molecular Kohn-Sham potential WITHOUT the nuclear attraction,
                     Vee = Vcoul + Vxc
    E             :  float, energy of continuum orbital

    Optional
    ========
    max_iter   :  int, maximum number of orbital correction
    thresh     :  float, the iteration is stopped as soon as the weight
                  of the orbital correction falls below this threshold:
                        2 
                     b_i  <= thresh
                       
    Returns
    =======
    phi        :  callable, phi(x,y,z) evaluates the improved orbital
    
    """
    # attraction potential between nuclei and electrons
    nuclear_potential = nuclear_potential_func(atomlist)
    # effective potential (electron-electron interaction + nuclei-electrons)
    def potential(x,y,z):
        return potential_ee(x,y,z) + nuclear_potential(x,y,z)
    
    print "orbital corrections..."

    phi = phi0
    for i in range(0, max_iter):
        residual = residual_ee_func(atomlist, phi, potential_ee, E)
        def source(x,y,z):
            return -residual(x,y,z)
        # orbital correction
        delta_phi = inhomogeneous_schroedinger(atomlist, potential, source, E)

        a,b = variational_mixture_continuum(atomlist, phi, delta_phi, potential_ee, E)

        # The error is estimated as the norm^2 of the orbital correction:
        #   error = <delta_phi|delta_phi>
        # For large radii, we cannot expect the solution to be correct since
        # the density of points is far too low. Therefore we exclude all points
        # outside the radius `rlarge` from the integration.
        rlarge = 10.0
        def error_density(x,y,z):
            err = abs(delta_phi(x,y,z))**2
            r = np.sqrt(x**2+y**2+z**2)
            err[rlarge < r] = 0.0
            return err
        error = integral(atomlist, error_density)
        print "  iteration i= %d   |dphi|^2= %e   b^2= %e  (threshold= %e )" % (i+1, error, b**2, thresh)
        # The solution is converged, when the weight of the
        # correction term b**2 is small enough.
        if b**2 < thresh:
            print "CONVERGED"
            break
        
        # next approximation for phi
        phi_next = add_two_functions(atomlist, phi, delta_phi, a, b)

        ### DEBUG
        import matplotlib.pyplot as plt
        plt.clf()
        # plot cuts along z-axis
        r = np.linspace(-40.0, 40.0, 10000)
        x = 0.0*r
        y = 0.0*r
        z = r

        plt.title("Iteration i=%d" % (i+1))
        plt.xlabel("z / bohr")
        plt.ylim((-0.6, +0.6))
        
        plt.plot(r, phi(x,y,z), label=r"$\phi$")
        plt.plot(r, delta_phi(x,y,z), ls="--", label=r"$\Delta \phi$")
        plt.plot(r, phi_next(x,y,z), label=r"$a \phi + b \Delta \phi$")
        plt.plot(r, residual(x,y,z), ls="-.", label="residual $(H-E)\phi$")
        plt.legend()
        plt.savefig("/tmp/iteration_%3.3d.png" % (i+1))
#        plt.show()

        ###
        
        # prepare for next iteration
        phi = phi_next
    else:
        msg = "Orbital corrections did not converge in '%s' iterations!" % max_iter
        print "WARNING: %s" % msg
        #raise RuntimeError(msg)    
    return phi

def relax_continuum_orbital(atomlist, phi0, potential_ee, E,
                            max_iter=100, thresh=1.0e-6):
    # DOES NOT WORK
    """
    improves an initial guess phi0 for a continuum orbital with energy
    E by iteratively solving the following linear equation

              2             2
         (H-E) dphi = -(H-E)  phi
                                 0
    for the orbital correction dphi. The improved orbital is phi0+dphi.
    The linear equation has the form

         A.x = b

    and can be solve using by modified Richardson iteration

          (k+1)     (k)              (k)
         x      =  x    +  w (b - A.x   )

    for some parameter w > 0.


    Parameters
    ==========
    atomlist      :  list of tuples (Zat,[x,y,z]) with atomic numbers
                     and positions
    phi0          :  callable, phi0(x,y,z) evaluates the initial orbital
    potential_ee  :  callable, potential_ee(x,y,z) evaluates the effective
                     molecular Kohn-Sham potential WITHOUT the nuclear attraction,
                     Vee = Vcoul + Vxc
    E             :  float, energy of continuum orbital

    Optional
    ========
    max_iter   :  int, maximum number of iterations
    thresh     :  float, the iteration is stopped as soon as the error
                  falls below this threshold: |b-A.x| < thresh
                       
    Returns
    =======
    phi        :  callable, phi(x,y,z) evaluates the improved orbital
    
    """
    w = 1.0
    # right hand side of A.x = b
    #  b = (H-E)^2 phi0
    residual1 = residual_ee_func(atomlist, phi0, potential_ee, E)
    residual2 = residual_ee_func(atomlist, residual1, potential_ee, E)
    b = residual2

    def null(x,y,z):
        return 0*x
    
    x = null
    for i in range(0, max_iter):
        # compute A.x^(k) = (H-E)^2 phi_k
        print "residuals..."
        residual1 = residual_ee_func(atomlist, x, potential_ee, E)
        residual2 = residual_ee_func(atomlist, residual1, potential_ee, E)
        Ax = residual2
        # b - A.x^(k)
        print "error..."
        bmAx = add_two_functions(atomlist, b, Ax, 1.0, -1.0)
        error = np.sqrt(overlap(atomlist, bmAx, bmAx))
        print "iteration i= %d  error |b-Ax|= %e" % (i, error)
        if error < thresh:
            print "CONVERGED"
            break
        
        # x^(k+1) = x^(k) + w * (b - A.x^(k))
        x_next = add_two_functions(atomlist, x, bmAx, 1.0, w)

        x = x_next

        phi = add_two_functions(atomlist, phi0, x, 1.0, 1.0)
        residual = residual_ee_func(atomlist, phi, potential_ee, E)

        import matplotlib.pyplot as plt
        plt.clf()
        plt.title("Iteration i=%d" % (i+1))
        plt.xlabel("z / bohr")
        plt.ylim((-0.6, +0.6))

        r = np.linspace(-10.0, 10.0, 1000)
        plt.plot(r, phi(0*r,0*r,r), label=r"$\phi$")
        plt.plot(r, residual(0*r,0*r,r), ls="-.", label="residual $(H-E)\phi$")

        plt.legend()
        plt.savefig("/tmp/relaxation_%3.3d.png" % (i+1))
        plt.show()
        
    else:
        msg = "Orbital relaxation did not converge in '%s' iterations!" % max_iter
        print "WARNING: %s" % msg
        #raise RuntimeError(msg)    
    return phi


def averaged_angular_distribution(atomlist, bound_orbs, continuum_orbs, E):
    nb = len(bound_orbs)
    nc = len(continuum_orbs)

    # compute dipole matrix elements between bound and free orbitals
    dipoles_bf = np.zeros((nb,nc,3))

    for i in range(0, nb):
        for j in range(0, nc):
            dipoles_bf[i,j,:] = electronic_dipole(atomlist, bound_orbs[i], continuum_orbs[j])
            print "dipole <bound %d|er|free %d> = %s" % (i+1,j+1, dipoles_bf[i,j,:])

    # wrappers
    class AtomicBasisFunction:
        def __init__(self, func):
            self.func = func
        def amp(self, x,y,z):
            return self.func(x,y,z)
        
    class ContinuumBasisSet:
        def __init__(self, bfs, E):
            self.bfs = bfs
            self.E = E

    bfs = [AtomicBasisFunction(orb) for orb in continuum_orbs]
    bs_free = ContinuumBasisSet(bfs, E)
    
    from DFTB.Scattering import PAD
    Rmax = 5000.0
    orientation_averaging = PAD.OrientationAveraging_small_memory(dipoles_bf, bs_free, Rmax, E, npts_euler=5, npts_theta=1000)

    # coefficients of HOMO in the MO basis
    mo_bound = np.zeros(nb)
    mo_bound[-1] = 1.0
    pad, betas = orientation_averaging.averaged_pad(mo_bound)

    print "betas"
    print betas


    
#############################################################
#
# Testing
#
#############################################################

def test_AO_basis(atomlist, bs, ps, E):
    """
    check that atomic continuum orbitals are solutions of the atomic
    Kohn-Sham equations

         (T + Vi - E) \psi_i = 0
    """
    import matplotlib.pyplot as plt
    r = np.linspace(-5.0, 5.0, 10000)
    x = 0.0*r
    y = 0.0*r
    z = r
    
    nb = len(bs.bfs)
    
    for i in range(0, nb):
        ui = bs.bfs[i]
        # effective potential of the atom to which the basis function belongs
        I = ui.atom_index
        poti = ps.pots[I]

        # compute residual (T+Vi-E)ui
        residual_i = residual_func([atomlist[I]], ui, poti, E)

        l, = plt.plot(r, ui(x,y,z), label=r"$\psi_{%d}$" % (i+1))
        plt.plot(r, residual_i(x,y,z), ls="-.", label=r"$(T+V_{%d}-E)\phi_{%d}$" % (i+1,i+1), color=l.get_color())

    plt.legend()
    plt.show()
    

def test_h2_continuum_orbital():
    """
    The sigma continuum orbital is approximated as a linear
    combination of two s continuum orbitals and is then
    improved iteratively
    """
    # bond length in bohr
    dist = 2.0
    # positions of protons
    posH1 = (0.0, 0.0, -dist/2.0)
    posH2 = (0.0, 0.0, +dist/2.0)
    
    atomlist = [(1, posH1),
                (1, posH2)]    

    # Set resolution of multicenter grid
    settings.radial_grid_factor = 20
    settings.lebedev_order = 23
    
    # energy of continuum orbital
    E = 1.0

    # same functional as used in the calculation of pseudo orbitals
    xc = XCFunctionals.libXCFunctional(Parameters.pseudo_orbital_x, Parameters.pseudo_orbital_c)
    dft = BasissetFreeDFT(atomlist, xc)

    print "initial orbital guess from DFTB calculation"
    orbitals = dft.getOrbitalGuess()
    
    norb = len(orbitals)
    # all orbitals are doubly occupied
    nelec = 2*norb

    bound_orbitals = dft.getOrbitalGuess()
    
    # electron density (closed shell)
    rho = density_func(bound_orbitals)
    # effective Kohn-Sham potential
    veff = effective_potential_func(atomlist, rho, xc,
                                    nelec=nelec,
                                    nuclear=True)
    # effective Kohn-Sham potential without nuclear attraction
    # (only electron-electron interaction)        
    veff_ee = effective_potential_func(atomlist, rho, xc,
                                       nelec=nelec,
                                       nuclear=False)
    
    ps = AtomicPotentialSet(atomlist)

    lmax = 0
    bs = AtomicScatteringBasisSet(atomlist, E, lmax=lmax)

    #test_AO_basis(atomlist, bs, ps, E)
    
    R = residual2_matrix(atomlist, veff, ps, bs)
    S = continuum_overlap(bs.bfs, E)
    print "continuum overlap"
    print S
    print "residual^2 matrix"
    print R        
        
    eigvals, eigvecs = sla.eigh(R, S)
    print eigvals
    print "eigenvector belonging to lowest eigenvalue"
    print eigvecs[:,0]

    # LCAO continuum orbitals
    continuum_orbitals = orbital_transformation(atomlist, bs.bfs, eigvecs)

    # improve continuum orbital by adding a correction term
    #
    #    phi = phi0 + dphi
    #
    # The orbital correction dphi is the solution of the inhomogeneous
    # Schroedinger equation
    #
    #   (H-E)dphi = -(H-E)phi0
    #
    print "orbital correction..."
    phi0 = continuum_orbitals[0]
    
    phi = improve_continuum_orbital(atomlist, phi0, veff_ee, E)


def test_hmi_continuum():
    """
    check that the continuum wavefunction of H2+ really are solutions
    of Schroedinger's equation, i.e. have (H-E)\phi = 0 everywhere

    starting from an LCAO guess for the continuum orbital, try to find
    the exact solution by adding orbital corrections iteratively
    """
    # First we compute the exact wavefunction of the hydrogen molecular ion.
    from DFTB.Scattering import HMI

    # The bond length and charges cannot be changed, since the
    # separation constants were calculated only for the H2+ ion at R=2!
    R = 2.0
    Za = 1.0
    Zb = 1.0

    # energy of continuum orbital
    E = 0.5

    ## sigma (m=0) orbital
    m = 0
    n = 0
    trig = 'cos'

    # separation constant
    Lsep = HMI.SeparationConstants(R, Za, Zb)
    Lsep.load_separation_constants()
    Lfunc = Lsep.L_interpolated(m,n)
    
    c2 = 0.5*E*R**2
    mL,nL,L = Lfunc(c2)
    
    parity = (-1)**(mL+nL)
    phi_exact = HMI.create_wavefunction(mL,L,R*(Za+Zb),0.0,R, c2,parity, trig)

    # Old implementation of H2+ wavefunctions, the wavefunction
    # looks indistinguishable from the exact wavefunction, but the
    # non-zero residue shows that is contains large errors.
    from DFTB.Scattering.hydrogen_molecular_ion import DimerWavefunctions
    wfn = DimerWavefunctions(R,Za,Zb, plot=False)

    delta, (Rfunc,Sfunc,Pfunc),wavefunction_exact = wfn.getContinuumOrbital(m,n,trig,E)
    def phi_exact_DW(x,y,z):
        return wavefunction_exact((x,y,z), None)

    
    # Set resolution of multicenter grid
    settings.radial_grid_factor = 10
    settings.lebedev_order = 41
    
    # Next we compute the wavefunction using the basis set free method
    atomlist = [(int(Za), (0.0, 0.0, -R/2.0)),
                (int(Zb), (0.0, 0.0, +R/2.0))]

    # no other electrons, only nuclear potential
    def potential(x,y,z):
        nuc = 0.0*x
        for Zi,posi in atomlist:
            ri = np.sqrt((x-posi[0])**2 + (y-posi[1])**2 + (z-posi[2])**2)
            nuc += -Zi/ri
        return nuc

    # electron-electron interaction 
    def potential_ee(x,y,z):
        return 0.0*x

    # Set resolution of multicenter grid
    settings.radial_grid_factor = 10
    settings.lebedev_order = 41
    
    # residual of exact wavefunction (should be zero)
    residual_exact = residual_func(atomlist, phi_exact, potential, E)
    residual_ee_exact = residual_ee_func(atomlist, phi_exact, potential_ee, E)
    residual_exact_DW = residual_func(atomlist, phi_exact_DW, potential, E)

    # Laplacian
    laplacian_exact = laplacian_func(atomlist, phi_exact)
    
    import matplotlib.pyplot as plt
    plt.clf()
    r = np.linspace(-15.0, 15.0, 5000)
    x = 0*r
    y = 0*r
    z = r

    # plot exact wavefunction
    plt.plot(r, phi_exact(x,y,z), label="$\phi$ exact")

    # 
    phi_exact_xyz = phi_exact(x,y,z)
    phi_exact_DW_xyz = phi_exact_DW(x,y,z)
    scale = phi_exact_xyz.max() / phi_exact_DW_xyz.max()
    plt.plot(r, scale * phi_exact_DW_xyz, label="$\phi$ exact (DimerWavefunction)")

    # and residual
    plt.plot(r, residual_exact(x,y,z), label=r"$(H-E)\phi$ (exact, old)")
    plt.plot(r, residual_exact_DW(x,y,z), ls="-.", label=r"$(H-E)\phi$ (exact, DimerWavefunction, old)")

    plt.plot(r, residual_ee_exact(x,y,z), ls="--", label=r"$(H-E)\phi$ (exact, new)")
    # kinetic energy 
    plt.plot(r, -0.5*laplacian_exact(x,y,z), ls="--", label=r"$-\frac{1}{2}\nabla^2 \phi$")
    # potential energy
    plt.plot(r, (potential(x,y,z)-E)*phi_exact(x,y,z), ls="--", label=r"$(V-E)\phi$")
    

    ## The initial guess for the \sigma continuum orbital
    ## is a regular Coulomb function centered on the midpoint
    ## between the two protons.
    #phi0 = regular_coulomb_func(E, +2, 0, 0, 0.0, center=(0.0, 0.0, 0.0))

    """
    ## The initial guess for the \sigma continuum orbital is
    ## the sum of two hydrogen s continuum orbitals
    bs = AtomicScatteringBasisSet(atomlist, E, lmax=0)
    
    phi0 = add_two_functions(atomlist,
                             bs.bfs[0], bs.bfs[1],
                             1.0/np.sqrt(2.0), 1.0/np.sqrt(2.0))
    """
    """
    ## start with exact wavefunction
    phi0 = phi_exact
    """
    
    # The initial guess for the \sigma continuum orbital is
    # a hydrogen continuum orbital in the center
    bs = AtomicScatteringBasisSet([(1, (0.0, 0.0, 0.0))], E, lmax=0)
    phi0 = bs.bfs[0]
    
    plt.plot(r, phi0(x,y,z), ls="-.", label="guess $\phi_0$")
    
    plt.legend()
    plt.show()
       

    #phi = improve_continuum_orbital(atomlist, phi0, potential_ee, E, thresh=1.0e-6)
    phi = relax_continuum_orbital(atomlist, phi0, potential_ee, E, thresh=1.0e-6)

    import matplotlib.pyplot as plt
    plt.clf()
    r = np.linspace(-15.0, 15.0, 5000)
    x = 0*r
    y = 0*r
    z = r

    phi_exact_xyz = phi_exact(x,y,z)
    phi_xyz = phi(x,y,z)
    # scale numerical phi such that the maxima agree
    scale = phi_exact_xyz.max() / phi_xyz.max()
    phi_xyz *= scale
    print "scaling factor  s = %s" % scale
    
    plt.plot(r, phi_exact_xyz, label="exact")
    plt.plot(r, phi_xyz, label="numerical")
    plt.legend()
    
    plt.show()
    
    
    
if __name__ == "__main__":
    test_hmi_continuum()
    #test_h2_continuum_orbital()

    
