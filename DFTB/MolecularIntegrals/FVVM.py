#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Finite-Volume Variational Method

References
----------
[1] H. Le Rouzo, G. Raseev, "Finite-volume variational method: First application to direct molecular photoionization",
    Phys. Rev. A, Vol. 29, Num. 3, pp. 1214-1223 (1984)
"""

import numpy as np
import numpy.linalg as la
import scipy.linalg as sla

from DFTB import AtomicData, Parameters
from DFTB.MolecularIntegrals.Ints1e import integral, overlap
from DFTB.MolecularIntegrals.MulticenterIntegration import atomlist2arrays, multicenter_gradient_product, select_angular_grid
from DFTB.MolecularIntegrals import settings

from DFTB.BasisSets import AtomicBasisSet
from DFTB.Scattering.SlakoScattering import AtomicScatteringBasisSet
from DFTB.SlaterKoster import XCFunctionals
from DFTB.MolecularIntegrals.BasissetFreeDFT import BasissetFreeDFT, density_func, effective_potential_func, orbital_transformation, energy


def vdw_sphere_radius(atomlist, fac=2.0):
    """
    find the radius r0 such that the van der Waals volume of the 
    molecule is entirely contained inside a sphere of radius r0
    around the origin

    Parameters
    ----------
    atomlist   :  list of tuples (Zat,[x,y,z]) with atomic numbers
                  and positions

    Optional
    --------
    fac        :  float, factor by which all vdW radii are increased

    Returns
    -------
    r0         :  float, radius of sphere
    """
    # The vdW volume is the intersection of spheres with atom-specific
    # vdW radii around each atom. Find the smallest radius r0
    # such that all atomic spheres are contained in the big molecular sphere.
    r0 = 0.0
    for Zi,posi in atomlist:
        ri_vdw = AtomicData.vdw_radii[AtomicData.atom_names[Zi-1]]
        ri = la.norm(posi)
        r0 = max(r0, ri + ri_vdw)

    return r0


def gradient_product(atomlist, f1, f2):
    """
    create function for evaluating the product of the gradients
    of two functions
           __      __
          (\/ f1).(\/ f2)

    Parameters
    ==========
    atomlist           : list of tuples (Zat,[xi,yi,zi]) with atomic numbers and nuclear coordinates
    f1,f2              : callable, f1(x,y,z) and f2(x,y,z) should evaluate the functions at the 
                         grid points specified by x = [x0,x1,...,xn], y = [y0,y1,...yn]
                         and z = [z0,z1,...,zn]

    Returns
    =======
    grad_prod_func     : callable, evaluates (grad f1).(grad f2)(x,y,z)
    """
    # Translate molecular geometry into the format understood
    # by the multicenter integration code
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    grad_prod_func = multicenter_gradient_product(f1, f2,
                                                  atomic_coordinates, atomic_numbers,
                                                  cusps_separate=True,
                                                  radial_grid_factor=settings.radial_grid_factor,
                                                  lebedev_order=settings.lebedev_order)

    return grad_prod_func


def fvvm_matrix_elements(atomlist, bfs, potential, E, r0):
    """
    compute matrix elements of eqns. (8) and (9) in Ref. [1] which are
    needed for the finite-volume variational method. The integration
    is limited to a sphere of radius r0 around the origin.

    Parameters
    ==========
    atomlist     :  list of tuples (Z,[x,y,z]) with atomic numbers
                    and positions
    bfs          :  list of callables, `n` real basis functions 
    potential    :  callable, potential(x,y,z) evaluates effective potential
    E            :  float, energy of continuum orbital
    r0           :  float, radius of sphere which defines the integration volume

    Returns
    =======
    A            :  n x n matrix
    D            :  n x n matrix
    S            :  n x n matrix with overlap in finite volume
    """
    # number of basis functions
    nbfs = len(bfs)
    print " %d basis functions" % nbfs

    # volume integral
    #       __     __
    # A  = <\/chi |\/chi >   +  2 <chi |(V-E)|chi >
    #  ij        i      j V           i          j V
    A = np.zeros((nbfs,nbfs))
    #
    # surface integral
    #  D  = <chi |chi >
    #   ij      i    j surface of V
    D = np.zeros((nbfs,nbfs))
    #
    # overlap inside volume
    #
    # S   = <chi |chi >
    #  ij       i    j V
    S = np.zeros((nbfs,nbfs))
    
    # angular grid for surface integral
    Lmax, (th,ph,angular_weights) = select_angular_grid(settings.lebedev_order)
    # cartesian coordinates of Lebedev points on sphere of
    # radius r0
    xS = r0*np.sin(th)*np.cos(ph)
    yS = r0*np.sin(th)*np.sin(ph)
    zS = r0*np.cos(th)
    
    for i in range(0, nbfs):
        # basis function i
        chi_i = bfs[i]
        for j in range(i, nbfs):
            # basis function j
            print " integrals between basis functions i= %d and j= %d" % (i+1,j+1)
            chi_j = bfs[j]

            # volume integral A_ij
            grad_prod_func = gradient_product(atomlist, chi_i, chi_j)
            def integrand(x,y,z):
                # (grad chi_i).(grad chi_j)
                integ = grad_prod_func(x,y,z)
                # 2 * chi_i * chi_j (V-E)
                integ += 2 * chi_i(x,y,z) * chi_j(x,y,z) * (potential(x,y,z) - E)
                # Outside the integration volume the integrand
                # is set to zero.
                r = np.sqrt(x*x+y*y+z*z)
                integ[r0 < r] = 0.0
                return integ

            A[i,j] = integral(atomlist, integrand)
            # A is symmetric
            A[j,i] = A[i,j]

            # surface integral D_ij
            D[i,j] = 4.0*np.pi * r0**2 * np.sum(angular_weights*chi_i(xS,yS,zS)*chi_j(xS,yS,zS))
            # D is symmetric
            D[j,i] = D[i,j]

            # overlap integral S_ij inside the volume
            def integrand(x,y,z):
                integ = chi_i(x,y,z)*chi_j(x,y,z)
                # Outside the integration volume the integrand
                # is set to zero.
                r = np.sqrt(x*x+y*y+z*z)
                integ[r0 < r] = 0.0
                return integ
            
            S[i,j] = integral(atomlist, integrand)
            # S is symmetric
            S[j,i] = S[i,j]
            
    return A, D, S

##################################################
#
# Testing
#
##################################################

def test_hydrogen_atom():
    """
    reproduce some results of Ref. [1] for the hydrogen atom 
    """
    # hydrogen atom
    atomlist = [(1, (0.0, 0.0, 0.0))]

    # Coulomb potential
    def potential(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        V = -1.0/r
        return V
    
    # Set resolution of multicenter grid
    settings.radial_grid_factor = 40
    settings.lebedev_order = 41

    # energy of continuum orbital
    E = 0.5

    # basis set consists of radial functions sin(k_i*r)
    # and Y_10(th,ph) angular terms.
    # The k_i's are equally spaced around the asymptotic value
    # k = sqrt(2*E)= 1 (see Ref.[1] table I)
    ks = np.linspace(0.3, 1.3, 6)

    def create_basis_function(ki):
        def basis_func(x,y,z):
            r = np.sqrt(x*x+y*y+z*z)
            Y10 = np.sqrt(3.0/(4.0*np.pi)) * z/r
            return 1.0/r * np.sin(ki*r)*Y10
        return basis_func

    # list of basis functions
    bfs = []
    for ki in ks:
        bfs.append( create_basis_function(ki) )

    # radius that separates inner from outer region
    r0 = 10.0

    A, D, S = fvvm_matrix_elements(atomlist, bfs, potential, E, r0)

    print "A matrix"
    print A
    print "D matrix"
    print D
    print "S matrix"
    print S
    
    # solve generalized eigenvalue problem
    #   A.C = b.D.C
    # Since D does not have full rank, M <= N, we need to transform
    # the eigenvalue problem (see eqn. A8 in Ref. [1])
    #    -1          -1 
    #  (A   D).C = -b   C
    #
    mbinv, C = sla.eigh(np.dot(la.inv(A),D))
    b = -1.0/mbinv
    
    print "b's"
    print b

    print "C's"
    print C
    
    # continuum orbitals in inner region
    psi_in = orbital_transformation(atomlist, bfs, C)

    # energy expectation values
    for i,bi in enumerate(b):
        en = energy(atomlist, psi_in[i], psi_in[i], potential)
        print " bi= %e  energy= %e" % (bi, en)
    
    import matplotlib.pyplot as plt
    r = np.linspace(0.0, 10.0, 1000)
    x = 0.0*r
    y = 0.0*r
    z = r

    for i,bi in enumerate(b):
        plt.plot(r, psi_in[i](x,y,z), label="b= %2.2e" % (bi))

    plt.legend()
    plt.show()
    
    
def test():    
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
    
    # effective potential
    rho = density_func(bound_orbitals)
    veff = effective_potential_func(atomlist, rho, xc, nelec=nelec)
    
    # radius that separates inner from outer region
    r0 = vdw_sphere_radius(atomlist, fac=2.0)
    # basis sets
    bs_core    = AtomicBasisSet(atomlist, orbital_set="core")
    bs_valence = AtomicBasisSet(atomlist, orbital_set="valence")
    bs_continuum = AtomicScatteringBasisSet(atomlist, E, lmax=1)
    # combine basis functions from all basis sets
    bfs = bs_core.bfs + bs_valence.bfs + bs_continuum.bfs
    
    A, D, S = fvvm_matrix_elements(atomlist, bfs, veff, E, r0)
    
    # solve generalized eigenvalue problem
    #   A.C = b.D.C
    b, C = sla.eigh(A, D)
    
    return b, C


if __name__ == "__main__":
    test_hydrogen_atom()
    #test()
    
            
