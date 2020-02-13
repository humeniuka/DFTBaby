#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
compute matrix elements of single particle operators such as
 - the overlap (a|b)
 - the kinetic energy (a|T|b)
 - the nuclear attraction energy (a|sum_I -Z(I)/|r-R(I)| |b)
 - the dipole operator (a|e*r|b)
using Becke's multicenter integration scheme
"""
import numpy as np
import numpy.linalg as la

from DFTB.MolecularIntegrals.MulticenterIntegration import multicenter_integration, multicenter_laplacian, atomlist2arrays
# default parameters controlling the resolution of multicenter grid
from DFTB.MolecularIntegrals import settings

def integral(atomlist, f):
    """
    compute the integral
       /
       | f(r) dV
       / 

    Parameters
    ==========
    atomlist       :  list of tuples (Z,[x,y,z]) with atomic numbers
                      and positions
    f              :  callable, f(x,y,z)

    Returns
    =======
    integ          :  float, volume integral
    """
    # Bring geometry data into a form understood by the module MolecularIntegrals
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    # compute integral on a multicenter grid
    integ = multicenter_integration(f, atomic_coordinates, atomic_numbers,
                                  radial_grid_factor=settings.radial_grid_factor,
                                  lebedev_order=settings.lebedev_order)

    return integ
    
def overlap(atomlist, bfA, bfB):
    """
    overlap between two basis functions
    
        (a|b)

    Parameters
    ==========
    atomlist       :  list of tuples (Z,[x,y,z]) with atomic numbers
                      and positions
    bfA, bfB       :  instances of AtomicBasisFunction or callables, 
                      e.g. bfA(x,y,z) etc.

    Returns
    =======
    Sab            :  float, overlap integral    
    """
    # Bring geometry data into a form understood by the module MolecularIntegrals
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # Now we compute the integrals numerically on a multicenter grid.
    # 1. define integrand s = a b
    def s_integrand(x,y,z):
        return bfA(x,y,z).conjugate() * bfB(x,y,z)
    #                                    
    # 2. integrate density on a multicenter grid
    Sab = multicenter_integration(s_integrand, atomic_coordinates, atomic_numbers,
                                  radial_grid_factor=settings.radial_grid_factor,
                                  lebedev_order=settings.lebedev_order)

    return Sab
    
    
def kinetic(atomlist, bfA, bfB):
    """
    matrix element of kinetic energy
                   __2
         (a| - 1/2 \/  |b)

    Parameters
    ==========
    atomlist       :  list of tuples (Z,[x,y,z]) with atomic numbers
                      and positions
    bfA, bfB       :  instances of AtomicBasisFunction or callables, 
                      e.g. bfA(x,y,z) etc.

    Returns
    =======
    Tab            :  float, kinetic energy integral
    """
    # Bring data into a form understood by the module MolecularIntegrals
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # Now we compute the integrals numerically on a multicenter grid.
    #                   __2
    # 1. find Laplacian \/ b
    lapB = multicenter_laplacian(bfB, atomic_coordinates, atomic_numbers,
                                 radial_grid_factor=settings.radial_grid_factor,
                                 lebedev_order=settings.lebedev_order)
    #                                __2
    # 2. define integrand t = -1/2 a \/ b
    def t_integrand(x,y,z):
        return -0.5 * bfA(x,y,z).conjugate() * lapB(x,y,z)
    #                                    
    # 3. integrate kinetic energy density on multicenter grid
    Tab = multicenter_integration(t_integrand, atomic_coordinates, atomic_numbers,
                                  radial_grid_factor=settings.radial_grid_factor,
                                  lebedev_order=settings.lebedev_order)

    return Tab

def nuclear(atomlist, bfA, bfB):
    """
    matrix element of nuclear attraction
                     (-Zk)
         (a| sum_k -------- |b)
                    |r-Rk|
    Parameters
    ==========
    atomlist       :  list of tuples (Z,[x,y,z]) with atomic numbers
                      and positions
    bfA, bfB       :  instances of AtomicBasisFunction or callables, 
                      e.g. bfA(x,y,z) etc.

    Returns
    =======
    Nab            :  float, integral of nuclear attraction
    """
    # Bring data into a form understood by the module MolecularIntegrals
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # Now we compute the integrals numerically on a multicenter grid.
    #                          
    # 1. define integrand n = sum_k (-Zk)/|r-Rk| * a(r) * b(r)
    def nuc_integrand(x,y,z):
        # nuclear attraction potential
        nuc = 0.0*x
        for Zk,Rk in atomlist:
            # position of k-th nucleus
            X,Y,Z = Rk
            # Zk is the atomic number of k-th nucleus
            nuc -= Zk / np.sqrt( (x-X)**2 + (y-Y)**2 + (z-Z)**2 )
        # product of bra and ket wavefunctions
        rhoAB = bfA(x,y,z).conjugate() * bfB(x,y,z)
        nuc *= rhoAB
        
        return nuc
    #                                    
    # 2. integrate nuclear attraction energy density on multicenter grid
    Nab = multicenter_integration(nuc_integrand, atomic_coordinates, atomic_numbers,
                                  radial_grid_factor=settings.radial_grid_factor,
                                  lebedev_order=settings.lebedev_order)

    return Nab

def nuclear_repulsion(atomlist):
    """
    repulsion energy between nuclei
                          Z(i) Z(j)
       V_rep = sum  sum  -------------
                i   j>i   |R(i)-R(j)|

    Parameters
    ==========
    atomlist        :  list of tuples (Z,[x,y,z]) for each atom,
                       molecular geometry

    Returns
    =======
    Vrep            :  float, repulsion energy
    """
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    Nat = len(atomic_numbers)
    
    Vrep = 0.0
    for i in range(0, Nat):
        Zi = atomic_numbers[i]
        Ri = atomic_coordinates[:,i]
        for j in range(i+1,Nat):
            Zj = atomic_numbers[j]
            Rj = atomic_coordinates[:,j]
            
            Vrep += Zi*Zj / la.norm(Ri-Rj)
            
    return Vrep

def electronic_dipole(atomlist, bfA, bfB):
    """
    electric dipole between two basis functions
    
        (a|e*r|b)

    Parameters
    ==========
    atomlist       :  list of tuples (Z,[x,y,z]) with atomic numbers
                      and positions
    bfA, bfB       :  instances of AtomicBasisFunction or callables, 
                      e.g. bfA(x,y,z) etc.

    Returns
    =======
    Dab            :  numpy array with components [Dx,Dy,Dz] of dipole
                      matrix elements (a|e*x|b), (a|e*y|b), (a|e*z|b)
    """
    # Bring geometry data into a form understood by the module MolecularIntegrals
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    # 1. define integrand, xyz=0,1,2 selects the x,y or z-component
    def dipole_density(x,y,z, xyz=0):
        er = [x,y,z]
        return bfA(x,y,z) * er[xyz] * bfB(x,y,z)
    #                                    
    # 2. integrate density on a multicenter grid
    Dab = np.zeros(3, dtype=complex)
    for xyz in [0,1,2]:
        Dab[xyz] = multicenter_integration(lambda x,y,z: dipole_density(x,y,z, xyz=xyz),
                                           atomic_coordinates, atomic_numbers,
                                           radial_grid_factor=settings.radial_grid_factor,
                                           lebedev_order=settings.lebedev_order)
    return Dab.real
    

###############################################################################################
#
# TESTING
#
################################################################################################

def test_integrals_h2():
    # Gaussian-type 1s function with exponent alpha
    # centered at position Rc=(xc,yc,zc)
    def gf_1s(alpha, xc,yc,zc, x,y,z):
        return (2*alpha/np.pi)**(0.75) * np.exp(-alpha*( (x-xc)**2 + (y-yc)**2 + (z-zc)**2 ))
    
    # contracted STO-3G basis function for hydrogen, see eqn. (3.225) in Szabo & Ostlund
    def cgf_1s(xc,yc,zc, x,y,z):
        wfn =   0.444635 * gf_1s(0.168856, xc,yc,zc, x,y,z) \
              + 0.535328 * gf_1s(0.623913, xc,yc,zc, x,y,z) \
              + 0.154329 * gf_1s(3.42525 , xc,yc,zc, x,y,z)
        return wfn
        
    # molecular integrals for STO-3G hydrogen molecule
    # with an internuclear distance of R=1.4 bohr
    R = 1.4
    atomlist = [(1, (0.0, 0.0, 0.0)),
                (1, (0.0, 0.0, R))]

    # The basis set consists of two functions phi_1 and phi_2
    # centered on each proton
    def phi_1(x,y,z):
        return cgf_1s(0.0, 0.0, 0.0, x,y,z)
    def phi_2(x,y,z):
        return cgf_1s(0.0, 0.0, R  , x,y,z)
    basis = [phi_1, phi_2]
    nbfs = len(basis)

    # overlap, kinetic and nuclear matrix elements
    S = np.zeros((nbfs,nbfs))
    T = np.zeros((nbfs,nbfs))
    N = np.zeros((nbfs,nbfs))
    
    for i in range(0, nbfs):
        for j in range(0, nbfs):
            S[i,j] = overlap(atomlist, basis[i], basis[j])
            T[i,j] = kinetic(atomlist, basis[i], basis[j])
            N[i,j] = nuclear(atomlist, basis[i], basis[j])

    print ""
    print "  The matrix elements for STO-3G H2 should be compared"
    print "  with those in chapter 3.5.2 of Szabo & Ostlund."
    print ""
    print "overlap matrix S, compare with eqn. (3.229)"
    print S
    print "kinetic energy T, compare with eqn. (3.230)"
    print T
    print "nuclear attraction V, compare with sums of eqns. (3.231) and (3.232)"
    print N

def test_dipole_integrals():
    """
    compare numerical integrals for dipole matrix elements with
    exact values for Gaussian orbitals
    """
    from DFTB.MolecularIntegrals import gtobasis, fchkfile, integrals

    # load test data with Gaussian AO basis
    filename = "test/h2o_hf_sv.fchk"
    res = fchkfile.G09ResultsDFT(filename)

    # geometry
    atomlist = []
    for i in range(0, res.nat):
        Z = res.atomic_numbers[i]
        posi = res.coordinates[i,:]
        atomlist.append( (Z, posi) )
    
    # Dipole matrix elements are computed from analytic expressions
    # for Gaussian orbitals
    dip_exact = integrals.basis_dipoles(res.basis)

    # Dipole matrix elements are computed by numerical integration

    # increase resolution of multicenter grid
    settings.radial_grid_factor = 20
    settings.lebedev_order = 41
    
    nbfs = res.basis.nbfs
    dip_numeric = np.zeros((3,nbfs,nbfs))
    for i in range(0, nbfs):
        # 
        orb_i = np.zeros(nbfs)
        orb_i[i] = 1.0
        # define wavefunctions for AO i
        def ao_i(x,y,z):
            return res.basis.wavefunction(orb_i, x,y,z)

        for j in range(0, nbfs):
            print "computing dipole matrix elements of AO pair %d-%d" % (i,j)
            # 
            orb_j = np.zeros(nbfs)
            orb_j[j] = 1.0
            
            # define wavefunctions for AO j
            def ao_j(x,y,z):
                return res.basis.wavefunction(orb_j, x,y,z)

            dip_numeric[:,i,j] = electronic_dipole(atomlist, ao_i, ao_j)

            print "  exact   = %s" % dip_exact[:,i,j]
            print "  numeric = %s" % dip_numeric[:,i,j]
            
    # comparison
    print "exact dipoles"
    print dip_exact

    print "numeric dipoles"
    print dip_numeric

    print "difference exact-numeric"
    print dip_exact-dip_numeric

    err = la.norm(dip_exact - dip_numeric)
    print "|Dip(exact)-Dip(numeric)|= %e" % err

    
if __name__ == "__main__":
    #test_integrals_h2()
    test_dipole_integrals()
