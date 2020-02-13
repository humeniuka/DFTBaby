#!/usr/bin/env python

from DFTB.MolecularIntegrals.MulticenterIntegration import multicenter_inhomogeneous_schroedinger, multicenter_operation, atomlist2arrays, print_grid_summary
from DFTB.MolecularIntegrals import settings
from DFTB.MolecularIntegrals.BasissetFreeDFT import BasissetFreeDFT, density_func, effective_potential_func, regular_coulomb_func, residual_func, add_two_functions, spherical_average_residual_func, phaseshift_lstsq, radial_wave_func
from DFTB.MolecularIntegrals.Ints1e import overlap
from DFTB import AtomicData

import numpy as np

def lithium_cation_continuum(l, m, k):
    """
    compute continuum orbital in the electrostatic potential of the Li^+ core

    Parameters
    ----------
    l,m         :  angular quantum numbers of asymptotic solution
                   e.g.  l=0,m=0  s-orbital
                         l=1,m=+1 px-orbital
    k           :  length of wavevector in a.u., the energy of the 
                   continuum orbital is E=1/2 k^2
    """
    # Li^+ atom
    atomlist = [(3, (0.0, 0.0, 0.0))]
    charge = +1

    # choose resolution of multicenter grids for bound orbitals
    settings.radial_grid_factor = 20      # controls size of radial grid  
    settings.lebedev_order = 25          # controls size of angular grid
    # 1s core orbitals for Li+^ atom
    RDFT = BasissetFreeDFT(atomlist, None, charge=charge)
    # bound_orbitals = RDFT.getOrbitalGuess()
    Etot, bound_orbitals, orbital_energies = RDFT.solveKohnSham()

    
    # choose resolution of multicenter grids for continuum orbitals
    settings.radial_grid_factor = 120      # controls size of radial grid  
    settings.lebedev_order = 41          # controls size of angular grid
    # show number of radial and angular points in multicenter grid
    print_grid_summary(atomlist,
                       settings.lebedev_order, settings.radial_grid_factor)
    
    print "electron density..."
    # electron density of two electrons in the 1s core orbital
    rho = density_func(bound_orbitals)
    print "effective potential..."
    # potential energy for Li nucleus and 2 core electrons
    potential = effective_potential_func(atomlist, rho, None, nelec=2)

    def v0(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        return -1.0/r

    def v1(x,y,z):
        return potential(x,y,z) - v0(x,y,z)
    
    # The continuum orbital is specified by its energy and asymptotic
    # angular momentum (E,l,m)
    
    # energy of continuum orbital
    E = 0.5 * k**2  
    # angular quantum numbers of asymptotic solution
    assert abs(m) <= l

    print " "
    print "Asymptotic continuum wavefunction"
    print "================================="
    print "  energy           E= %e Hartree  ( %e eV )" % (E, E*AtomicData.hartree_to_eV)
    print "  wavevector       k= %e a.u." % k
    print "  angular moment   l= %d m= %+d" % (l,m)
    print " "
    
    # asymptotically correct solution for V0 = -1/r (hydrogen)
    Cf = regular_coulomb_func(E, charge, l, m, 0.0)
    phi0 = Cf
    
    # right-hand side of inhomogeneous Schroedinger equation
    def source(x,y,z):
        return -v1(x,y,z) * phi0(x,y,z)

    #
    #  solve  (H0 + V1 - E) dphi = - V1 phi0
    # for orbital correction dphi

    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)

    print "Schroedinger equation..."
    dphi = multicenter_inhomogeneous_schroedinger(potential, source,  E,
                                                  atomic_coordinates, atomic_numbers,
                                                  radial_grid_factor=settings.radial_grid_factor,
                                                  lebedev_order=settings.lebedev_order)

    # Combine asymptotically correct solution with correction 
    #   phi = phi0 + dphi
    phi = add_two_functions(atomlist, phi0, dphi, 1.0, 1.0)

    # residual for phi0    R0 = (H-E)phi0
    residual0 = residual_func(atomlist, phi0, potential, E)
    # residual for final solution  R = (H-E)phi
    residual  = residual_func(atomlist, phi, potential, E)

    # spherical average of residual function
    residual_avg = spherical_average_residual_func(atomlist[0], residual)

    # The phase shift is determined by matching the radial wavefunction
    # to a shifted and scaled Coulomb function at a number of radial
    # sampling points drawn from the interval [rmin, rmax].
    # On the one hand rmin < rmax should be chosen large enough,
    # so that the continuum orbital approaches its asymptotic form,
    # on the other hand rmax should be small enough that the accuracy
    # of the solution due to the sparse r-grid is still high enough.
    # A compromise has to be struck depending on the size of the radial grid.

    # The matching points are spread over several periods,
    # but not more than 30 bohr.
    wavelength = 2.0*np.pi / k
    print "wavelength = %e" % wavelength
    rmin = 70.0 
    rmax = rmin + max(10*wavelength, 30.0)
    Npts = 100

    # determine phase shift and scaling factor by a least square
    # fit the the regular Coulomb function
    scale, delta = phaseshift_lstsq(atomlist, phi, E, charge, l, m, rmin, rmax, Npts)

    print "scale factor (relative to Coulomb wave) = %s" % scale
    print "phase shift (relative to Coulomb wave) = %e " % delta

    # normalize wavefunction, so that 1/scale phi(x,y,z) approaches
    # asymptotically a phase-shifted Coulomb wave
    phi_norm = multicenter_operation([phi], lambda fs: fs[0]/scale, 
                                     atomic_coordinates, atomic_numbers,
                                     radial_grid_factor=settings.radial_grid_factor,
                                     lebedev_order=settings.lebedev_order)

    # The continuum orbital should be orthogonal to the bound
    # orbitals belonging to the same Hamiltonian. I think this
    # should come out correctly by default.
    print " "
    print " Overlaps between bound orbitals and continuum orbital"
    print " ====================================================="
    for ib,bound_orbital in enumerate(bound_orbitals):
        olap_bc = overlap(atomlist, bound_orbital, phi_norm)
        print "  <bound %d| continuum> = %e" % (ib+1, olap_bc)
    print ""

    # shifted regular coulomb function
    Cf_shift = regular_coulomb_func(E, charge, l, m, delta)
    
    # save radial wavefunctions and spherically averaged residual
    # radial part of Coulomb wave without phase shift
    phi0_rad = radial_wave_func(atomlist, phi0, l, m)
    # radial part of shifted Coulomb wave
    Cf_shift_rad = radial_wave_func(atomlist, Cf_shift, l, m)
    # radial part of scattering solution
    phi_norm_rad = radial_wave_func(atomlist, phi_norm, l, m)

    print ""
    print "# RADIAL_WAVEFUNCTIONS"
    print "# Asymptotic wavefunction:"
    print "#   charge            Z= %+d" % charge
    print "#   energy            E= %e  k= %e" % (E, k)
    print "#   angular momentum  l= %d m= %+d" % (l, m)
    print "#   phase shift       delta= %e rad" % delta
    print "# "
    print "#     R/bohr           Coulomb           Coulomb       radial wavefunction   spherical avg. residual"
    print "#                                        shifted           R_{l,m}(r)           <|(H-E)phi|^2>"
    import sys
    # write table to console
    r = np.linspace(1.0e-3, 100, 1000)
    data = np.vstack((r, phi0_rad(r), Cf_shift_rad(r), phi_rad(r), residual_avg(r))).transpose()
    np.savetxt(sys.stdout, data, fmt="    %+e    ")
    print "# END"
    print ""
    
if __name__ == "__main__":
    import sys
    import os.path

    args = sys.argv[1:]
    if len(args) < 3:
        usage = """
        Usage:
                %s  l   m  k

           compute the continuum orbital in the presence of the Li^+ core

        Parameters:
             l,m       -  integers, -l <= m <= l, angular quantum numbers
                          of asymptotic solution
             k         -  float, energy of continuum orbital is E = 1/2 k^2
        
        """ %  os.path.basename(sys.argv[0])

        print usage
        exit(-1)

    l = int(args[0])
    m = int(args[1])
    k = float(args[2])

    lithium_cation_continuum(l, m, k)

    
