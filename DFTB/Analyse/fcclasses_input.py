#!/usr/bin/env python
"""
generates auxiliary formatted input file F01 or F02 for FCclasses with information
about the state (geometry, vibrational modes and frequencies) based on semiempirical 
calculations with MOPAC-7 or DFTBaby.
"""
import numpy as np

from DFTB import XYZ, AtomicData
from DFTB.Dynamics import HarmonicApproximation

if __name__ == "__main__":
    import sys
    import os.path
    
    usage = """
    Usage:

        python %s   opt.xyz    hessian.dat    state.dat

    generate auxiliary state file for FCclasses from semiemperical frequency calculation
    
    Input Files:
       opt.xyz          -  contains optimized geometry of ground or excited state
       hessian.dat      -  contains matrix with Hessian of energy of state of interest
    Output File:
       state.dat        -  file with input for state1 or state2 for FCclasses
                           (see section 6.2 of FCclasses-2.1 manual)
       masses.dat       -  column vector with atomic masses in amu

    """ % os.path.basename(sys.argv[0])

    if len(sys.argv) < 4:
        print usage
        exit(-1)
        
    args = sys.argv[1:]
    
    # input files
    # ... optimized geometry
    xyz_file = args[0]
    # ... Hessian of energy
    hessian_file = args[1]
    # output file
    state_file = args[2]
    
    # load minimum geometry and Hessian 
    atomlist = XYZ.read_xyz(xyz_file)[-1]
    hess = np.loadtxt(hessian_file)

    print "optimized geometry read from '%s'" % xyz_file
    print "Hessian read from '%s'" % hessian_file
    
    # compute normal modes and frequencies
    xopt = XYZ.atomlist2vector(atomlist)
    masses = AtomicData.atomlist2masses(atomlist)
    # setting zero_threshold to -0.1 makes sure that all frequencies are returned
    # even if some of them are imaginary
    vib_freq, vib_modes = HarmonicApproximation.vibrational_analysis(xopt, hess, masses, \
                                                                     zero_threshold=-0.1, is_molecule=True)

    # sort frequencies in ascending order
    vib_freq = vib_freq.real
    sort_index = np.argsort(abs(vib_freq)**2)
    vib_freq = vib_freq[sort_index]
    vib_modes = vib_modes[:,sort_index]
    
    # remove lowest 6 vibrations (3 translation + 3 rotation) assuming the molecule
    # is not linear
    nat = len(atomlist)
    nvib = 3*nat-6
    # Modes which have frequencies exactly equal to 0 (or with imaginary parts)
    # are already removed inside `vibrational_analysis()`. Another `nzero`
    # low frequency modes have to be removed, so that only 3*nat-6 vibrational frequencies remain
    nzero = len(vib_freq) - nvib
    assert nzero >= 0
    vib_freq = vib_freq[nzero:]
    vib_modes = vib_modes[:,nzero:]

    # Vibrational modes are in mass weighted coordinates
    #  q_i = sqrt(m_i) * (delta x_i)
    #
    sqiM = np.diag(1.0/np.sqrt(masses))
    # FCclasses needes cartesian displacement vectors (delta x_i)
    vib_modes_cart = np.dot(sqiM, vib_modes)
    
    # write state file
    fh = open(state_file, "w")

    
    # The auxiliary state file for FCclasses (F01 or F02) contains
    # one long vector of floats:
    #
    
    # 1) cartesian coordinates in Angstrom
    #    GEO(I), I=1,3N
    for Zi,posi in atomlist:
        # convert from bohr to Angstrom
        posi = np.array(posi) * AtomicData.bohr_to_angs
        # x,y,z coordinates
        print>>fh, "  %+20.15e" % posi[0]
        print>>fh, "  %+20.15e" % posi[1]
        print>>fh, "  %+20.15e" % posi[2]
        
    # 2) normal mode displacements in cartesian coordinates
    #    T(I,J), J=3*N-6, I=1,3*N, J is the faster index
    for i in range(0,3*nat):
        for j in range(0,nvib):
            print>>fh, "   %+20.15e" % vib_modes_cart[i,j]

    # 3) vibrational frequencies in cm^-1
    #    G(J), J=1,3*N-6
    # convert from Hartree to cm^-1
    vib_freq *= AtomicData.hartree_to_wavenumbers
    for j in range(0, nvib):
        print>>fh, "%20.15f" % vib_freq[j]

    fh.close()

    print "State file written to '%s'" % state_file

    masses = np.array(AtomicData.atomlist2masses(atomlist))
    masses_amu = masses[0::3] * AtomicData.aumass2amu

    np.savetxt("masses.dat", masses_amu, fmt="     %15.8f  ")

    print "Masses in amu written to 'masses.dat'"
    
    
        
