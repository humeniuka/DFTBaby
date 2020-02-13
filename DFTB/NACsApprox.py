#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
approximate non-adiabatic coupling vectors are obtained from
Mulliken transition charges
"""

import numpy as np
import numpy.linalg as la

from DFTB import XYZ, AtomicData


def coupling_vector(atomlist, transition_charges, en_excI, nuclear_charge="screened"): #"full"):
    """
    approximate NACs for transition from S_I to S_0.

    compute approximate non-adiabatic coupling vector from transition charges
    {q_A}_{A=1,...,Nat} and cartesian positions of atoms {R_A}_{A=1,...,Nat}
    which have nuclear charges {Z_A}_{A=1,...,Nat}:

                  < S_0 | d(V_{nuc-elec})/dR_A | S_I >
          NAC_A = -------------------------------------
                           E(n) - E(0)
        
                      /        R_A - r
                = Z_A | d^3r ----------   rho(r)
                      /      |R_A - r|^3     I,0

    rho_{I,0}(r) is the transition density between the ground state S_0 and the
    excited state S_I.
    
    Approximating the transition density by a sum of delta-functions centered
    on the atoms,

          rho_{I,0}(r) = sum_A  q_A delta(R_A - r)

    the non-adiabatic coupling vector becomes approximately

                     Z_A                  R_A - R_B
          NAC_A = ----------  sum   q_B -------------
                   E(I)-E(0)   B!=A      |R_A - R_B|^3


    Parameters
    ----------
    atomlist           : list of tuples (Z,[x,y,z]) with atomic numbers and atomic coordinates
    transition_charges : numpy array with Mulliken transition charges for state I
    en_excI            : excitation energy of state I (in Hartree)

    Optional
    --------
    nuclear_charge     : In the semiempirical approximation the nuclear charge is either
                         taken as the full charge ('full') or as the negative of the number
                         of valence electrons ('screened')

    Returns
    -------
    nac        : numpy array of shape (3,Nat),
                 non-adiabatic coupling vector for transition from I-th excited
                 state to ground state, 
                 nac[:,A] is the vector on the A-th atom
    """
    Nat = len(atomlist)
    
    nac = np.zeros((3,Nat))
    for A,(ZA,posA) in enumerate(atomlist):
        for B,(qB,(ZB,posB)) in enumerate(zip(transition_charges, atomlist)):
            if A == B:
                continue
            rBA = np.array(posA) - np.array(posB)

            # Since only valence electrons are considered, we maybe need to use the 
            # charges of the core (atomic number - nr. cor electrons) instead
            # of the full nuclear charges? But I am not sure.
            if nuclear_charge == "full":
                nac[:,A] += ZA/en_excI * qB * rBA / la.norm(rBA)**3
            elif nuclear_charge == "screened":
                ZcoreA = AtomicData.valence_electrons[AtomicData.atom_names[ZA-1]]
                nac[:,A] += ZcoreA/en_excI * qB * rBA / la.norm(rBA)**3
            else:
                raise ValueError("'nuclear_charge' option should be 'full' or 'screened'")
            
    return nac

        
    
if __name__ == "__main__":
    import sys
    import os.path
    from optparse import OptionParser
    
    usage = """
       Usage: %s    <.chg transition charges>  <.nac output file>  <energy diff. in eV>

    compute the non-adiabatic coupling vector 

                                < S0 | ( d/dR V_{nuc-elec} ) | S1 >
            < S0 | d/dR S1 > = ------------------------------------
                                         E(1) - E(0)

    approximately using transition charges.

    The energy difference in the denominator has to be provided as the 3rd argument (in eV).

    Input files:
      - file with molecular geometry in xyz format and transition charges in 4th column
    Output file:
      - file with coupling vector as a matrix of shape (num.atoms,3)

    """ % os.path.basename(sys.argv[0])

    parser = OptionParser(usage)
    parser.add_option("--nuclear_charge", dest="nuclear_charge", type=str, default="full", help="In the semiempirical approximation the nuclear charge Z is either taken as the full charge ('full') or as the negative of the number of valence electrons ('screened') [default: %default]")
    (opts, args) = parser.parse_args()
    
    if len(args) < 3:
        parser.print_help()
        exit(-1)
        
    # input files
    chg_file = args[0]
    # output file
    nac_file = args[1]
    # S0-S1 excitation energy (in eV)
    en_diff = float(args[2]) / AtomicData.hartree_to_eV
    
    # load molecular geometry and transition charges
    atomlist, transition_charges = XYZ.read_charges(chg_file)

    # compute NAC vector
    nac = coupling_vector(atomlist, transition_charges, en_diff, nuclear_charge=opts.nuclear_charge)

    # and save it to a file
    fh_nac = open(nac_file, "w")
    # write header for NAC vector
    print>>fh_nac, "# non-adiabatic coupling vector ( in bohr^-1 ) "
    print>>fh_nac, "#    X             Y             Z             "
    
    # print non-adiabatic coupling vector as (Nat,3) matrix
    Nat = len(atomlist)
    np.savetxt(fh_nac, nac.transpose(), fmt="%+e")

    print "Approximate non-adiabatic coupling vector written to %s" % nac_file
