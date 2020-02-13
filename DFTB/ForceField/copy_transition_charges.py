#!/usr/bin/env python
"""
copy the transition charges and excitation energies of the monomer to the monomeric units in the stack or tube.
"""

from DFTB import XYZ
import DFTB.ForceField.PeriodicForceField as PFF
from DFTB.ForceField.fit_magnetic_dipoles import local_axes

import numpy as np

import sys
import os.path

def rotate_magnetic_dipoles(atomlist, transition_charges):
    """
    rotate magnetic dipoles so that they are defined relative to the global coordinate frame
    and not the local molecular frame (principle axes of rotation)
    """
    nat = len(atomlist)
    R = np.zeros((3,nat))
    for i,(Zi,posi) in enumerate(atomlist):
        R[:,i] = posi
    # 
    axes = local_axes(R)
    transition_charges = np.array(transition_charges)
    _nat, nst4 = transition_charges.shape
    assert _nat == nat
    nst = nst4/4
    assert nst*4 == nst4

    transition_charges_rot = np.zeros((nat,4*nst))
    for i in range(0, nat):
        for j in range(0, nst): 
            mvec = transition_charges[i,4*j+1:4*j+4]
            mvec_rot = np.dot(axes, mvec)
            transition_charges_rot[i,4*j] = transition_charges[i,4*j]  # copy charge
            transition_charges_rot[i,4*j+1:4*j+4] = mvec_rot # rotated magnetic dipole
            
    return transition_charges_rot.tolist()
            
if __name__ == "__main__":
    from optparse import OptionParser
    
    usage  = "Usage: %s <monomer.chromo file> <oligomer.xyz file> <oligomer.chromo file>\n" % os.path.basename(sys.argv[0])
    usage += "  copy the transition charges and magnetic dipoles  of the monomer (in monomer.chromo)\n"
    usage += "  to all monomeric units in oligomer.xyz and save the result to <oligomer.chromo>\n"
    
    parser = OptionParser(usage)
    parser.add_option("--rotate_dipoles", dest="rotate_dipoles", type=int, default=0, help="Rotate magnetic dipoles from molecular frame into global frame, so that they can be visualized with the correct orientation [default: %default].")

    (opts, args) = parser.parse_args()
    if len(args) < 3:
        print usage
        exit(-1)

    # input files
    monomer_chromo = args[0]
    oligomer_xyz = args[1]
    # output file
    oligomer_chromo = args[2]

    # read monomer 
    atomlist_mono = XYZ.read_xyz(monomer_chromo)[0]
    for chromophore in PFF.read_transition_charges(monomer_chromo):
        indeces, excitation_energies, transition_charges = chromophore
        break
    nat = len(atomlist_mono)
    
    # read oligomer
    atomlist_oligo = XYZ.read_xyz(oligomer_xyz)[0]

    # number of units
    nunit = len(atomlist_oligo) / len(atomlist_mono)
    assert nunit * len(atomlist_mono) == len(atomlist_oligo), "Number of atoms must be multiple of number of atoms in monomer"

    mode = "w"
    for i in range(0, nunit):
        atomlist = atomlist_oligo[i*nat:(i+1)*nat]
        if opts.rotate_dipoles > 0:
            transition_charges_copy = rotate_magnetic_dipoles(atomlist, transition_charges)
        else:
            transition_charges_copy = transition_charges
        PFF.write_transition_charges(oligomer_chromo, atomlist,
                                     indeces, excitation_energies, transition_charges_copy,
                                     index_offset=i*nat,
                                     mode=mode)
        mode = "a"
    
