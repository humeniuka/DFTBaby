#!/usr/bin/env python
"""
copy the atom type assignments in the monomer to the monomeric units in the stack or tube.
"""

from DFTB import XYZ
import DFTB.ForceField.PeriodicForceField as PFF

import sys
import os.path

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: %s <monomer.ff file> <oligomer.xyz file> <oligomer.ff file>" % os.path.basename(sys.argv[0])
        print "  copy the atom type assignments of the monomer (in monomer.ff) to all"
        print "  monomeric units in oligomer.xyz and save the result to <oligomer.ff>"
        print ""
        exit(-1)

    # input files
    monomer_ff = sys.argv[1]
    oligomer_xyz = sys.argv[2]
    # output file
    oligomer_ff = sys.argv[3]

    atomlist_mono, atomtypes_mono, partial_charges_mono, lattice_vectors = PFF.read_force_field(monomer_ff)
    nat = len(atomlist_mono)
    
    atomlist_oligo = XYZ.read_xyz(oligomer_xyz)[0]

    # number of units
    nunit = len(atomlist_oligo) / len(atomlist_mono)
    assert nunit * len(atomlist_mono) == len(atomlist_oligo), "Number of atoms must be multiple of number of atoms in monomer"

    atomtypes_oligo = []
    partial_charges_oligo = []
    for i in range(0, nunit):
        atomtypes_oligo += atomtypes_mono
        partial_charges_oligo += partial_charges_mono

    PFF.write_force_field(oligomer_ff, atomlist_oligo, atomtypes_oligo, partial_charges_oligo, lattice_vectors=[])

    
