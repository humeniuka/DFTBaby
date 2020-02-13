#!/usr/bin/env python
"""
enlarge unit cell by replicating it along all non-zero lattice vectors
"""

import DFTB.ForceField.PeriodicForceField as PFF

if __name__ == "__main__":
    import sys
    from optparse import OptionParser
    
    usage =  "Usage: python %s unit_cell.ff enlarged.ff  nmax \n" % sys.argv[0]
    usage += "   enlarge unit cell in 'unit_cell.ff' by replicating it `nmax` times it along all\n"
    usage += "   non-zero lattice vectors in positive and negative directions. The enlarged unit\n"
    usage += "   cell is written to 'enlarged.ff'.\n"

    parser = OptionParser(usage)
    (opts,args ) = parser.parse_args()

    if len(args) < 3:
        print usage
        exit(-1)

    ff_file = args[0]
    ff_enlarged_file = args[1]
    nmax = int(args[2])

    # read force field definition
    atomlist, atomtypes, partial_charges, lattice_vectors = PFF.read_force_field(ff_file)
    # enlarge unit cell
    atomlist_enlarged, lattice_vectors_enlarged = PFF.enlarged_unitcell(atomlist, lattice_vectors,
                                                                    nmax=nmax)
    # copy atom types and partial charges for replicas
    nunit = len(atomlist_enlarged) / len(atomlist)
    atomtypes_enlarged = []
    partial_charges_enlarged = []
    for i in range(0, nunit):
        atomtypes_enlarged += atomtypes
        partial_charges_enlarged += partial_charges

    # save enlarged unit cell
    PFF.write_force_field(ff_enlarged_file, atomlist_enlarged, atomtypes_enlarged,
                          partial_charges_enlarged, lattice_vectors_enlarged)
