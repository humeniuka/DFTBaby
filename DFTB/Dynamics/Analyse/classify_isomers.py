#!/usr/bin/env python
"""
Given a set of isomer geometries, the time steps in a MD simulation are classified.
This gives a simplified representation of the dynamics  where the trajectory seems 
to hop between discrete isomer states. The assignment of isomer states is based on
the connectivity matrix: If a permutation of identical atoms exist such that the 
adjacency matrices of two geometries are identical the geometries are assumed to 
belong to the same isomer.

The classification is written to a table. The first column contains the time 
step, the second column constains the index of the isomer (starting from 0) 
which the current geometry resembles. If the molecule breaks up into several
fragments each fragment is classified separately so that the second column may
contain several indeces. 

Example:
========
Suppose the file 'dynamics.xyz' contains the geometries of a dynamics simulation.
The file 'isomers.xyz' contains several geometries which define the isomer states.
The command

   classify_isomers.py dynamics.xyz isomers.xyz --out_file=isomer_classification.dat

produces a table in 'isomer_classification.dat' that could look like this:

# GEOMETRY   ISOMER
 0           2        # In the first time step the geometry is in isomer state 2
 1           2
 ...
 100        5,6       # In the 100-th time step the geometry has fragmented. The first
                      # fragment resembles isomer 5 and the second one isomer 6.
 ...


"""
from DFTB.Analyse import MolecularGraph
from DFTB import XYZ

import sys
import optparse
import os.path

import numpy as np

if __name__ == "__main__":
    
    usage  = "python %s <dynamics.xyz file> <.xyz file with isomers>\n" % os.path.basename(sys.argv[0])
    usage += " For each time step in the trajectory the geometry is compared with the isomer geometries.\n"
    usage += " A classification of each geometry as one of the isomers is attempted based on the connectivity matrix.\n"
    usage += " If no isomer can be assigned the step is left out.\n"
    usage += " Type --help to see all options.\n"

    parser = optparse.OptionParser(usage)
    
    parser.add_option("--step", dest="step", help="Only every N-th geometry in the trajectory is extracted and classified [default: %default]", type=int, default=1)
    parser.add_option("--out_file", dest="out_file", help="The table with the isomer assignments for each geometry is written to this file [default: %default]", type=str, default="isomer_classification.dat")
    
    (opts, args) = parser.parse_args()
    if len(args) != 2:
        print usage
        exit(-1)

    dynamics_file = args[0]
    isomer_file = args[1]

    isomers = XYZ.read_xyz(isomer_file)
    # determined canonical connectivities of isomers
    isomer_connectivities = []
    for isomer in isomers:
        isomer_can, A_can = MolecularGraph.morgan_ordering(isomer, hydrogen_bonds=True)
        isomer_connectivities.append( A_can )
    # all isomers should have different connectivities
    n = len(isomers)
    for i in range(0, n):
        for j in range(i+1,n):
            if np.sum((isomer_connectivities[i] - isomer_connectivities[j])**2) == 0:
                print "WARNING: Isomer %d and %d in file '%s' have the same adjacency matrices!" % (i+1,j+1, isomer_file)
    
    fh = open(opts.out_file, "w")

    print>>fh, "# isomer indeces refer to the geometries in '%s'" % isomer_file
    print>>fh, "# GEOMETRY        ISOMER(S)"

    print "classify geometries by comparison with isomers"
    for i_step,atomlist in enumerate(XYZ.read_xyz_it(dynamics_file)):

        if i_step % opts.step != 0:
            continue

        # Each fragment has to be analyzed separately separately
        fragments = MolecularGraph.disconnected_fragments(atomlist, hydrogen_bonds=True)

        # The geometry is assigned the index of the isomer which it resembles.
        # If the molecule consists of several fragments each fragment is assigned
        # separately.

        if len(fragments) == 1:
            # disconnected_fragments sometimes messes up, ;(. If there is only one fragment,
            # it is better to use the original ordering of the atoms.
            fragments = [atomlist]
        classification = []
        for i_frag,fragment in enumerate(fragments):
            #print "Fragment %d" % i_frag
            # canonically ordered atoms and adjacency matrix
            atomlist_can, A_can = MolecularGraph.morgan_ordering(fragment, hydrogen_bonds=True)
        
            for i in range(0, n):
                # loop over isomers and compare with connectivity matrix of current geometry
                if A_can.shape != isomer_connectivities[i].shape:
                    # fragment and isomer do not contain the same number of atoms
                    continue
                if np.sum((A_can - isomer_connectivities[i])**2) == 0:
                    print "%d-th fragment of %d-th geometry classified as isomer %d" % (i_frag, i_step, i)                    
                    classification.append( i )
                    break
            else:
                print "%d-th fragment of %d-th geometry could not be classified" % (i_frag, i_step)
                continue
        # print table
        #
        #   GEOMETRY I        ISOMER(FRAGMENT 1) ISOMER(FRAGMENT 2), ....
        #
        if len(classification) > 0:
            print>>fh, "%d              %s" % (i_step, "   ".join(map(str, classification)))

    fh.close()
    print "Table with classifications written to %s" % opts.out_file
