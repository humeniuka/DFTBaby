#!/usr/bin/env python
"""
This script allows to keep or remove certain fragments from an XYZ file. You can specify the atoms
of the fragments that should be retained as a string in the canonical order according to the
Morgan algorithm. 
Example:

    filter_fragments.py pyrene_and_stuff.xyz hcchcchchcchchchcchchcchcc --filter_mode=keep

will remove all atoms that are not part of a pyrene molecule, whose string representation is 
'hcchcchchcchchchcchchcchcc'

"""

from DFTB.Analyse import MolecularGraph as MG
from DFTB import XYZ

def filter_fragments(atomlist_full, selected_fragments=[], filter_mode="keep"):
    """

    Parameters:
    ===========
    atomlist_full: geometry that should be filtered
    selected_fragments: list of fragment labels that should be kept or deleted
    filter_mode: 'keep' - keep selected fragments, 'delete' - delete selected fragments

    Returns:
    ========
    atomlist_filtered: geometry with atoms that passed the filter
    """
    #
    fragments = MG.disconnected_fragments(atomlist_full)
    fragment_labels = [MG.identifier( *MG.morgan_ordering(atomlist) ) for atomlist in fragments]
    atomlist_filtered = []
    #
    unique_labels = set(fragment_labels)
    print "The following fragments were found:"
    print "==================================="
    for flabel in unique_labels:
        print " %s" % flabel
    print ""
    #
    for frag,flabel in zip(fragments, fragment_labels):
        if flabel in selected_fragments:
            if filter_mode == "keep":
                atomlist_filtered += frag
        else:
            if filter_mode == "delete":
                atomlist_filtered += frag

    return atomlist_filtered


if __name__ == "__main__":
    import sys
    import optparse
    usage = "Usage: python %s <xyz-file> <list of fragment names>\n" % sys.argv[0]
    usage += "  filters geometry so that only the listed fragments are retained.\n"
    usage += "  To see a list of all available fragments run the program without the list parameter.\n"
    usage += "  Type --help to see all options."
    parser = optparse.OptionParser(usage)
    parser.add_option("--filter_mode", dest="filter_mode", help="Determines whether the filter should 'keep' or 'delete' the selected fragments [default: %default]", default="keep")
    parser.add_option("--out_xyz", dest="out_xyz", help="Filtered geometry is written to this file [default: %default]", default="filtered.xyz")

    (opts,args) = parser.parse_args()
    if len(args) < 1:
        print usage
        exit(-1)
    xyz_in = args[0]
    selected_fragments = args[1:]
    
    atomlist_full = XYZ.read_xyz(xyz_in)[-1]
    atomlist_filtered = filter_fragments(atomlist_full, selected_fragments, filter_mode=opts.filter_mode)

    if len(selected_fragments) > 0:
        XYZ.write_xyz(opts.out_xyz, [atomlist_filtered])
        print "filtered geometry written to %s" % opts.out_xyz
