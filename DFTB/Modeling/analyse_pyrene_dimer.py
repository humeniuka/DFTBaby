#!/usr/bin/env python
"""
computes the horizontal (Rx and Ry) and vertical (Rz) distances between two pyrene monomers
for each geometry in a trajectory
"""

from DFTB.Modeling import MinimalEnclosingBox
from DFTB import utils
from DFTB import QMMM

def getQMatoms(qmmm_partitioning=None):
    """
    QM/MM.qmmm_partitioning: If the file contains additional molecules (e.g. MM molecules), you should provide a list of the indeces belonging to the pyrene dimer of interest.
    """
    if qmmm_partitioning == None:
        return None
    if "-" in qmmm_partitioning:
        # index list contains ranges such as "9-14"
        indeces = QMMM.parseAtomTags(inner_indeces)
    else:
        indeces = list(eval(qmmm_partitioning))
    # start at 0
    indeces = [i-1 for i in indeces]
    return indeces
    
if __name__ == "__main__":
    import sys
    usage  = "Usage: python %s <.xyz file> <.dat output file>\n" % sys.argv[0]
    usage += "  The xyz-file should contain the geometries of the two pyrene units.\n"
    usage += "  A table with Rx,Ry and Rz will be written to the output file.\n"
    usage += "  Type --help to see all options.\n"

    parser = utils.OptionParserFuncWrapper([getQMatoms], usage)

    (opts,args) = parser.parse_args(getQMatoms)
    if len(args) < 2:
        print usage
        exit(-1)

    xyz_traj_file = args[0]
    out_file = args[1]

    MinimalEnclosingBox.pyrene_dimer_analysis(xyz_traj_file, out_file, qmmm_partitioning=getQMatoms(**opts))


        
