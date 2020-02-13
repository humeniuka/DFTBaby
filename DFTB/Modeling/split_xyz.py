#!/usr/bin/env python
"""
split an xyz file that contains many structures into separate files
"""
from DFTB import XYZ

if __name__ == "__main__":
    import sys
    import optparse
    usage =  "Usage: python %s <xyz file with many structures> <output prefix>\n" % sys.argv[0]
    usage += "  splits the input file into output prefix_####.xyz"

    parser = optparse.OptionParser(usage)
    parser.add_option("--step", dest="step", help="Only every N-th geometry is extracted [default: %default]", type=int, default=1)
    (opts, args) = parser.parse_args()
    
    if len(args) < 2:
        print usage
        exit(-1)
    xyz_file = sys.argv[1]
    prefix = sys.argv[2]
    structures = XYZ.read_xyz(xyz_file)
    for i, atomlist in enumerate(structures):
        if (i+1) % opts.step == 0 or i == 0:
            XYZ.write_xyz(prefix + "%.4d.xyz" % (i+1), [atomlist])
