#!/usr/bin/env python
"""
convert partial charges in Turbomole format to xyz-format
"""
from DFTB.Formats import Turbomole
from DFTB import XYZ
import sys

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: python %s <TM file> <xyz output file>" % sys.argv[0]
        sys.exit(-1)
    tmfile = sys.argv[1]
    xyzfile = sys.argv[2]

    Data = Turbomole.parseTurbomole(tmfile)
    point_charges = Data["point_charges"]
    XYZ.write_xyz(xyzfile, [point_charges], title="Converted from Turbomole %s" % tmfile)

