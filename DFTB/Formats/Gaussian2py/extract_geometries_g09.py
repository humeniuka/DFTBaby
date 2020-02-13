#!/usr/bin/env python
#
# extract all geometries and SCF energies from Gaussian log-file
#

import sys
import os.path

from DFTB import XYZ
from DFTB.Formats.Gaussian2py import Gaussian

if len(sys.argv) < 3:
    print ""
    print "Usage: %s  log-file  xyz-file" % os.path.basename(sys.argv[0])
    print ""
    print "   extract all geometries and SCF energies from Gaussian log-file"
    print "   and write them to an xyz-file."
    print ""
    exit(-1)

# input file
log_file = sys.argv[1]
# output file
xyz_file = sys.argv[2]

atomlists = Gaussian.read_geometries_it(log_file)
energies = Gaussian.read_scf_energies_it(log_file)

titles = ["ENERGY= %f" % en for en in energies]

XYZ.write_xyz(xyz_file, atomlists, title=titles)

print "geometries written to '%s'" % xyz_file
