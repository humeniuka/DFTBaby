#!/usr/bin/env python
"""
convert repulsive potential of hotbit or DFTB+ into format
directly readable with python.
"""
import sys
from optparse import OptionParser

from RepulsivePotential import read_repulsive_potential, write_repulsive_potential

usage  = "Usage: %s <.par or.skf reppot-file1> <.py reppot-file2>\n" % sys.argv[0]
usage += "Read first repulsive potential in .par or .skf format and write it\n"
usage += "to the second file as python data.\n"
usage += "Type --help to see all options.\n"

parser = OptionParser(usage)
parser.add_option("--scaling", dest="scaling", type=float, default=1.0, help="Scale the r-axis of the repulsive potential with a factor S to make it more repulsive (S > 1.0) or less repulsive (S < 1.0) [default: %default]")

(opts, args) = parser.parse_args()
if len(args) < 2:
    print usage
    exit(-1)

print "CAUTION: You cannot reuse repulsive potentials unless they are based on the same parametrization of the electronic part (same r0, same Hubbard U, same xc potential) as the Slater-Koster tables."
mod = read_repulsive_potential(sys.argv[1])
# scale d-axis
mod.d *= opts.scaling
#
title = "converted from %s" % sys.argv[1]
if opts.scaling != 1.0:
    title += " and scaled, Vrep'(r) =  Vrep(%s*r)" % opts.scaling
write_repulsive_potential(mod, sys.argv[2], title=title)

