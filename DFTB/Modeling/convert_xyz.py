#!/usr/bin/env python
"""
Convert coordinates in an .xyz or .in file from bohr to angstrom
"""
from DFTB import XYZ
from DFTB import AtomicData

unit_conversion = {"b2a": (AtomicData.bohr_to_angs, "Angstrom"),
                   "a2b": (1.0/AtomicData.bohr_to_angs, "bohr")}

if __name__ == "__main__":
   import sys
   import os.path
   from optparse import OptionParser
   usage = "Usage: python %s <xyz input> <xyz output>\n" % os.path.basename(sys.argv[0])
   usage += "  converts the coordinates in the input file from bohr to angstrom or vice versa and writes the result into <xyz output>"

   parser = OptionParser(usage)
   parser.add_option("--convert", dest="convert", help="Type of conversion: 'b2a' - from bohr to Angstrom, 'a2b' - from Angstrom to bohr [default: %default]",  default="b2a")
   (opts, args) = parser.parse_args()
   if len(args) < 2:
      print usage
      exit(-1)
   xyz_in = args[0]
   xyz_out = args[1]
   if xyz_in[-2:] == "in":
      # initial conditions file ###.in
      coords, vels = XYZ.read_initial_conditions(xyz_in, units="")
      structures = [coords]
   else:
      structures = XYZ.read_xyz(xyz_in, units="")
   unit_fac, units = unit_conversion[opts.convert]
   for i, atomlist_in in enumerate(structures):
      vec_in = XYZ.atomlist2vector(atomlist_in)
      vec_out = unit_fac * vec_in
      atomlist_out = XYZ.vector2atomlist(vec_out, atomlist_in)
      if i == 0:
         mode = "w"
      else:
         mode = "a"
      XYZ.write_xyz(xyz_out, [atomlist_out], units="", mode=mode)
