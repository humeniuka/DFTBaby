#!/usr/bin/env python
"""
print xyz-coordinates in a format that can be pasted into python scripts
"""
from DFTB import XYZ
import pprint

if __name__ == "__main__":
   import sys
   usage = "Usage: python %s <xyz input>\n" % sys.argv[0]

   if len(sys.argv) < 2:
      print usage
      exit(-1)
   xyz = sys.argv[1]
   atomlist = XYZ.read_xyz(xyz)[0]

   pprint.pprint(atomlist)
