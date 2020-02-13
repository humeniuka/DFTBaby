#!/usr/bin/env python
"""
scan the value of an internal coordinate (bond length, angle, dihedral)
in order to create a sequence of geometries 
"""

from DFTB.Optimize.InternalCoords import InternalValenceCoords
from DFTB import XYZ, AtomicData

import numpy as np

if __name__ == "__main__":
    import os.path
    import sys
    from optparse import OptionParser

    usage = """
       %s  initial.xyz  scan.xyz

    scans the value of an internal coordinate (bond length, angle or dihedral)
    starting with the initial geometry in `initial.xyz` and writes the 
    displaced geometries to `scan.xyz`.

    The scan coordinate and the number of steps and increment have to specified
    using the --scan=... option.

Example:

       %s  initial.xyz  scan.xyz  --scan="(1,2, 10, 0.1)"

    scans the distance between atoms 1 and 2 in 10 steps of 0.1 Angstroms.

    """ % (os.path.basename(sys.argv[0]), os.path.basename(sys.argv[0]))

    parser = OptionParser(usage)
    parser.add_option("--scan", dest="scan", type=str, default="[]", help="Define internal coordinate to be scanned. The coordinate is incremented from its initial value by `nsteps` steps of size `incr`. The internal coordinate whose value should be changed is specified by 2, 3 or 4 atom indicies followed by the number of steps and the increment, which is in Angstrom for bond lengths and degrees for angles. The format is \"(I,J, nsteps, incr)\" for scanning a bond length, \"(I,J,K, nsteps, incr)\" for scanning a valence angle and \"(I,J,K,L, nsteps, incr)\" for scanning a torsion.  For example \"(1,2, 5, 0.1)\" will scan the bond between atom 1 and 2 in 5 steps of 0.1 Angstrom and \"(1,2,3, 9, 10.0)\" will scan the angle 1-2-3 in 9 steps of 10.0 degrees. Dihedral angles are limited to the range [0,180], so if the angle is close to 180 degrees a negative increment should be used, otherwise the scan will stop at 180 degrees. [default: %default]")
    parser.add_option("--freeze", dest="freeze", type=str, default="[]", help="Freeze internal coordinates. Modifying one internal coordinate usually affects other coordinates, too. The internal coordinates that should be kept at their current value during the scan are specified as a list of tuples of atom indices (starting at 1). Each tuple may contain 2, 3 or 4 atom indices, (I,J) - bond between atoms I and J, (I,J,K) - valence angle I-J-K, (I,J,K,L) - dihedral angle between bonds I-J, J-K and K-L. For example \"[(1,2), (4,5,6)]\" freezes the bond between atoms 1 and 2 and the angle 4-5-6. The atom indices do not necessarily have to correspond to a 'physical' bond, angle or dihedral. So, for instance, you can also freeze the distance between two atoms that are not connected. [default: %default]")
    parser.add_option("--verbose", dest="verbose", type=int, default=0, help="Request additional output about internal coordinate transformations, 0 (no output), > 0 (lots of output)")
    
    (opts, args) = parser.parse_args()
    if len(args) < 2:
        parser.print_help()
        exit(-1)

    xyz_in = args[0]
    xyz_out = args[1]
    
    # read initial geometry
    atomlist0 = XYZ.read_xyz(xyz_in)[-1]
    x0 = XYZ.atomlist2vector(atomlist0)

    # freezing of internal coordinates
    freeze = []
    for IJKL in eval(opts.freeze):
        # Indices on the command line start at 1, but internally
        # indices starting at 0 are used.
        IJKL = tuple([I-1 for I in IJKL])
        freeze.append(IJKL)

    # scan coordinate
    scan = eval(opts.scan)

    # Determine type of internal coordinate and convert the units
    # accordingly
    if len(scan) == 4:
        # scan bond length
        I,J, nsteps, incr = scan
        # convert increment from Angstrom to bohr
        incr /= AtomicData.bohr_to_angs
        IJKL = (I,J)
        coord_name = "BOND(%d-%d)" % (I,J)
    elif len(scan) == 5:
        # scan valence angle
        I,J,K, nsteps, incr = scan
        # convert angle from degrees to radians
        incr *= np.pi/180.0
        IJKL = (I,J,K)
        coord_name = "ANGLE(%d-%d-%d)" % (I,J,K)
    elif len(scan) == 6:
        # scan dihedral angle
        I,J,K,L, nsteps, incr = scan
        # convert angle from degrees to radians
        incr *= np.pi/180.0
        IJKL = (I,J,K,L)
        coord_name = "DIHEDRAL(%d-%d-%d-%d)" % (I,J,K,L)
    else:
        raise ValueError("Format of scan '%s' not understood!" % scan)

    # shift indices to programmer's style (starting at 0)
    IJKL = map(lambda I:I-1, IJKL)
    
    IC = InternalValenceCoords(atomlist0, freeze=freeze, verbose=opts.verbose)
    
    # cartesian coordinates along the scan
    scan_geometries = [atomlist0]
    # value of internal coordinate IJKL along the scan
    val0 = IC.coordinate_value(x0, IJKL)
    scan_coords = [val0]
    
    x1 = x0
    for i in range(0, nsteps):
        # cartesian coordinates at displaced geometry
        x1 = IC.internal_step(x1, IJKL, incr)
        # new value of internal coordinate, should be approximately
        #  val0 + incr
        val1 = IC.coordinate_value(x1, IJKL)

        atomlist1 = XYZ.vector2atomlist(x1, atomlist0)
        scan_geometries.append(atomlist1)
        scan_coords.append(val1)
        print "step %d  %s = %10.6f" % (i+1, coord_name, val1)
        
    XYZ.write_xyz(xyz_out, scan_geometries, mode="w")
    print "displaced geometries were saved to '%s'" % xyz_out
