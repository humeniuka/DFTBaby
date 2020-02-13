#!/usr/bin/env python
"""
create intermediate geometries between two points by linear
interpolation in cartesian or internal coordinates
"""
from DFTB.Optimize.InternalCoords import InternalValenceCoords
from DFTB import XYZ

import numpy as np


if __name__ == "__main__":
    import sys
    import os.path
    from optparse import OptionParser
    
    usage = """

        %s  initial.xyz  final.xyz  N  path.xyz

    The geometry is interpolated linearly between the initial
    and the final structures using N-2 intermediate points.
    The interpolated geometries are written to 'path.xyz'."

    """ % os.path.basename(sys.argv[0])
    
    parser = OptionParser(usage)
    parser.add_option("--coord_system", dest="coord_system", help="The geometry can be interpolated in 'cartesian' or 'internal' coordinates. [default: %default]", type=str, default="cartesian")
    parser.add_option("--explicit_bonds", dest="explicit_bonds", help="Inserts artificial bonds between pairs of atoms. The bonds are specified as a list of tuples (I,J) of atom indices (starting at 1). default: %default]", type=str, default="[]")

    (opts, args) = parser.parse_args()

    if len(args) < 4:
        parser.print_help()
        exit(-1)
    
    xyz0 = args[0]
    xyz1 = args[1]
    N = int(args[2])
    xyz_interp = args[3]
    
    # read initial geometry
    atomlist0 = XYZ.read_xyz(xyz0)[0]
    x0 = XYZ.atomlist2vector(atomlist0)
    
    # read final geometry
    atomlist1 = XYZ.read_xyz(xyz1)[0]
    x1 = XYZ.atomlist2vector(atomlist1)

    # interpolation parameter
    rs = np.linspace(0.0, 1.0, N)
    
    geometries_interp = []

    if opts.coord_system == "cartesian":
        for r in rs:
            xr = x0 + r*(x1-x0)
            geometries_interp.append( XYZ.vector2atomlist(xr, atomlist0) )
    elif opts.coord_system == "internal":
        # explicit bonds, shift atom indices from 1- to 0-indexing
        explicit_bonds = [(I-1,J-1) for (I,J) in eval(opts.explicit_bonds)]
        
        IC = InternalValenceCoords(atomlist0, explicit_bonds=explicit_bonds)
        # initial and final geometry in internal coordinates
        q1 = IC.cartesian2internal(x1)
        q0 = IC.cartesian2internal(x0)
        for r in rs:
            qr = q0 + r*(q1-q0)
            xr = IC.internal2cartesian(qr)
            geometries_interp.append( XYZ.vector2atomlist(xr, atomlist0) )      
    else:
        raise ValueError("Coordinate system '%s' not understood, valid options are 'internal' and 'cartesian'" % opts.coord_system)
            
    XYZ.write_xyz(xyz_interp, geometries_interp)
    print "Interpolated geometries written to %s" % xyz_interp

