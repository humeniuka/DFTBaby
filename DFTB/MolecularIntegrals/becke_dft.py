#!/usr/bin/env python
"""
Becke's basis set free DFT calculation
"""

from DFTB import XYZ
from DFTB.SlaterKoster import XCFunctionals
from DFTB.MolecularIntegrals.BasissetFreeDFT import BasissetFreeDFT
from DFTB.MolecularIntegrals import settings

import os.path
import sys

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: %s xyz-file  rfac  Lmax" % os.path.basename(sys.argv[0])
        print "  compute DFT ground state"
        exit(-1)

    # load geometry
    xyz_file = sys.argv[1]
    atomlist = XYZ.read_xyz(xyz_file)[0]
    # set resolution of multicenter grid
    rfac = int(sys.argv[2])
    Lmax = int(sys.argv[3])

    settings.radial_grid_factor = rfac      # controls size of radial grid 
    settings.lebedev_order = Lmax           # controls size of angular grid
    
    charge = 0
    
    # PBE functional
    xc = XCFunctionals.libXCFunctional('gga_x_pbe', 'gga_c_pbe')
    
    
    RDFT = BasissetFreeDFT(atomlist, xc, charge=charge)
    RDFT.solveKohnSham_new(thresh=1.0e-8)

    
