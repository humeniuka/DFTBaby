#!/usr/bin/env python
"""
Scan excitation energies along the geometries in 
a xyz-file
"""
from DFTB.PES import PotentialEnergySurfaces
from DFTB import XYZ, utils, AtomicData
from DFTB.Timer import GlobalTimer as T
from DFTB.Solver import ExcitedStatesNotConverged, ExcitedStatesError # exceptions

import numpy as np
import numpy.linalg as la

if __name__ == "__main__":
    import sys
    import os.path
    
    if len(sys.argv) < 4:
        print "Usage: %s <xyz-file> <nr. states> <dat-file with energies>" % os.path.basename(sys.argv[0])
        print "  scans the excitation energies along the geometries in the xyz-file"
        print "  and writes a column with the energy (in eV) for each state"
        print "  (including ground state) to the dat-file"
        print "  type --help to see all options"
        print "  to reduce the amount of output add the option --verbose=0"
        exit(-1)
    
    #
    xyz_file = sys.argv[1]    # path to xyz-file
    Nst = int(sys.argv[2])
    dat_file = sys.argv[3]    # output file
    # Read the geometry from the xyz-file
    atomlists = XYZ.read_xyz(xyz_file) 
    atomlist = atomlists[0]
    # read the charge of the molecule from the comment line in the xyz-file
    kwds = XYZ.extract_keywords_xyz(xyz_file)
    # initialize the TD-DFTB calculator
    pes = PotentialEnergySurfaces(atomlist, Nst=Nst, **kwds)

    fh = open(dat_file, "w")
    for i,atomlist in enumerate(atomlists):
        print "SCAN geometry %d of %d" % (i+1, len(atomlists))
        # convert geometry to a vector
        x = XYZ.atomlist2vector(atomlist)
        try:
            if Nst == 1:
                ens = pes.getEnergy_S0(x)
            else:
                ens = pes.getEnergies(x)
        except ExcitedStatesError as e:
            print "WARNING: %s" % e
            print "%d-th point is skipped" % i
            continue
        if i == 0:
            E0 = ens[0]
        # subtract ground state energy at first geometry
        ens -= E0
        # convert to eV
        ens *= AtomicData.hartree_to_eV
        print>>fh, "%d  " % i,
        for ei in ens:
            print>>fh, "%3.7f " % ei,
        print>>fh, ""
        fh.flush()
    fh.close()

    # timing
    print T
    
