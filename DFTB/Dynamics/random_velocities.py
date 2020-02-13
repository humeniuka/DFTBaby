#!/usr/bin/env python
"""
initial conditions with random velocities in the format expected by the surface hopping program
"""
from DFTB import AtomicData
from DFTB import XYZ
from DFTB.Dynamics.HarmonicApproximation import dynamics_in_format

import sys
import os.path
import numpy as np
from optparse import OptionParser

if __name__ == "__main__":
    usage = "Usage: python %s <.xyz file> <.in output file>\n" % os.path.basename(sys.argv[0])
    usage += "  convert the xyz-file to initial conditions with random velocities that correspond to a certain temperature.\n"
    usage += "  Type --help to see all options.\n"
    parser = OptionParser(usage)
    parser.add_option("--temperature", dest="temperature", default=0.0, type=float, help="Set the temperature T (in Kelvin) of the molecule such that Ekin = 1/2 * nvib * kB * T  [default: %default]")
    (opts, args) = parser.parse_args()
    if len(args) < 2:
        print usage
        exit(-1)
    xyz_file = args[0]
    in_file = args[1]
    # temperature
    T = opts.temperature
    atomlist = XYZ.read_xyz(xyz_file)[-1]
    masses = np.array(AtomicData.atomlist2masses(atomlist))
    nat = len(atomlist)
    # initial positions are read from file
    q = XYZ.atomlist2vector(atomlist)
    # draw random numbers for initial momenta
    p = np.random.rand(3*nat)
    # scale momenta so that the kinetic energy equals 1/2 nvib * kB * T
    if nat < 2:
        raise Exception("Since vibrational and rotational degrees are removed, a single atom cannot have any meaningful temperature.")
    elif nat == 2:
        nvib = 1
    elif nat > 2:
        nvib = 3*nat-6
    print "number of atoms: %s" % nat
    print "number of vibrational modes: %s" % nvib
    ekin0 = np.sum( p**2 / (2.0 * masses) )
    T0 = ekin0 / (0.5 * nvib * AtomicData.kBoltzmann)
    # scale momenta
    p *= np.sqrt(T/T0)
    # check 
    ekin = np.sum( p**2 / (2.0 * masses) )
    Tnew = ekin / (0.5 * nvib * AtomicData.kBoltzmann)
    print "temperature: %.5f K" % Tnew
    
    txt = dynamics_in_format(atomlist, q, p, "")
    fh = open(in_file, "w")
    fh.write(txt)
    fh.close()
    print "Wrote initial conditions to %s" % in_file

