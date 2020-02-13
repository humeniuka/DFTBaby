#!/usr/bin/env python
"""
create a sequence of structures with different bond lengths between the two dimer atoms.
"""
from DFTB import XYZ
from DFTB.AtomicData import bohr_to_angs, atomic_number
from numpy import linspace

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 7:
        print "\nUsage: %s <atom 1> <atom 2> Rmin Rmax NR <dimer geometries>\n" % sys.argv[0]
        print "create dimer geometries with NR different bond lengths in the interval [Rmin,Rmax]."
        print "Rmin and Rmax are in Angstrom."
        exit(-1)
    ZA = int(atomic_number(sys.argv[1]))
    ZB = int(atomic_number(sys.argv[2]))
    Rmin, Rmax, NR = float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]) 
    dimer_file = sys.argv[6]
    def dimer(ZA, ZB, bond_length):
        atomA = (ZA, [0.0, 0.0, 0.0])
        atomB = (ZB, [0.0, 0.0, bond_length])
        return [atomA, atomB]
    bond_lengths = linspace(Rmin,Rmax,NR)
    XYZ.write_xyz(dimer_file, [dimer(ZA,ZB,bl/bohr_to_angs) for bl in bond_lengths], title="dimer trajectory")

