#!/usr/bin/env python
"""
check that the initial velocity distribution has the right temperature
"""
import glob
import sys
import numpy as np

from DFTB import XYZ, AtomicData

def plot_histogram(T):
    from matplotlib import pyplot as plt
    hist, bin_edges = np.histogram(T, normed=True)
    plt.xlabel("Temperature / K")
    plt.ylabel("P(T) dT")
    plt.plot((bin_edges[1:] + bin_edges[:-1])/2.0, hist)
    plt.show()
    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: python %s <pattern for initial conditions, e.g. structure_>" % sys.argv[0]
        sys.exit(-1)
    infiles = glob.glob(sys.argv[1] + "*.in")
    ekins = []
    for f in infiles:
        # compute kinetic energy for all initial conditions
        atomlist, vellist = XYZ.read_initial_conditions(f, units="bohr")
        masses = np.zeros(len(atomlist))
        for i,(Zi,posi) in enumerate(atomlist):
            masses[i] = AtomicData.atom_masses[AtomicData.atom_names[Zi-1]]
        vel = XYZ.atomlist2vector(vellist)
        ekin = sum(0.5 * masses * (pow(vel[0::3],2) + pow(vel[1::3],2) + pow(vel[2::3],2)))
        ekins.append(ekin)
    # convert kinetic energy to temperature
    # ekin = nvib * kB * T
    ekins = np.array(ekins)
    nat = len(atomlist)
    if nat < 2:
        raise Exception("Since vibrational and rotational degrees are removed, a single atom cannot have any meaningful temperature.")
    elif nat == 2:
        nvib = 1
    elif nat > 2:
        nvib = 3*nat-6
    print "number of atoms: %s" % nat
    print "number of vibrational modes: %s" % nvib
    T = ekins / (0.5 * nvib * AtomicData.kBoltzmann)
    print "average temperature: %.5f K" % np.mean(T)
    print "standard deviation:  %.5f K" % np.std(T)
    # make a histogram
    plot_histogram(T)
