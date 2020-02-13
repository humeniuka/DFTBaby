#!/usr/bin/env python
"""
plot adiabatic energies along a trajectory and highlight the current state
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":
    usage = "python %s <state.dat file> <... list of energy_##.dat files>\n"
    usage += " show adiabatic energies along a trajectory and highlight current state\n"

    if len(sys.argv) < 3:
        print usage
        exit(-1)

    st = np.loadtxt(sys.argv[1])
    times = st[:,0]
    states = st[:,1].astype(int)
    energies = []
    for f in sys.argv[2:]:
        if f[-3:] == "dat":
            en = np.loadtxt(f)
            assert np.sum(abs(en[:,0] - times)) < 1.0e-10
            energies.append( en[:,1] )
    energies = np.array(energies) * 27.211

    Nst, Nt = energies.shape
    print "%s electronic states" % Nst
    print "%s time steps" % Nt

    # energy of active state
    en_active = np.zeros(Nt)
    for t in range(0, Nt):
        en_active[t] = energies[states[t],t]

    plt.rcParams["xtick.labelsize"] = 15
    plt.rcParams["ytick.labelsize"] = 15

    plt.xlabel("Time / fs", fontsize=15)
    plt.ylabel("Energy / eV", fontsize=15)

    for i in range(0, Nst):
        plt.plot(times, energies[i,:], lw=2)
    plt.plot(times, en_active, "o", lw=4, color="red", alpha=0.5, label="active state")

    plt.legend()
    plt.show()
    
