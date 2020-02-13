#!/usr/bin/env python
"""
This script plots the trajectory populations
"""
import numpy as np

if __name__ == "__main__":
    pop = np.loadtxt("pop.dat")
    
    import matplotlib.pyplot as plt

    plt.rcParams["xtick.labelsize"] = 15
    plt.rcParams["ytick.labelsize"] = 15

    plt.xlabel("Time / fs", fontsize=15)
    plt.ylabel("Populations", fontsize=15)

    for i in range(0, Nst):
        plt.plot(pop[:,0], pop[:,1+i], lw=2)

    plt.legend()
    plt.show()

