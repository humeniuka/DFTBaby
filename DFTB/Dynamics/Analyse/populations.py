#!/usr/bin/env python
"""
This script reads the 'state.dat' files of the individual trajectories
and computes the average state populations.
"""

import numpy as np
import sys
import os

if __name__ == "__main__":
    usage  = "%s <list of state.dat files>\n" % os.path.basename(sys.argv[0])
    usage += " The populations are saved to the file 'pop.dat'"

    if len(sys.argv) < 2:
        print usage
        exit(-1)
    paths = sys.argv[1:]

    tmin=0.0
    tmax=0.0
    Nst=0
    Nt=0
    data_list = []
    for i,state_file in enumerate(paths):
        data = np.loadtxt(state_file)
        Nst = max(data[:,1].max()+1, Nst)
        tmin = min(data[:,0].min(), tmin)
        tmax = max(data[:,0].max(), tmax)
        Nt = max(len(data[:,0]), Nt)
        data_list.append(data) 
    Nst = int(Nst)
    Nt = int(Nt)
    print "%d electronic states" % Nst
    print "%d time steps between %s and %s" % (Nt, tmin, tmax)
    pop = np.zeros((Nt,Nst+1))
    pop[:,0] = np.linspace(tmin,tmax,Nt) # time axis in fs
    Ntraj = [0 for t in range(0, Nt)]  # count trajectories available at each time step
    for i,data in enumerate(data_list):
        # only consider trajectories that finished nicely
        Nt_i = data.shape[0]
        if Nt_i != Nt:
            print "Trajectory %d has only %d time steps" % (i,Nt_i)
        for t in range(0, Nt_i):
            st = int(data[t,1])
            pop[t,st+1] += 1
            Ntraj[t] += 1.0
    print "%s trajectories" % Ntraj[0]
    # divide populations by the number of trajectories
    for t in range(0, Nt):
        pop[t,1:] /= float(Ntraj[t])
    fh = open("pop.dat", "w")
    print>>fh, "# TIME / fs               ",
    for st in range(0, Nst):
        print>>fh, "POP_%d                   " % st,
    print>>fh,""
    np.savetxt(fh, pop)
    fh.close()
    print "Populations saved to file 'pop.dat'"

