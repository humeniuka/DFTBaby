#!/usr/bin/env python

import sys
import numpy as np

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: python %s <list of displacement.dat files>" % sys.argv[0]
        print "  averages dislacements over all trajectories and writes"
        print "  the result to the file 'displacements_avg.dat'."
        exit(-1)
    displacements = []
    for dat_file in sys.argv[1:]:
        data = np.loadtxt(dat_file)
        displacements.append( data )
    # find shortest trajectory
    nsteps = min([data.shape[0] for data in displacements])
    # average over all trajectories
    ntraj = len(displacements)
    displacements_avg = np.zeros((nsteps,4))
    for i in range(0, ntraj):
        displacements_avg += displacements[i][:nsteps]
    displacements_avg /= float(ntraj)

    fh = open("displacements_avg.dat", "w")
    print>>fh, "# TSTEP     R_X / Angstrom    R_Y / Angstrom   R_Z / Angstrom"
    np.savetxt(fh, displacements_avg)
    fh.close()
    
