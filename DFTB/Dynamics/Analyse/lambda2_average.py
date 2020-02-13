#!/usr/bin/env python
"""
average current Lambda2 over all trajectories
"""

import numpy as np

if __name__ == "__main__":
    import sys
    import os.path
    
    usage = """
  Usage:  %s  <list of lambda2_current.dat files>

    compute the trajectory average of lambda2 and save it to the file 'lambda2_average.dat'
    """ % os.path.basename(sys.argv[0])

    if len(sys.argv) < 2:
        print usage
        exit(-1)
        
    paths = sys.argv[1:]

    tmin=0.0    
    tmax=0.0    
    Nt=0         # number of time steps

    data_list = []
    for i,lambda2_file in enumerate(paths):
        data = np.loadtxt(lambda2_file)
        tmin = min(data[:,0].min(), tmin)
        tmax = max(data[:,0].max(), tmax)
        Nt = max(len(data[:,0]), Nt)
        data_list.append(data) 

    # averaged Lambda2
    lambda2_avg = np.zeros((Nt,2))
    lambda2_avg[:,0] = np.linspace(tmin,tmax,Nt) # time axis in fs
    Ntraj = [0 for t in range(0, Nt)]  # count trajectories available at each time step

    for i,data in enumerate(data_list):
        # warn about trajectories that did not finish nicely
        Nt_i = data.shape[0]
        if Nt_i != Nt:
            print "Trajectory %d has only %d time steps" % (i,Nt_i)
        for t in range(0, Nt_i):
            if data[t,1] >= 0.0:
                #lambda2 < 0 indicates invalid data (lambda2 is undefined on S0)
                lambda2_avg[t,1] += data[t,1]
                Ntraj[t] += 1.0
    print "%s trajectories" % Ntraj[0]
    # divide lambda2_avg by the number of trajectories
    for t in range(0, Nt):
        lambda2_avg[t,1] /= float(Ntraj[t])
    fh = open("lambda2_average.dat", "w")
    print>>fh, "# TIME / fs               LAMBDA2"
    np.savetxt(fh, lambda2_avg)
    fh.close()
    print "Averaged Lambda2 values saved to file 'lambda2_average.dat'"

