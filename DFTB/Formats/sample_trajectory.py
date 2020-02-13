#!/usr/bin/env python
"""
sample initial conditions from a trajectory in the dynamics.out
"""

from DFTB.Formats import DynamicsOut

if __name__ == "__main__":
    import sys
    from os.path import join
    from optparse import OptionParser
    usage = "Usage: python %s <path to dynamics.out> <output directory>\n" % sys.argv[0]
    usage += "  sample initial conditions (geometries + velocities) from\n"
    usage += "  the ground state trajectory given in <dynamics.out> and\n"
    usage += "  save the initial conditions to dynamics_####.in in the output directory\n"
    usage += "  Type --help to see all options.\n"

    parser = OptionParser(usage)
    parser.add_option("--skip", type=int, dest="skip", help="Skip the first time steps. [default: %default]",  default=200)
    parser.add_option("--freq", type=int, dest="freq", help="Extract geometries and velocities from every FREQ-th time step [default: %default]", default=100)
    (opts, args) = parser.parse_args()

    if len(args) < 2:
        print usage
        exit(-1)

    dynout = args[0]
    inidir = args[1]

    i = 0 # count initial conditions
    for it,(t,atomlist,velocities) in enumerate(DynamicsOut.t_geometries_velocities(dynout)):
        if it > opts.skip and it % opts.freq == 0:
                ini_file = join(inidir, "dynamics_%.4d.in" % i)
                print "time step %s fs   =>  %s" % (t, ini_file)
                DynamicsOut.write_initial_conditions(atomlist, velocities, ini_file)
                i += 1
