#!/usr/bin/env python
"""
This script classifies the separate molecular fragments along a trajectory
and generates a table with the fractions of the molecules that have dissociated
via different reaction channels.

The program obminimize (from the openbabel project) is used to optimize snapshots along the trajectory 
so that the large fluctuation in bond length that can lead to wrong naming of a hot molecule is minimized.
"""
from DFTB import XYZ
from DFTB.Analyse import MolecularGraph

import numpy as np
import sys
import os
from os.path import dirname, join, isfile
import optparse



if __name__ == "__main__":
    usage  = "python %s <list of dynamics.xyz files>\n" % sys.argv[0]
    usage += " The trajectories are classified by reaction channels and the abundance of molecular fragments.\n"
    usage += " Type --help to see all options.\n"

    parser = optparse.OptionParser(usage)
    parser.add_option("--only_final", dest="only_final", default=1, type=int, help="This flag controls whether the whole trajectory should be analysed (0) or only the final geometry (1) [default: %default]")
    parser.add_option("--optimize", dest="optimize", default=0, type=int, help="In order to minimize bond length fluctuations in hot molecules that can lead to the wrong identification of fragments (e.g. CH3+H instead of methane), the geometries are optimized with a force field [default: %default]")
    parser.add_option("--overwrite", dest="overwrite", default=1, type=int, help="If a trajectory has already been optimized, the result is overwritten (1) or reused (0) [default: %default]")
    parser.add_option("--plot", dest="plot", default=0, type=int, help="Plot fractions of fragments and channels as a function of time [default: %default]")
    parser.add_option("--verbose", dest="verbose", default=0, type=int, help="verbose > 0 writes out messages from Babel [default: %default]")
    parser.add_option("--save_fragments", dest="save_fragments", default=1, type=int, help="If this flag is set (1) the fragments, that are discovered are saved to xyz-files with names fragment_#.xyz . 2D formulae are also drawn using 'babel' and are saved to png-files with names fragment_#.png . [default: %default]")
    
    (opts, args) = parser.parse_args()
    if len(args) < 1:
        print usage
        exit(-1)
    paths = args

    # A fragment is a part of a molecule that is separate from the rest by a large distance
    fragment_fractions = {}
    fragment_geometries = {}
    # A reaction channel is defined by a set of products that are created during the reaction
    channel_fractions = {}
    
    ntraj = np.array([])               # number of trajectories for each time step
    for i,dynamic_file in enumerate(paths):
        dyndir = dirname(dynamic_file)
        if opts.optimize == 1:
            # optimize each time step with a force field to remove bond length fluctuations in hot molecules
            opt_file = join(dyndir, "dynamics_opt_%.4d.xyz" % i)
            if (not isfile(opt_file)) or (opts.overwrite == 1):
                if opts.verbose == 0:
                    redirect = " 2> /dev/null"
                else:
                    redirect = ""
                print "Optimizing geometries in %s with UFF force field => %s" % (dynamic_file, opt_file)
                ret = os.system("obminimize -ff UFF %s %s | sed -e \"/^WARNING:/ d\" >  %s" % (dynamic_file, redirect, opt_file))
                assert ret == 0, "Optimization with obminimize failed!"
            else:
                print "found optimized %s" % opt_file
            dynamic_file = opt_file
        
        print "Identifying and classifying fragments in %s" % dynamic_file
        try:
            geometries = XYZ.read_xyz(dynamic_file)
        except IOError:
            print "WARNING: could not open %s. - skipped" % dynamic_file
            continue
        geometries = [ geometries[-1] ]
        fragtraj, fragdic = MolecularGraph.fragment_trajectory( geometries )
        fragment_geometries.update(fragdic)
        # length of the i-th trajectory
        nt = max(len(fragtraj), len(ntraj))
        # extend fraction list so that it has as many time steps as the longest trajectory
        if len(ntraj) < nt:
            nadd = nt-len(ntraj)
            ntraj = np.hstack((ntraj, np.zeros(nadd)))
        ntraj[:len(fragtraj)] += 1
        # fragments
        for k,fracs in fragment_fractions.iteritems():
            if len(fracs) < nt:
                nadd = nt-len(fracs)
                fragment_fractions[k] = np.hstack((fracs, np.zeros(nadd)))
        # channels
        for k,fracs in channel_fractions.iteritems():
            if len(fracs) < nt:
                nadd = nt-len(fracs)
                channel_fractions[k] = np.hstack((fracs, np.zeros(nadd)))

        # At each time step in a trajectory the fragments are identified and
        # the counts for the respective fragments and channels are increased, accordingly.
        for t,frags in enumerate(fragtraj):
            # fragments
            w = 1.0/len(frags)
            for f in frags:
                if fragment_fractions.has_key(f):
                    fragment_fractions[f][t] += w
                else:
                    fracs = np.zeros(nt)
                    fracs[t] = w
                    fragment_fractions[f] = fracs
            # channels
            frags.sort()
            channel = "+".join(tuple(frags))
            if channel_fractions.has_key(channel):
                channel_fractions[channel][t] += 1
            else:
                fracs = np.zeros(nt)
                fracs[t] = 1
                channel_fractions[channel] = fracs
                
        #
        if i == 0:
            # try to find the time step
            try:
                data = np.loadtxt(join(dyndir, "state.dat"))
                time = data[:,0]
                dt_fs = time[1]-time[0] # time step in fs
            except (IOError, IndexError) as e:
                print e
                dt_fs = 0.1 # assume a time step of 0.1 fs

    time = np.array(range(0,len(ntraj)))*dt_fs
    # FRAGMENTS
    # A table is saved with the fragment fractions as a function of time
    # and a table with the final fragment distribution at the last time step.
    data = [time]
    col_titles = "# TIME  "
    for k,f in fragment_fractions.iteritems():
        col_titles += "%s  " % k
        data += [f]
    data = np.vstack(data).transpose()
    # check
    for t in range(0, nt):
        assert abs(np.sum(data[t,1:]) - ntraj[t]) < 1.0e-5, "t = %s ntraj = %s   |w| = %s" % (t, ntraj[t], np.sum(data[t,1:]))
    # 
    fh = open("fragments.dat", "w")
    print>>fh, "# Fragment fractions "
    print>>fh, col_titles
    np.savetxt(fh, data)
    fh.close()
    print "Fragment fractions written to fragments.dat"

    print "          Final Fragment Distribution"
    print "          ==========================="
    final_fragment_distribution = []
    for k,v in fragment_fractions.iteritems():
        final_fragment_distribution.append( (k,v[-1]) )
    # sort distribution in descending order
    final_fragment_distribution.sort(key=lambda (k,v): -v)
    for (k,v) in final_fragment_distribution:
        print " %20.20s        %5.5f" % (k, v)

    # CHANNELS 
    data = [time]
    col_titles = "# TIME  "
    for k,f in channel_fractions.iteritems():
        col_titles += "%s  " % k
        data += [f]
    data = np.vstack(data).transpose()
    # check
    for t in range(0, nt):
        assert abs(np.sum(data[t,1:]) - ntraj[t]) < 1.0e-5, "t = %s ntraj = %s   |w| = %s" % (t, ntraj[t], np.sum(data[t,1:]))
    # 
    fh = open("channels.dat", "w")
    print>>fh, "# Channel fractions "
    print>>fh, col_titles
    np.savetxt(fh, data)
    fh.close()
    print "Channel fractions written to channels.dat"

    print "         Final Channel Distribution"
    print "         =========================="
    final_channel_distribution = []
    for k,v in channel_fractions.iteritems():
        final_channel_distribution.append( (k,v[-1]) )
    # sort distribution in descending order
    final_channel_distribution.sort(key=lambda (k,v): -v)
    for (k,v) in final_channel_distribution:
        print " %20.20s        %5.5f" % (k, v)

    if opts.save_fragments == 1:
        # generate png images of fragments
        png_files_str = ""  # a string holding the images sorted by abundance
        for (k,v) in final_fragment_distribution:
            frag_geom = fragment_geometries[k]
            frag_file = "fragment_%s.xyz" % k
            png_file  = "fragment_%s.png" % k
            XYZ.write_xyz(frag_file, [frag_geom], title="%s (%6.3f)" % (k, fragment_fractions[k][-1]))
            ret = os.system("obabel -ixyz %s -opng -O %s" % (frag_file, png_file))
            assert ret == 0
            print "saved png-file for fragment %s to %s" % (k, png_file)
            png_files_str += " %s" % png_file
        # use image magick to create a map between the labels and the formulae
        os.system("montage -mode concatenate -tile 4x4 %s map_fragments.png" % png_files_str)
        
    if opts.plot == 1:
        from matplotlib import pyplot as plt

        # plot fragment fractions
        plt.xlabel("time / fs")
        plt.ylabel("numbers of trajectories")
        vsum = None
        for k,v in fragment_fractions.iteritems():
            if v[-1] > 0.0:
                plt.plot(time, v, lw=2, label=k)
            if vsum == None:
                vsum = v
            else:
                vsum += v
        plt.plot(time, vsum, lw=2, ls="-.", label="Sum")
        plt.legend()
        plt.savefig("fragments.png")
        plt.show()

        # plot channel fractions
        vsum = None
        for k,v in channel_fractions.iteritems():
            if v[-1] > 0.0:
                plt.plot(time, v, lw=2, label=k)
            if vsum == None:
                vsum = v
            else:
                vsum += v
        plt.plot(time, vsum, lw=2, ls="-.", label="Sum")
        plt.legend()
        plt.savefig("channels.png")
        plt.show()
