#!/usr/bin/env python

import sys
import os.path
import optparse

usage="""
Usage: %s <list of fragment specifications and xyz-files>

sum additional numerical data in columns 5,6,... for atoms in the
same fragment and average the result over all xyz-files.

Fragments are specified as quoted strings containing the indeces 
of atoms as ranges (e.g. 10-20) and/or comma-separated lists:

  "1-10"              fragment contains atoms 1 to 10
  "1,2,3,10-20"       fragments contains atoms 1,2,3 and 10-20

Atom indeces start at 1. File names and coordinate specifications 
can be interspersed on the command line.
""" % os.path.basename(sys.argv[0])

from DFTB.QMMM import parseAtomTags

import numpy as np

def read_atomic_data_it(filename):
    """
    read additional atomic data from columns 5,6,...,ncol of 
    an xyz-file.

    The iterator returns one block of data with shape (nat,ncol)
    at a time which can be processed in a pipeline. The current
    time (in fs) is extracted from the comment line.

    Parameters:
    ===========
    filename: path to xyz-file

    Returns:
    ========
    iterator to tuples (time, data) 
    """
    fh = open(filename)
    igeo = 1
    while 1:
        line = fh.readline()
        if not line:
            # end of file reached
            break
        words = line.strip().split()
        try:
            nat = int(words[0])
        except ValueError as e:
            print e
            raise Exception("Probably wrong number of atoms in xyz-file '%s'" % filename)
        # extract time of current step from title line
        title = fh.readline()
        words = title.strip().split()
        assert words[2] == "fs"
        time = float(words[1])
        # read numerical data from columns 5,6,..., for `nat` atoms
        data = []
        for i in xrange(nat):
            line = fh.readline()
            words = line.split()
            cols = map(float,words[4:])
            data.append(cols)
        igeo += 1
        data = np.array(data)
        yield time,data
        
if __name__ == "__main__":
    parser = optparse.OptionParser(usage)
    
    (opts,args) = parser.parse_args()
    
    if len(args) == 0:
        print usage
        exit(-1)

    # extract list of xyz-files from the command line
    xyz_files = []
    # list of list with atomic indices belonging to each fragment
    fragments = []
    fragment_labels = []
    for arg in args:
        try:
            # indices of atoms "1,2,3" or "9-14,15,16" or "1-10,15-20"
            atom_ids = parseAtomTags(arg)
            # internally indices are 0-based
            atom_ids = np.array(atom_ids) -1
            fragments.append(atom_ids)
            fragment_labels.append(arg)
        except ValueError:
            # probably this argument is an xyz-file
            assert ".xyz" in arg, "%s does not have .xyz suffix nor is it a list of integer" % arg
            xyz_files.append(arg)

    if len(xyz_files) < 1:
        print "At least one xyz-file should be specified!"
        exit(-1)
    data_trajs = []
    for f in xyz_files:
        data_timeseries = []
        times = []
        for time_fs,data_atoms in read_atomic_data_it(f):
            times.append(time_fs)
            data_frags = []
            for atom_ids in fragments:
                # sum data over all atoms in one fragment
                data_frags.append( np.sum(data_atoms[atom_ids,:], axis=0) )
            data_timeseries.append(data_frags)
        data_trajs.append( np.array(data_timeseries) )

    ncols = data_trajs[0].shape[-1]
    nfrags = len(fragments)    
    ntrajs = len(xyz_files)
    nsteps = min(map(len, data_trajs))     # find length of shorted trajectories
    print "additional columns in xyz-file : %d" % ncols
    print "number of fragments            : %d" % nfrags
    print "number of trajectories         : %d" % ntrajs
    print "time steps                     : %d" % nsteps
    #
    data = np.zeros((ntrajs,nsteps,nfrags,ncols), dtype=float)
    for i in range(0, ntrajs):
        data[i,:,:,:] = data_trajs[i][:nsteps,:,:]
    # average over trajectories
    data_avg = np.sum(data, axis=0) / float(ntrajs)
    
    # write file with averages for each fragment
    for k,fragment_label in enumerate(fragment_labels):
        fname = "fragment_%s.dat" % fragment_label
        atom_ids = fragments[k]
        with open(fname, "w") as fh:
            print>>fh, "# Atom indices belonging to this fragment:"
            coord_str = "-".join(map(str,np.array(atom_ids)+1))
            print>>fh, "# %s" % coord_str
            print>>fh, "# Time / fs              Data"
            data_avg_frag = np.zeros((nsteps,ncols+1))
            # insert time axis
            data_avg_frag[:,0] = times[:nsteps]
            data_avg_frag[:,1:] = data_avg[:,k,:]
            np.savetxt(fh, data_avg_frag)
            print "trajectory average for fragment %s written to file '%s'" % (fragment_label, fname)
