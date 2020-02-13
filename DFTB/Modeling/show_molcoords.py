#!/usr/bin/env python

import sys
import os.path
import optparse

usage="""
Usage: %s <list of coordinate specifications and xyz-files>

shows bond lengths, angles and dihedrals for selected combinations
of atoms in an xyz-file. 

Coordinates are specified as quoted strings containing the indeces
of
  2 atoms A and B     for the distance A-B (e.g. "1 2")
  3 atoms A,B and C   for the angle B-A-C (e.g. "4 10 3")
  4 atoms A,B,C and D for the dihedral angle between the planes
                      spanned by A->B,B->C and B->C,C->D.

Atom indeces start at 1. Distances are in Angstrom and angles in degrees.

The selected internal coordinates are printed for each geometry in the 
file. If there are several xyz-files the coordinates are averaged over all of them. 
File names and coordinate specifications can be interspersed on the command line.

Examples:
   "1 2" water.xyz "1 2 3"    
          -> shows distance between atoms 1 and 2 and the angle 1-2-3.
   traj1.xyz traj2.xyz traj3.xyz "3 4 5"
          -> shows the bond angle 3-4-5 for each time step averaged over all 3 trajectories.
""" % os.path.basename(sys.argv[0])

from DFTB.Modeling import MolecularCoords as MolCo
from DFTB import XYZ, AtomicData

import numpy as np

def get_coord(atomlist, atom_ids):
    n = len(atom_ids)
    if n == 2:
        # bond length
        A,B = atom_ids
        coord = MolCo.distance(atomlist[A], atomlist[B])
        coord *= AtomicData.bohr_to_angs
    elif n == 3:
        # angle
        A,B,C = atom_ids
        coord = MolCo.angle(atomlist[A], atomlist[B], atomlist[C])
        coord *= 180.0/np.pi
    elif n == 4:
        # dihedral
        A,B,C,D = atom_ids
        coord = MolCo.dihedral_angle(atomlist[A], atomlist[B], atomlist[C], atomlist[D])
        coord *= 180.0/np.pi
    else:
        raise ValueError("Internal coordinate should be specified by 2, 3 or 4 atom indeces!")
    return coord
        
if __name__ == "__main__":
    parser = optparse.OptionParser(usage)
    parser.add_option("-n", "--no_header", dest="print_header", default=True,
                      action="store_false", help="suppress header with atom labels")
    
    (opts,args) = parser.parse_args()
    
    if len(args) == 0:
        print usage
        exit(-1)

    # extract list of xyz-files from the command line
    xyz_files = []
    coord_specs = []
    for arg in args:
        try:
            atom_ids = map(lambda i: int(i)-1, arg.split())
            coord_specs.append(atom_ids)
        except ValueError:
            # probably this argument is an xyz-file
            assert ".xyz" in arg, "%s does not have .xyz suffix nor is it a list of integer" % arg
            xyz_files.append(arg)

    if len(xyz_files) < 1:
        print "At least one xyz-file should be specified!"
        exit(-1)
    coords_trajs = []
    for f in xyz_files:
        coords_timeseries = []
        for geom in XYZ.read_xyz_it(f):
            coords = []
            for atom_ids in coord_specs:
                coords.append( get_coord(geom, atom_ids) )
            coords_timeseries.append(coords)
        coords_trajs.append( np.array(coords_timeseries) )
        
    ntrajs = len(xyz_files)
    nsteps = min(map(len, coords_trajs))     # find length of shorted trajectories
    ncoords = len(coord_specs)
    #
    coords_data = np.zeros((ntrajs,nsteps,ncoords), dtype=float)
    for i in range(0, ntrajs):
        coords_data[i,:,:] = coords_trajs[i][:nsteps,:]
    # average over trajectories
    coords_avg = np.sum(coords_data, axis=0) / float(ntrajs)
    # write header unless suppressed
    if opts.print_header == True:
        print "#",
        for atom_ids in coord_specs:
            coord_str = "-".join(map(str,np.array(atom_ids)+1))
            print ("   %s  " % coord_str).ljust(24),
        print ""
    np.savetxt(sys.stdout, coords_avg)
    
