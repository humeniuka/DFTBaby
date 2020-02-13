#!/usr/bin/env python
import sys
import os.path

usage="""

perturb equilibrium geometry in different ways to create paths for 
fitting repulsive potentials. 

    %s <equilibrium geom. .xyz> <fit path .xyz> <command ...>

The first argument should point to a xyz-file with the equilibrium geometry. The geometries
of the fit paths are saved to the filename in the second argument. The remaining arguments
specify the type of perturmations.

The following perturbation commands are available:

      'scale' smin smax N

            scale all position vectors by a factor s=smin(1),...,smax(N)

      'dislocate' at r nshells N

            displace atom at on `nshells` equidistant shells in a sphere of radius `r` (Angstrom)
            around its equilibrium position. For each shell N random displacements are generated.

Example:

     %s  h2.xyz h2_stretched.xyz    scale  0.5 2.0 10

creates a fit path `h2_stretched.xyz` with 10 geometries where the bond length is changed in equidistant
steps from 0.5*rHH to 2.0*rHH. rHH is the equilibrium bond length in `h2.xyz`.

""" % (os.path.basename(sys.argv[0]), os.path.basename(sys.argv[0]))

from DFTB import XYZ, AtomicData

import numpy as np
import numpy.linalg as la

import copy

class PerturbedGeometries:
    def __init__(self, atomlist0):
        """
        Parameters
        ----------
        atomlist: equilibrium geometry
        """
        self.atomlist0 = atomlist0

    def scale(self, smin=0.6, smax=2.0, N=20):
        """
        scale all position vectors by a factor smin < s < smax
        
        Parameters:
        ===========
        smin: scale factor for smallest structure
        smax: scale factor for largest structure
        N: number of points in the interval [smin, smax]
        """
        scaled_geoms = []
        for s in np.linspace(smin, smax, N):
            scaled_atomlist = []
            for (Zi,[xi,yi,zi]) in self.atomlist0:
                scaled_atomlist.append( (Zi, (xi*s,yi*s,zi*s)) )
            scaled_geoms.append(scaled_atomlist)
        return scaled_geoms

    def dislocate(self, atom=1, radius=0.5, nshells=4, N=5):
        """

        """
        # convert radius to bohr
        radius /= AtomicData.bohr_to_angs
        # selected atom
        Z_sel,pos_sel = self.atomlist0[atom-1]
        pos_sel = np.array(pos_sel)

        # equilibrium geometry 
        disloc_geoms = [self.atomlist0]

        for s in range(1,nshells+1):
            rs = s*radius/float(nshells)
            for i in range(0, N):
                # random unit vector
                uvec = np.random.rand(3)
                uvec /= la.norm(uvec)
                # position of selected atom is displaced along random vector
                # on shell s
                pos_disloc = pos_sel + rs * uvec
                #
                disloc_atomlist = copy.deepcopy(atomlist0)
                disloc_atomlist[atom-1] = (Z_sel, pos_disloc)

                disloc_geoms.append( disloc_atomlist )

        return disloc_geoms

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print usage
        exit(-1)

    args = sys.argv[1:]

    if len(args) < 2:
        print "Missing arguments:"
        print "    xyz_input xyz_output"
        exit(-1)
        
    xyz_in = args[0]
    xyz_out = args[1]

    if len(args) < 3:
        print "possible commands: 'scale', 'dislocate'"
        exit(-1)
    cmd = args[2]
    # the remaining command line arguments depend
    cmd_args = map(eval, args[3:])

    atomlist0 = XYZ.read_xyz(xyz_in)[-1]
    PG = PerturbedGeometries(atomlist0)

    kwds = XYZ.extract_keywords_xyz(xyz_in)
    
    if cmd == "scale":
        if len(cmd_args) < 3:
            print "Arguments for 'scale':"
            print "    smin smax N"
            exit(-1)
        geometries = PG.scale(*cmd_args)
    elif cmd == "dislocate":
        if len(cmd_args) < 4:
            print "Arguments for 'dislocate':"
            print "    atom radius nshells N"
            exit(-1)
        geometries = PG.dislocate(*cmd_args)
    else:
        raise ValueError("Command '%s' not understood" % cmd)

    XYZ.write_xyz(xyz_out, geometries, title="charge=%d" % kwds.get("charge",0))
    print "fit path written to '%s'" % xyz_out

    
