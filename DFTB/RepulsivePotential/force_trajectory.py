"""
compute gradients along a trajectory of geometries using Gaussian
"""
from DFTB import XYZ
import Gaussian

if __name__ == "__main__":
    import sys
    import os
    if len(sys.argv) < 3:
        print "Usage: %s <geometry file> <gradient file>" % sys.argv[0]
        print "Compute gradients for all geometries in <geometry file> and write them to <gradient file>"
        exit(-1)
    geom_file = sys.argv[1]
    grad_file = sys.argv[2]
    tmp_dir = "/scratch/humeniuka/dftb/"
    os.system("module load g09")
    for atomlist in XYZ.read_xyz(geom_file):
        tmp_com = tmp_dir + "gaussian.com"
        tmp_out = tmp_dir + "gaussian.log"
        Gaussian.write_input(tmp_com, atomlist, \
                        route="# PBE/6-311G(3df,3pd)++ Force", \
                        title="Gradients for fitting repulsive potential")
        # compute forces with gaussian
        os.system("g09 < %s > %s" % (tmp_com, tmp_out))
        try:
            forces = Gaussian.read_forces(tmp_com)
        except Gaussian.FormatError as e:
            print e
        XYZ.write(grad_file, [forces], title="forces", units="forces", mode='a')
