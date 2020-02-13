#!/usr/bin/env python
"""
Analyse an MD trajectory in terms on bond indeces. 
The bond index of an atom is the sum over the distances to all other
atoms
"""
from DFTB import XYZ
from DFTB import AtomicData

from numpy import array
from numpy.linalg import norm
import sys

def bond_indeces(atomlist):
    """
    The bond index of an atom i is the sum of the
    distances to all other atoms
    """
    bindeces = []
    for i,(Zi,posi) in enumerate(atomlist):
        bi = 0.0 # bond index of atom i
        for j,(Zj,posj) in enumerate(atomlist):
            if i!=j:
                rij = norm(array(posi)-array(posj))
                bi += rij
        bindeces.append(bi)
    return bindeces

def plot_bond_indeces(atomlist, out_file):
    """
    show a plot of the calculated bond indeces 
    """
    from matplotlib.pyplot import plot, xlabel, ylabel, legend, show, savefig
    from numpy import loadtxt
    b = loadtxt(out_file)
    (ntraj,ndist) = b.shape
    xlabel("MD steps")
    ylabel("bond index $b_i = \\sum_j r_{ij}$") 
    for n in range(0, ndist):
        at = AtomicData.atom_names[atomlist[n][0]-1]
        plot(b[:,n], label="%s%d" % (at,n))
    legend()
    savefig(out_file + ".png")
    show()
                
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: %s <xyz-file with list of nuclear geometries> <output file with bond indeces>" % sys.argv[0]
        exit(-1)
    xyz_file = sys.argv[1]
    out_file = sys.argv[2]
    fh = open(out_file, 'w')
    for atomlist in XYZ.read_xyz_it(xyz_file):
        for bi in bond_indeces(atomlist):
            print>>fh, "%.7f " % bi, 
        print>>fh, ""
    fh.close()
    plot_bond_indeces(atomlist, out_file)
