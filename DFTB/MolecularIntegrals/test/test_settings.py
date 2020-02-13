#!/usr/bin/env python
"""
This script demonstrates how to change the resolution of the multicenter 
grid.
"""

# routine for integration
from DFTB.MolecularIntegrals.Ints1e import kinetic
from DFTB.BasisSets import AtomicBasisSet

if __name__ == "__main__":
    # default settings for grid size
    from DFTB.MolecularIntegrals import settings
    # increase radial grid
    settings.radial_grid_factor = 10
    settings.lebedev_order = 23
    
    atomlist = [(6, (0,0,0))]
    basis = AtomicBasisSet(atomlist)

    # matrix elements of kinetic energy computed on a fine radial grid
    for a,bfA in enumerate(basis.bfs):
        for b,bfB in enumerate(basis.bfs):
            t = kinetic(atomlist, bfA, bfB)
            print "(%d|T|%d) = %e" % (a,b, t)
    
