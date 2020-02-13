#!/usr/bin/env python
from DFTB.MolecularIntegrals import settings
from DFTB.MolecularIntegrals.BasissetFreeDFT import BasissetFreeDFT
import numpy as np

def hmi_continuum(l, m, E):
    """
    compute continuum orbitals of the hydrogen molecular ion H2+

    Parameters
    ----------
    l,m         :  angular quantum numbers of asymptotic solution
                   e.g.  l=0,m=0  s-orbital
                         l=1,m=+1 px-orbital
    E           :  energy (in a.u.) of continuum orbital, E = 1/2 k^2
    """
    # H2^+
    # bond length in bohr
    R = 2.0
    atomlist = [(1, (0.0, 0.0, -R/2.0)),
                (1, (0.0, 0.0, +R/2.0))]
    
    # choose resolution of multicenter grids for continuum orbitals
    settings.radial_grid_factor = 120      # controls size of radial grid  
    settings.lebedev_order = 25          # controls size of angular grid

    RDFT = BasissetFreeDFT(atomlist, None, charge=+1)

    # This is a one-electron system, so there are no other occupied orbitals
    def rho(x,y,z):
        return 0*x
    def homo(x,y,z):
        return 0*x
    
    delta, phi = RDFT.solveScatteringProblem(rho, homo, E, l, m)

    
if __name__ == "__main__":
    import sys
    import os.path

    args = sys.argv[1:]
    if len(args) < 3:
        usage = """
        Usage:
                %s  l   m  E

           compute the continuum orbital of H2^+ (hydrogen molecular ion)

        Parameters:
             l,m       -  integers, -l <= m <= l, angular quantum numbers
                          of asymptotic solution
             E         -  float, energy of continuum orbital is E = 1/2 k^2
        
        """ %  os.path.basename(sys.argv[0])

        print usage
        exit(-1)

    l = int(args[0])
    m = int(args[1])
    E = float(args[2])

    hmi_continuum(l, m, E)

    
