# -*- coding: utf-8 -*-
"""
uncontracted cartesian Gaussian-type orbitals
"""

import numpy as np
import numpy.linalg as la

from DFTB.MolecularIntegrals.integrals import norm

def cartesian_shell_ordering(LType):
    """
    ordering of cartesian basis function in one shell
    according to the Gaussian 09 source code in file 
    'utilnz.F'  lines  70053-70064

    Parameters
    ----------
    LType  :  L value of shell, 0: s, 1: p, 2: 6D, etc.
    
    Returns
    -------
    iterator to labels of cartesian basis functions
    """
    # s,p,6d and 10f shells
    ordering = [
            ["S"],                              # S
            ["X","Y","Z"],                      # P
            ["XX","YY","ZZ", "XY","XZ","YZ"],   # 6D
            ["XXX","YYY","ZZZ","XYY","XXY","XXZ","XZZ","YZZ","YYZ","XYZ"]]  # 10F
    if LType < 4:
        for XLab in ordering[LType]:
            yield XLab
        raise StopIteration
    # higher shells: g, h, i, ...
    XLab = ['' for i in range(0, 100)]
    for L in range(0, LType+1):
        for M in range(0, LType-L+1):
            for I in range(1, L+1):
                XLab[I] = 'X'
            for I in range(1, M+1):
                XLab[L+I] = 'Y'
            for I in range(1, LType-L-M+1):
                XLab[L+M+I] = 'Z'
            yield "".join(XLab[1:LType+1])

class UncontractedBasisSet:
    def __init__(self, nbfs, exponents, shell_types, shell_coords, shell_atoms):
        """
        initialize basis for `nbfs` primitive Gaussian functions

        Parameters
        ----------
        nbfs        : int, number of basis functions
        exponents   : 1d numpy array of floats, list of exponents of Gaussians
        shell_types : 1d numpy array of integers, 
                      angular momentum of each shell (0: s, 1: p, 2: d, ...)
        shell_coords: centers of basis functions, 
                      shell_coords[:,i] are the cartesian coordinates for the i-th shell
        shell_atoms : map between shell index and atom index
        """
        assert np.all(shell_types >= 0), "Only cartesian basis functions (6D, 10F, ...) are supported!"
        
        self.nbfs = nbfs
        self.exponents = np.zeros(nbfs,     dtype=float)
        self.powers    = np.zeros((3,nbfs), dtype=int)
        self.centers   = np.zeros((3,nbfs), dtype=float)

        # index of atomic center to which the basis function belongs
        self.atom_ids  = np.zeros(nbfs,     dtype=int)
        # index of angular momentum shell
        self.shell_index = np.zeros(nbfs,   dtype=int)
        
        ibf = 0   # `ibf` counts cartesian basis functions
        # `ish` counts the shells
        for ish, (alpha, ltype) in enumerate(zip(exponents, shell_types)):
            # enumerate cartesian basis functions in this shell in the same order
            # used by Gaussian 09
            for powstr in cartesian_shell_ordering(ltype):
                self.exponents[ibf] = alpha
                self.powers[:,ibf]  = powstr.count('X'), powstr.count('Y'), powstr.count('Z')
                self.centers[:,ibf] = shell_coords[:,ish]
                #
                self.atom_ids[ibf] = shell_atoms[ish]
                self.shell_index[ibf] = ish
                ibf += 1
        assert ibf == nbfs
                
    def wavefunction(self, orb, x,y,z):
        """
        Parameters
        ----------
        orb  : MO coefficients, orb[i] is the coefficient of the i-th basis function
        x,y,z: numpy grids with x-,y- and z-coordinates at which the orbital should
               be evaluated

        Returns
        -------
        wfn:   amplitude of orbital at the grid points
        """
        wfn = 0j*x
        assert len(orb) == self.nbfs
        for i in range(0, self.nbfs):
            wfn += orb[i] * wavefunction(self.exponents[i], self.powers[:,i], self.centers[:,i], x,y,z)
        return wfn

def wavefunction(alpha, (l,m,n), A, x,y,z):
    """evaluate a cartesian Gaussian basis function on a grid"""
    nrm = norm(alpha,(l,m,n))
    r2 = (x-A[0])**2 + (y-A[1])**2 + (z-A[2])**2
    wfn = nrm*(x-A[0])**l * (y-A[1])**m * (z-A[2])**n * np.exp(-alpha*r2)
    return wfn
    
if __name__ == "__main__":
    # print ordering of basis functions in each shell
    for LType in [0,1,2,3,4,5,6,7]:
        print "Shell Type: %s" % LType
        print list(cartesian_shell_ordering(LType))

