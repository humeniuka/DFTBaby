#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
read data from formatted checkpoint file (.fchk) produced by Gaussian 09
"""

from DFTB.MolecularIntegrals.gtobasis import UncontractedBasisSet
from DFTB.MolecularIntegrals.integrals import basis_overlap

import numpy as np
import numpy.linalg as la

def read_fchk(filename):
    """
    read data from formatted checkpoint file

    Parameters
    ----------
    filename   :  path to .fchk file

    Returns
    -------
    data       :  dictionary with field names as keys and numerical data as values

    References
    ----------
    The format of fchk-file is described in 
    http://wild.life.nctu.edu.tw/~jsyu/compchem/g09/g09ur/f_formchk.htm
    """
    fh = open(filename)
    # first two lines contain description of the job
    title = fh.readline()
    calc_type, method, basis = fh.readline().split()
    # dictionary is filled with data
    data = {}
    while True:
        line = fh.readline()
        if line == "":
            # end of file reached
            break
        # name of data block
        name = line[:43].strip()
        # size and type of data ( I = Integer, R = Real)
        words = line[43:].split()
        # type of data 
        if (words[0] == "I"):
            dtype = int
            # integer arrays are written in chunks of 6 integers per line
            chunk_size = 6
        elif (words[0] == "R"):
            dtype = float
            # float arrays are written in chunks of 5 floats per line
            chunk_size = 5
        elif (words[0] == "L"):
            dtype = bool
            chunk_size = 6
        else:
            dtype = str
            chunk_size = 1
              
        # size of data
        if len(words) == 3:
            size = int(words[2])
        else:
            size = 1

        # read numerical data
        if size == 1 and len(words) == 2:
            # data of length 1 is written on the same line as the name of the section
            data[name] = dtype(words[1])
            continue
        else:
            arr = []
            # number of lines to read
            nl = int(np.ceil(size/float(chunk_size)))
            # 
            for i in range(0, nl):
                # append chunks of numerical data read in this line
                arr += map(dtype, fh.readline().split())
            # convert list to numpy array
            data[name] = np.array(arr, dtype=dtype)
            assert len(data[name]) == size
            
    return data


class G09ResultsDFT:
    """
    stores results from a DFT calculation with Gaussian 09
    """
    def __init__(self, fchk_file):
        """read results from formatted checkpoint file"""
        print "reading data from formatted checkpoint file '%s'" % fchk_file
        fchk = read_fchk(fchk_file)
        assert fchk['Largest degree of contraction'] == 1, "Basis set should be fully uncontracted!"

        self.nat = fchk['Number of atoms']
        # 
        self.coordinates = np.reshape(fchk['Current cartesian coordinates'],
                                      (self.nat, 3)).transpose()

        self.atomic_numbers = fchk['Atomic numbers']
        
        self.nelec_alpha = fchk['Number of alpha electrons']
        self.nelec_beta  = fchk['Number of beta electrons']
        self.nelec = self.nelec_alpha + self.nelec_beta
        
        self.nbfs = fchk['Number of basis functions']
        self.nmo = fchk['Number of independent functions']

        self.orbe_alpha = fchk['Alpha Orbital Energies']
        if 'Beta Orbital Energies' in fchk:
            self.orbe_beta  = fchk['Beta Orbital Energies']
        else:
            # closed shell - energies for alpha and beta orbitals are the same
            self.orbe_beta  = fchk['Alpha Orbital Energies']

        # In the fchk-file the M.O. coefficients are stored in column-major order
        self.orbs_alpha = np.reshape(fchk['Alpha MO coefficients'],
                                     (self.nmo, self.nbfs)).transpose()
        # Now orbs_alpha[:,i] contains the coefficients for the i-th molecular orbital
        if 'Beta MO coefficients' in fchk:
            self.orbs_beta  = np.reshape(fchk['Beta MO coefficients'],
                                         (self.nmo, self.nbfs)).transpose()
        else:
            # closed shell - copy coefficients
            self.orbs_beta  = np.reshape(fchk['Alpha MO coefficients'],
                                         (self.nmo, self.nbfs)).transpose()

        # initialize the uncontracted basis set
        exponents = fchk['Primitive exponents']
        shell_types = fchk['Shell types']        
        shell_coords = np.reshape(fchk['Coordinates of each shell'],
                                  (len(exponents), 3)).transpose()
        shell_atoms = fchk['Shell to atom map']-1
        
        self.basis = UncontractedBasisSet(self.nbfs,
                                          exponents, shell_types, shell_coords, shell_atoms)
        
    def consistency_checks(self):
        """
        perform some consistency checks on the loaded data
        """
        print "consistency checks..."
        # check that the M.O.s are orthogonal
        # overlap between A.O.s
        Sao = basis_overlap(self.basis)
        Smo_alpha = np.dot(self.orbs_alpha.transpose(), np.dot(Sao, self.orbs_alpha))
        Smo_beta  = np.dot(self.orbs_beta.transpose(),  np.dot(Sao, self.orbs_beta))
        # accuracy for passing checks
        eps = 1.0e-4
        err = la.norm(Smo_alpha - np.eye(self.nmo))
        assert err < eps, "alpha M.O.s not normalized, |Ca^T.S.Ca - Id|= %e" % err
        err = la.norm(Smo_beta - np.eye(self.nmo))
        assert err < eps, "beta M.O.s not normalized, |Cb^T.S.Cb - Id|= %e" % err

        
if __name__ == "__main__":
    import sys
    import os.path
    
    if len(sys.argv) < 2:
        print "Usage: %s  <formatted checkpoint file>" % os.path.basename(sys.argv[0])
        print "  read data from formatted checkpoint file"
        exit(-1)
        
    filename = sys.argv[1]
    fchk = read_fchk(filename)

    for key, value in fchk.iteritems():
        print key
        print value

    res = G09ResultsDFT(filename)
    res.consistency_checks()
