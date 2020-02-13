#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The Mulliken approximation for electron repulsion integrals fails if different
orbitals are located on the same center, since in this case the overlap
matrix is diagonal. Therefore these integrals have to be calculated exactly
by numerical integration. This has to be done only once, since the integrals
do not depend on the molecular geometry.

The pseudo-orbitals a and b are both located on the same atom, 
whose neutral density is rho0. The onsite integrals are
the Coulomb-like integral

   (aa|(1/r12 + fxc[rho0])|bb)

and the exchange-like integral

   (ab|(1/r12 + fxc[rho0])|ab)
"""    

from DFTB.MolecularIntegrals.ERIs import electron_repulsion_integral, AtomicDensitySuperposition, angmom_to_xyz
from DFTB.BasisSets import AtomicBasisSet
from DFTB.SlaterKoster import XCFunctionals
from DFTB.AtomicData import atom_names, atomic_number, hartree_to_eV

import numpy as np

def onsite_integrals_unique(atom_name='h', confined=True):
    """
    compute the 5 unique on-site electron integrals for an atom.
    Only s- and p-valence orbitals are considered. The 5 integrals
    are

      (ss|ss)
      (sp|sp)   = (s px|s px) = (s py|s py) = (s pz|s pz)
      (ss|pp)   = (s s |pxpx) = (s s |pypy) = (s s |pzpz)
      (pp|pp)   = (pxpx|pxpx) = (pypy|pypy) = (pzpz|pzpz)
      (pp'|pp') = (pxpy|pxpy) = (pxpz|pxpz) = (pypz|pypz)

    Optional
    --------
    confined      :   controls where confined atoms (True) or free atoms (False)
                      are used

    Returns
    -------
    unique_integrals  :  numpy array with unique on-site integrals
                         [(ss|ss), (sp|sp), (ss|pp), (pp|pp), (pp'|pp')]

    References
    ----------
    [1] http://openmopac.net/manual/1c2e.html
    """
    print "computing unique on-site electron integrals between valence orbitals of '%s'" % atom_name
    print "only s- and p-functions are considered"
    Zat = atomic_number(atom_name)
    atomlist = [(Zat, (0,0,0))]
    basis = AtomicBasisSet(atomlist, confined=confined)
    density = AtomicDensitySuperposition(atomlist, confined=confined)
    #xc_functional = XCFunctionals.libXCFunctional("lda_x", "lda_c_pw")
    xc_functional = XCFunctionals.libXCFunctional("gga_x_pbe", "gga_c_pbe")

    unique_integrals = np.zeros(5)
    for a,bfA in enumerate(basis.bfs):
        for b,bfB in enumerate(basis.bfs):
            # choose index into unique_integrals to which the integrals is written
            if (bfA.l == 0 and bfB.l == 0):
                # (ss|ss)
                i = 0
            elif (bfA.l == 0 and bfB.l == 1 and bfB.m == 0):
                # (sp|sp) and (ss|pp)
                i = 1
            elif (bfA.l == 1 and bfA.m == 0 and bfB.l == 1 and bfB.m == 0):
                # (pp|pp)
                i = 3
            elif (bfA.l == 1 and bfA.m == 0 and bfB.l == 1 and bfB.m == 1):
                # (pp'|pp')
                i = 4
            else:
                continue

            # compute Coulomb integral    (ab|(1/r12 + fxc[rho0])|ab)
            eri_abab = electron_repulsion_integral(atomlist, bfA, bfB, bfA, bfB, density, xc_functional)
            unique_integrals[i] = eri_abab

            if (i == 1): 
                # in addition to (sp|sp) we also need to calculated (ss|pp)
                eri_aabb = electron_repulsion_integral(atomlist, bfA, bfA, bfB, bfB, density, xc_functional)
                unique_integrals[2] = eri_aabb
            elif (i == 4):
                # in addition to (pp'|pp') we also compute (pp|p'p') and verify the identity
                #  (pp|p'p') = (pp|pp) - 2 (pp'|pp')
                eri_aabb = electron_repulsion_integral(atomlist, bfA, bfA, bfB, bfB, density, xc_functional)
                assert (eri_aabb - (unique_integrals[3] - 2*unique_integrals[4])) < 1.0e-5
                
    print "unique 1-center integrals (in eV) for atom '%s'" % atom_name
    print "       (ss|ss)     (sp|sp)     (ss|pp)      (pp|pp)     (pp'|pp')"
    print "    %+.7f  %+.7f  %+.7f  %+.7f  %+.7f" % tuple(unique_integrals*hartree_to_eV)
                
    return unique_integrals
            
            
def onsite_integrals(atom_name='h', confined=True):
    """
    compute on-site Coulomb and exchange integrals for an atom

    Parameters
    ----------
    atom_name      :  name of atom, e.g. 'h', 'c', etc.

    Optional
    --------
    confined      :   controls where confined atoms (True) or free atoms (False)
                      are used
    """
    print "computing on-site Coulomb and exchange integrals between valence orbitals of '%s'" % atom_name
    Zat = atomic_number(atom_name)
    atomlist = [(Zat, (0,0,0))]
    basis = AtomicBasisSet(atomlist, confined=confined)
    density = AtomicDensitySuperposition(atomlist, confined=confined)
    #xc_functional = XCFunctionals.libXCFunctional("lda_x", "lda_c_pw")
    xc_functional = XCFunctionals.libXCFunctional("gga_x_pbe", "gga_c_pbe")
    
    for a,bfA in enumerate(basis.bfs):
        stra = "%d%s" % (bfA.n, angmom_to_xyz[(bfA.l,bfA.m)])  # name of valence orbital, e.g. 1s, 2px, etc.
        for b,bfB in enumerate(basis.bfs):
            strb = "%d%s" % (bfB.n, angmom_to_xyz[(bfB.l,bfB.m)])
            # Coulomb integral    (aa|(1/r12 + fxc[rho0])|bb)
            eri_aabb = electron_repulsion_integral(atomlist, bfA, bfA, bfB, bfB, density, xc_functional)
            print "Coulomb  (%s,%s|{1/r12+f_xc[rho0]}|%s,%s)= %e" % (stra,stra,strb,strb,eri_aabb)
            if a != b:
                # exchange integral    (ab|(1/r12 + fxc[rho0])|ab)
                eri_abab = electron_repulsion_integral(atomlist, bfA, bfB, bfA, bfB, density, xc_functional)
                print "exchange (%s,%s|{1/r12+f_xc[rho0]}|%s,%s)= %e" % (stra,strb,stra,strb,eri_abab)
            

def tabulate_onsite_integrals(atom_names, filename, confined=True):
    """
    compute unique 1-center electron integrals for the list of atoms and save them to 
    a python file that can be imported

    Parameters
    ----------
    atom_names  :   list of atom names, e.g. ['h','c',...]
    filename    :   path to a python file, where the integrals will be stored in a
                    human-readable form

    Optional
    --------
    confined      :   controls where confined atoms (True) or free atoms (False)
                      are used
    """
    import pprint
    import sys
    
    # onsite integrals
    #  (ss|ss)     (sp|sp)     (ss|pp)      (pp|pp)     (pp'|pp')
    # for each atom
    onsite_integrals_dic = {}
    for atom_name in atom_names:
        onsite_integrals_dic[atom_name] = onsite_integrals_unique(atom_name=atom_name, confined=confined)

    fh = open(filename, "w")
    np.set_printoptions(threshold=sys.maxint)
    pp = pprint.PrettyPrinter(depth=10)
    print>>fh, "# This file has been generated automatically by %s." % sys.argv[0]
    print>>fh, "from numpy import array"
    print>>fh, ""
    print>>fh, "# unique on-site electron integrals (ab|1/r12+f_xc[rho0]|cd) in Hartree for each atom"
    print>>fh, "# in the order   (ss|ss)     (sp|sp)     (ss|pp)      (pp|pp)     (pp'|pp')"
    print>>fh, "onsite_integrals = \\\n%s" % pp.pformat(onsite_integrals_dic)

    print "unique on-site integrals for atoms %s written to '%s'" % (" ".join(atom_names), filename)
    fh.close()
    
    
    
if __name__ == "__main__":

    # compute on-site integrals for these atoms
    atom_names = ['h', 'c', 'n', 'o']
    tabulate_onsite_integrals(atom_names, "/tmp/onsite_electron_integrals.py", confined=True)

    
