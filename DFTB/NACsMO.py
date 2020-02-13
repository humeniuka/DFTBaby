#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
approximate non-adiabatic coupling vectors are obtained from eqn. (19) in Ref. [1].
The matrix elements of the coupling vector in the MO basis are approximately
given by
                1                *              0            en(i) + en(j)
     d   = ------------ sum    C    C   [ - d(H    )/dR  + --------------- d(S   )/dR  ]
      ij    en(i)-en(j)    m,n  i,m  j,n       m,n                2           m,n

Coupling vectors between excited states and the ground state are obtained 
by contracting the single-particle matrix elements with the TD-DFTB 
excitation coefficients.

This still is an approximation because the dependence of partial charges and the 
TD-DFTB coefficients on the nuclear geometry is neglected.


References
----------
[1] E. Abad, P. Lewis, V. Zobac, P. Hapala, P. Jelinek, J. Ortega,
    "Calculation of non-adiabatic coupling vectors in a local-orbital basis set"
    J. Chem. Phys. (2013), 138, 154106.

"""
import numpy as np

from DFTB.ExcGradients import gradients_H0andS

def coupling_vectors_kohn_sham(orbeB, orbeK, orbsB, orbsK, gradS, gradH0):
    """
    approximate coupling vectors d_ij between two sets of Kohn-Sham (KS) orbitals

                  d
       d   = <i |---- j>           i=1,...,NorbB       j=1,...,NorbK
        ij        dR

    Parameters
    ----------
    orbeB, orbeK   : Kohn-Sham orbital energies of bra orbitals and ket orbitals
    orbsB, orbsK   : orbsB[:,i] are the MO coefficients of i-th KS orbital
                     in the bra, orbsK[:,j] is are the MO ceofficients of the j-th
                     KS orbital in the ket orbitals
    gradS, gradH0  : gradients of overlap and 0-th order Hamiltonian in AO basis
                     gradS[:,m,n] is the vector d(<m|n>)/dR of length 3*Nat and 
                     similarly for gradH0

    Returns
    -------
    nacvMO         : matrix of shapre (3*Nat,NorbB,NorbK) with non-adiabatic coupling
                     vectors in MO basis between the bra and ket orbitals,
                     nacv[:,i,j] = <i|d/dR|j>
    """
    # transform dH0/dR and dS/dR from AO to MO basis
    #   gradSmo[:,i,j] = sum_m sum_n c^*_{m,i} c_{n,j} gradS[:,m,n]
    gradSmo = np.tensordot(np.tensordot(
        gradS,  orbsB, axes=[1,0]).swapaxes(1,2),
                        orbsK, axes=[2,0])
    #   gradH0mo[:,i,j] = sum_m sum_n c^*_{m,i} c_{n,j} gradH0[:,m,n]
    gradH0mo = np.tensordot(np.tensordot(
        gradH0, orbsB, axes=[1,0]).swapaxes(1,2),
                        orbsK, axes=[2,0])
    # 0.5*(en(i)+en(j))
    avg_orbe = 0.5 * np.add.outer(orbeB, orbeK)
    # en(i)-en(j)
    dif_orbe = np.subtract.outer(orbeB,orbeK)
    # remove elements where en(i) = en(j)
    idif_orbe = 1.0/dif_orbe
    idif_orbe[dif_orbe == 0.0] = 0.0
    # eqn. (19) in ref. [1]
    # 1/(ei-ej) sum_{m,n} c^*_{i,m} c_{j,n} ( - gradH0[:,m,n] + (ei+ej)/2 * gradS[:,m,n] )
    nacvMO = - idif_orbe[np.newaxis,:,:] * gradH0mo  \
             + (idif_orbe * avg_orbe)[np.newaxis,:,:] * gradSmo

    return nacvMO
    
class NACsMO:
    def __init__(self, tddftb):
        self.tddftb = tddftb
        self.dftb = tddftb.dftb2
        
    def activeOrbitals(self):
        orbe = self.dftb.getKSEnergies()
        orbs = self.dftb.getKSCoefficients()
        occ, virt = self.tddftb.getActiveOrbitals()
        nocc = len(occ)
        nvirt = len(virt)
        orbe_occ = orbe[occ]
        orbe_virt = orbe[virt]
        orbs_occ = orbs[:,occ]
        orbs_virt = orbs[:,virt]
        return orbe_occ, orbe_virt, orbs_occ, orbs_virt
    
    def coupling_vectors(self):
        """ 
        coupling vectors between all excited states and the ground state

        Returns
        -------
        nacv    :   nacv[:,j] = <0|d/dR j> is the non-adiabatic coupling vector
                    of length 3*Nat between the excited state j and the ground state
        """
        if self.dftb.verbose > 0:
            print "  gradients of S and H0 from Slater-Koster rules"
        atomlist = self.dftb.getGeometry()
        valorbs = self.dftb.valorbs
        gradS, gradH0 = gradients_H0andS(atomlist, valorbs,
                                         self.dftb.SKT, self.dftb.Mproximity)
        if self.dftb.verbose > 0:
            print "  non-adiabatic coupling vectors in MO basis"
        # NACV's are only computed between active occupied and virtual orbitals
        orbe_occ, orbe_virt, orbs_occ, orbs_virt = self.activeOrbitals()
        nacvMO = coupling_vectors_kohn_sham(orbe_occ, orbe_virt, orbs_occ, orbs_virt,
                                             gradS, gradH0)
        
        # NACV's of many-body electronic states, contract with CIS coefficients
        nacv = np.tensordot(nacvMO, self.tddftb.Cij, axes=([1,2],[1,2]))

        return nacv


