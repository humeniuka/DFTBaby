"""
Natural transition orbitals (NTOs)  

References
----------
[1] Richard Martin, "Natural Transition Orbitals", 
    J. Chem. Phys. (2003), 118, 4775.
"""

import numpy as np
import numpy.linalg as la

def natural_transition_orbitals(T, orbs, occs, virts):
    """
    compute the MO coefficients of the natural transition orbitals

    Parameters
    ----------
    T    :   numpy array of shape (nocc,nvirt),
             transition density matrix for a single excited state I
             in the MO basis:
                             +
                T[i,a] = <I|c_i c_a|0>

    orbs :   Kohn-Sham orbitals, orbs[:,k] are the coefficients of 
             the k-th orbital
    occs :   list of indeces of occupied orbitals
    virts:   list of indeces of virtual orbitals

    Returns:
    --------
    orbs_nto : coefficients of natural transition orbital
                               (  occupied NTOs  for 0 < k < nocc
               orbs_nto[:,k] = {
                               (  virtual NTOs   for nocc <= k
    orbe_nto : eigenvalues of NTOs, occupied and virtual NTOs with
               the same eigenvalue form pairs
    """
    nocc,nvirt = T.shape
    # unitary transformation of occupied orbitals
    # T.T^dagger
    TTt = np.dot(T,T.transpose())
    wu,U = la.eigh(TTt)
    # unitary transformation of virtual orbitals
    # T^dagger.T
    TtT = np.dot(T.transpose(), T)
    wv,V = la.eigh(TtT)
    # sort in decreasing order
    V = V[:,::-1]
    wv = wv[::-1]
    
    # consistency checks:
    # 1) All eigenvalues should be in the range
    #   0 <= wu(i) <= 1
    #   0 <= wv(i) <= 1
    eps = 1.0e-3
    assert np.all((-eps <= wu) & (wu <= 1+eps))
    assert np.all((-eps <= wv) & (wv <= 1+eps))
    # 3) We should be able to form `nocc` pairs of occupied and virtual
    # eigenvalues which are the same:
    #    0 <= wu(i) = wv(nocc-i) <= 1
    err = la.norm(wu[-min(nocc,nvirt):] - wv[:min(nocc,nvirt)][::-1])
    assert err < 1.0e-10
    # For TD-DFT the following sum rule is only approximately fulfilled
    # because of deexcitation 
    #    __nocc-1
    # 4) \        wu(i) = 1
    #    /_i=0
    sum_rule = np.sum(wu)
    print "sum of NTO eigenvalues: %s  (should be approximately 1.0)" % sum_rule

    nao,nmo = orbs.shape
    orbs_nto = np.zeros((nao,nocc+nvirt))
    orbe_nto = np.zeros(nocc+nvirt)
    # transform occupied orbitals
    orbs_nto[:, :nocc] = np.dot(orbs[:,occs], U)
    orbe_nto[:nocc] = wu
    # transform virtual orbitals
    orbs_nto[:, nocc:] = np.dot(orbs[:,virts], V)
    orbe_nto[nocc:] = wv

    return orbs_nto, orbe_nto


        
        
