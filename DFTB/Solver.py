"""
solve the non-Hermitian eigenvalue problem

  ( A  B ) (X)     ( 1  0 ) (X)
  |      |*| | = w |      |*| |
  ( B  A ) (Y)     ( 0 -1 )*(Y)

that appears in linear-response TD-DFT.

# using the iterative algorithm
#proposed by Stratmann, Scuseria and Frisch in J. Chem. Phys. (1998), 109, 19 
# "An efficient implementation of time-dependent density-functional theory
#  for the calculation of excitation energies of large molecules"
"""
import numpy as np
import numpy.linalg as la
import scipy.linalg as sla

# The solver for the TD-DFTB equations with the long-range correction can use a Fortran
# or a python implementation for the matrix multiplication (A+-B).v. The Fortran implementation
# is not faster than the python implementation if numpy is compiled with Intel's MKL.
SOLVER_IMPLEMENTATION = "PY" #"F90" # "PY"

import time
from DFTB.Timer import GlobalTimer as T

class ExcitedStatesNotConverged(Exception):
    pass

class ExcitedStatesError(Exception):
    # Any error that can occur during the solution of the TD-DFTB eigenvalue problem
    # These errors can be caught with an exception handler
    pass

@T.timer
def buildAB(gamma, gamma_lr,\
                qtrans_oo, qtrans_vv, qtrans_ov, \
                omega, omega_shift, df, nocc, nvirt, multiplicity="S"):
    """
    build the matrices A and B. These are the same as for TDHF, but the 2e-integrals are approximated
    by transition charges and the Coulomb potential is screened for matrix elements corresponding to
    particle exchange. The TDHF formulas can be found in 
    Poul Jorgensen, "Molecular and atomic applications of time-dependent Hartree-Fock theory"
    Reference: Jorgensen,P. Annu. Rev. Phys. Chem. (1975), 26, 359-380. equation (58)

      for Singlets:

       A^S_(o,v;o',v') = delta_(o,o') delta_(v,v') (en_v - en_o) + 2*(vo|o'v') - (vv'|o'o)_lr
       B^S_(o,v;o',v') = 2*(ov|o'v') - (vo'|v'o)_lr

      and for Triplets:

       A^T_(o,v;o',v') = delta_(o,o') delta_(v,v') (en_v - en_o) - (vv'|o'o)_lr
       B^T_(o,v;o',v') = - (vo'|v'o)_lr

    So, the triplet matrices only differ by the absence of the electrostatic Coulomb interaction.

      Note that for real orbitals, (vo|o'v') = (ov|o'v') and (vv'|o'o) = (oo'|vv')
      and (vo'|v'o) = (o'v|ov')
    """
    # long range Hartree-Fock contribution to coupling matrix 
    # ... for A matrix, (ij|ab)
    K_lr_A = np.tensordot(qtrans_oo, np.tensordot(gamma_lr, qtrans_vv, axes=(1,0)),axes=(0,0))
    # got K_ij_ab, but we need K_ia_jb
    K_lr_A = np.swapaxes(K_lr_A, 1, 2)
    # ... and for B matrix, (ia|jb)
    K_lr_B = np.tensordot(qtrans_ov, np.tensordot(gamma_lr, qtrans_ov, axes=(1,0)),axes=(0,0))
    # got K_ia_jb but we need K_ib_ja  
    K_lr_B = np.swapaxes(K_lr_B, 1, 3)

    # construct A and B matrices
    K_A = - K_lr_A
    K_B = - K_lr_B

    if multiplicity == "S":
        # coupling matrix 
        # electrostatic Coulomb interaction 2*(ov|o'v'), only needed for Singlets
        K_singlet = 2.0*np.tensordot(qtrans_ov, np.tensordot(gamma, qtrans_ov, axes=(1,0)),axes=(0,0))
        K_A += K_singlet
        K_B += K_singlet

    KcouplingA = np.reshape(K_A, (nocc*nvirt, nocc*nvirt))
    KcouplingB = np.reshape(K_B, (nocc*nvirt, nocc*nvirt))

    dfhalf = np.diag(np.reshape(df/2.0, nocc*nvirt))
    omega = np.diag(np.reshape(omega, nocc*nvirt))

    A = np.dot(dfhalf, omega) + np.dot(dfhalf, np.dot(KcouplingA, dfhalf))          
    B = np.dot(dfhalf, np.dot(KcouplingB, dfhalf))
    
    # add shift from long-range charge transfer correction
    omega_shift = np.diag(np.reshape(omega_shift, nocc*nvirt))
    A += omega_shift
    B += omega_shift

    return A, B

@T.timer
def Casida(gamma, gamma_lr,\
               qtrans_oo, qtrans_vv, qtrans_ov, \
               omega, omega_shift, df, nocc, nvirt, multiplicity="S"):
    """
    solves the non-Hermitian eigenvalue problem for the general case
    where A-B can be diagonal or not

    Parameters:
    ===========
    K: coupling matrix 
    omega: differences of orbital energies, omega[i,a] = en_a - en_i
    omega_shift: shift of diagonal elements of coupling matrix from CT correction
    df: difference in orbital occupation numbers, df[i,a] = f_i - f_a
    nocc, nvirt: number of (active) occupied and virtual orbitals
    multiplicity: "S" for Singlets or "T" for Triplets

    Returns:
    ========
    Omega: 1D numpy array with eigen energies
        Omega[I] is the energy of eigen state I
    Cij: 3D numpy array with coefficients for single excitations,
        Cs[I, i,j] is the coefficient for the excitation from orbital i to orbital j
        in state I.
    """
    A, B = buildAB(gamma, gamma_lr,\
                qtrans_oo, qtrans_vv, qtrans_ov, \
                omega, omega_shift, df, nocc, nvirt, multiplicity=multiplicity)
    
    # check whether A-B is diagonal
    AmB = A-B
    offdiag = la.norm(np.diag(np.diag(AmB)) - AmB)
#    assert offdiag < 1.0e-10, "A-B is not diagonal, error=%s" % offdiag
    AmB_diagonal = False
    if offdiag < 1.0e-10:
#        print "A-B is diagonal"
        sqAmB = np.diag(np.sqrt(np.diag(A-B)))
    else:
#        print "A-B is NOT diagonal"
#        print "calculating matrix square-root  (A-B)^(1/2)"
        sqAmB = sla.sqrtm(A-B)

    # construct hermitian eigenvalue problem
    # see eqn. 16 in J. Chem. Phys., 109, 19 (1998)
    #  (A-B)^(1/2) (A+B) (A-B)^(1/2) F = Omega^2 F
    R = np.dot(sqAmB, np.dot(A+B, sqAmB))
    # R should be Hermitian
    assert la.norm(R-R.conjugate().transpose()) < 1.0e-10, "TD-DFTB matrix not hermitian, probably the SCF cycle did not converge to the right solution, you might try to dissable the DIIS mixer and play around with the option --linear_mixing_coefficient"
    Omega2, F = sla.eigh(R)
    assert la.norm(np.dot(R,F[:,0]) - Omega2[0]*F[:,0]) < 1.0e-10, "maybe BUG in numpy's eigenvalue solver eigh?"
    # sort excitations by energy
    sort_indx = np.argsort(Omega2)
    Omega2 = Omega2[sort_indx]
    F = F[:,sort_indx]
    Omega = np.sqrt(Omega2)

    # compute X-Y and X+Y
    # X+Y = 1/sqrt(Omega) * (A-B)^(1/2).F
    # X-Y = 1/Omega * (A+B).(X+Y)
    XpY = np.dot(sqAmB,F) / np.sqrt(Omega)
    ApB = A+B
    XmY = np.dot(ApB, XpY) / Omega

    err = abs(np.sum(XpY[:,0]*XmY[:,0]) - 1.0)
    assert err < 1.0e-10, "err = %s" % err
    err = np.sum(abs(np.dot(ApB, XpY) - Omega*XmY))
    assert err < 1.0e-5, "err = %s" % err


    # C = (A-B)^(-1/2).((X+Y) * sqrt(Omega))
    # so that C^T.C = (X+Y)^T.(A-B)^(-1).(X+Y) * Omega
    #               = (X+Y)^T.(X-Y)
    # since (A-B).(X-Y) = Omega * (X+Y)
    C = la.solve(sqAmB, XpY * np.sqrt(Omega))
    err = abs(np.sum(C[:,0]*C[:,0]) - 1.0)
    assert err < 1.0e-10, "err = %s" % err

#
    shape = (nocc*nvirt, nocc, nvirt)
    XmY = np.reshape(XmY.transpose(), shape)
    XpY = np.reshape(XpY.transpose(), shape)
    C = np.reshape(C.transpose(), shape)

    return Omega, C, XmY, XpY

@T.timer
def TDA(gamma, gamma_lr,\
            qtrans_oo, qtrans_vv, qtrans_ov, \
            omega, df, nocc, nvirt):
    """
    perform the Tamm-Dancoff approximation and solve

      A*X = w*X

    instead of the non-Hermitian eigenvalue problem.
    """
    raise Exception("Revise the TDA code!!!!")
    # build coupling matrix
    K_singlet = np.tensordot(qtrans_ov, np.tensordot(gamma, qtrans_ov, axes=(1,0)),axes=(0,0))
    # long range Hartree-Fock contribution to coupling matrix
    K_lr = np.tensordot(qtrans_oo, np.tensordot(gamma_lr, qtrans_vv, axes=(1,0)),axes=(0,0))
    # got K_ij_ab, but we need K_ia_jb
    K_lr = np.swapaxes(K_lr, 1, 2)
    K = K_singlet - K_lr
    Kcoupling = np.reshape(K, (nocc*nvirt, nocc*nvirt))

    dfhalf = np.diag(np.reshape(df/2.0, nocc*nvirt))
    omega = np.diag(np.reshape(omega, nocc*nvirt))
    # construct A matrix
    #        A = omega + Kcoupling
    A = np.dot(dfhalf, omega) + np.dot(dfhalf, np.dot(Kcoupling, dfhalf))          
    # see formula (83) in Jorgensen (1975) "Molecular and Atomic Applications of Time-Dependent Hartree-Fock Theory"
    Omega, X = la.eigh(A)
    sort_indx = np.argsort(Omega)
    # CIS excitation coefficients for excitation (occupied i) -> (virtual j)
    Cij = np.reshape(X.transpose(),(nocc*nvirt,nocc,nvirt))
    Cij = Cij[sort_indx,:,:]
    assert abs(np.sum(Cij[0,:,:]*Cij[0,:,:])-1.0) < 1.0e-10

    return Omega, Cij


"""
implementation of Stratmann's algorithm for solving the non-hermitian
TD-DFT eigenvalue problem iteratively

see  
     Stratmann,~E., Scuseria,~G., Frisch,~M.
    "An efficient implementation of time-dependent density-functional theory 
     for the calculation of excitation energies of large molecules."
     J. Chem. Phys., 109, 8218, (1998).
"""
from DFTB.extensions import tddftb

def initial_expansion_vectors(omega, lmax):
    """
    The initial guess vectors are the lmax lowest energy
    single excitations
    """
    nocc,nvirt = omega.shape
    bs = np.zeros((nocc,nvirt,lmax))
    omega_flat = np.ravel(omega)
    sort_indx = np.argsort(omega_flat)
    for l in range(0, lmax):
        idx = sort_indx[l]
        # row - occupied index
        i = int(idx / nvirt)
        # col - virtual index
        a = idx % nvirt
        assert i*nvirt+a == idx
        bs[i,a,l] = 1.0 
#        print "%d -> %d   en=%s   en=%s" % (i,a, omega_flat[idx], omega[i,a])
#        assert omega_flat[idx] == omega[i,a]
    return bs

def linear_independent_cols(X, singular_threshold):
    """
    extract a set of columns of a matrix X
    such that the subset of columns is linear independent

    Parameters:
    ===========
    X: input matrix
    singular_threshold: threshold of treating small numbers as zero

    Returns:
    ========
    indx: list of N indeces into columns of X
    """
    # Q: unitary matrix
    # R: upper triangular matrix with non-increasing diagonal elements
    # P: permutation matrix
    # X = Q R P^T
    Q,R,P = sla.qr(X, pivoting=True)
    # index of smallest diagonal element larger than zero_thresh
    N = np.where(abs(np.diag(R)) > singular_threshold)[0][-1]+1
    indx = P[:N]
    indx.sort()

    return indx

def norm(v):
    # parallelized dot in Fortran
#    v2 = tddftb.tddftb.dot(v,v)
    # python
    v2 = np.tensordot(v,v,axes=([0,1],[0,1]))
    return np.sqrt(v2)
def dot(v,w):
    # parallelized dot in Fortran
#    return tddftb.tddftb.dot(v,w)
    # python
    return np.tensordot(v,w,axes=([0,1],[0,1]))

def gram_schmidt_2d(W):
    """
    orthogonalize a set of vectors vs using the Gram-Schmidt procedure.
    W holds all the vectors and has the shape (nocc,nvirt,nvectors)
    A single vector w has the shape (nocc,nvirt)

    Parameters:
    ===========
    W: set of linear independent vectors
    
    Returns:
    ========
    V: set of orthonormalized vectors
    """
    def olap(W):
        dim = W.shape[-1]
        # vectors W are not normalized
        norms = np.zeros(dim)
        for i in range(0, dim):
            norms[i] = np.sqrt(dot(W[:,:,i],W[:,:,i]))
        S = np.zeros((dim,dim))
        for i in range(0, dim):
            for j in range(0, dim):
                S[i,j] = dot(W[:,:,i],W[:,:,j])/(norms[i]*norms[j])
        return S
    ## check that expansion vectors are linearly independent
    S = olap(W)
    detS = la.det(S)
    if abs(detS) == 0.0:
        print "Trying to remove linear dependent expansion vectors ...",
        indx = linear_independent_cols(S, 1.0e-18)
        print "done"
        W = W[:,:,indx]
        S = olap(W)
#    assert la.det(S) > 0.0, "Gram-Schmidt orthogonalizer got vectors which are not linearly independent!" 
    ##
    dim = W.shape[-1]
    # V holds the orthonormalized vectors
    V = np.zeros(W.shape)
    V[:,:,0] = W[:,:,0]/norm(W[:,:,0])
    for n in range(1, dim):
        V[:,:,n] = np.copy(W[:,:,n])
        for m in range(0,n):
#            V[:,:,n] -= dot(V[:,:,m],W[:,:,n])*V[:,:,m] # numerically unstable classical Gram-Schmidt 
            V[:,:,n] -= dot(V[:,:,m],V[:,:,n])*V[:,:,m] # numerically stable modified Gram-Schmidt
        V[:,:,n] /= norm(V[:,:,n])
    # check that all vectors are orthonormal
    S = np.zeros((dim,dim))
    for i in range(0, dim):
        for j in range(0, dim):
            S[i,j] = dot(V[:,:,i],V[:,:,j])
    err = np.sum(abs(S - np.eye(dim)))
#    assert err < 1.0e-10, "Gram-Schmidt orthogonalization failed, error = %s" % err
    if err > 1.0e-10:
        print "WARNING: Gram-Schmidt orthogonalization failed, error = %s" % err
    return V

def reorder_vectors_Lambda2(Oia, w2, T, L2_thresh=0.4):
    """
    reorder the expansion vectors so that those with Lambda2 values
    above a certain threshold come first
    """
    # compute Lambda2
    nocc,nvirt,nst = T.shape
    L2 = np.zeros(nst)
#    print "Lambda2"
    for n in range(0, nst):
        L2[n] = np.tensordot(T[:,:,n]**2,Oia, axes=([0,1],[0,1]))
#        print "%d      %2.5f" % (n, L2[n])
    # indeces of
    good = np.where(L2 > L2_thresh)[0]
    bad = np.where(L2 <= L2_thresh)[0]
    Tnew = np.zeros((nocc,nvirt,nst))
    w2new = np.zeros(nst)
    # place good states at the beginning
    for i,ig in enumerate(good):
        Tnew[:,:,i] = T[:,:,ig]
        w2new[i] = w2[ig]
    # bad states follow
    ng = len(good)
    for i,ib in enumerate(bad):
        Tnew[:,:,i+ng] = T[:,:,ib]
        w2new[i+ng] = w2[ib]
        
    return w2new, Tnew

def _reorder_vectors(w2, T, good, bad, labels, nstates):
    """
    reorder eigenvalues and eigenvectors so that those with indeces
    in good come first
    """
    Tnew = np.zeros(T.shape)
    w2new = np.zeros(w2.shape)
    labels_reordered = []
    #place good states at the beginning
    for i,ig in enumerate(good):
        Tnew[...,i] = T[...,ig]
        w2new[i] = w2[ig]
        labels_reordered.append( labels[ig] )
    #bad states follow
    ng = len(good)
    for i,ib in enumerate(bad):
        Tnew[...,i+ng] = T[...,ib]
        w2new[i+ng] = w2[ib]
        labels_reordered.append( labels[ib] )
        
    nst = len(good)
    if nst < nstates:
        nst += len(bad)
    #print "ngood = %d   nbad = %s" % (len(good), len(bad))
    #print "nst = %s   nstates = %s" % (nst, nstates)
    return w2new[...,:nst], Tnew[...,:nst], labels_reordered[:nst]

@T.timer
def HermitianDavidson(g0, q_ov, \
                omega, omega_shift, nocc, nvirt, \
                XmYguess, XpYguess, \
                Oia, \
                selector=None, \
                multiplicity="S", \
                nstates=4, ifact=1, conv=1.0e-14, \
                maxiter=10, \
                L2_thresh=0.5, verbose=1):
    """
    If A-B is diagonal the TD-DFT equations can be made hermitian

      (A-B)^(1/2).(A+B).(A-B)^(1/2).T = Omega^2 T
                   R               .T = Omega^2 T

    """
    if verbose > 0:
        print "dimension of full vector space: %d" % (nocc*nvirt)
    omega2 = omega**2
    omega_sq = np.sqrt(omega)
    omega_sq_inv = 1.0/omega_sq
    # wq_ov = sqrt(orbe_a - orbe_i) * q_A^(ia)
    wq_ov = q_ov * omega_sq
    # diagonal elements of R
    om = omega2+2*omega*omega_shift
    # matrix-vector product
    def Rv(vs):
        lmax = vs.shape[2]
        us = np.zeros(vs.shape)
        for l in range(0, lmax):
            v = vs[:,:,l]
            # matrix product u = sum_jb (A-B)^(1/2).(A+B).(A-B)^(1/2).v
            # 1st term in (A+B).v  - KS orbital energy differences
            u = om*v
            # 2nd term - Coulomb
            tmp = np.tensordot(wq_ov, v, axes=([1,2],[0,1]))
            tmp2 = np.dot(g0, tmp)
            u += 4*np.tensordot(wq_ov, tmp2, axes=(0,0))
            us[:,:,l] = u
        return us

    def writeIteration(it, norms, labels, w, l):
        dt = time2-time1
        print "Iteration %d: ( %d expansion vectors )   time: %2.6f s" % (it, l, dt)
        nst = len(norms)
        for n in range(0, nst):
            conv_str =     "not converged"
            if norms[n] < conv:
                conv_str = "converged    "
            print "  Eigenvalue %3d      %2.7f    %s   res. norm = %e     %s" % (n+1, w[n], conv_str, norms[n], labels[n])


    # initial number of expansion vectors
    kmax = nocc*nvirt  # at most there are nocc*nvirt excited states
    lmax = min(ifact*nstates,kmax)  # nocc*nvirt
    #
    if type(XpYguess) == type(None):
        # If there are many low-lying CT states, then using
        # the orbital energy differences to select the initial vectors
        # is a bad choice. In this case the guess energies need to be 
        # corrected
        omega_guess = np.sqrt(om)
        bs = initial_expansion_vectors(omega_guess, lmax)
    else:
        if verbose > 0:
            print "start with vectors from previous calculation"
        # initialize expansion vectors from guess vectors
        # If there are degenerate states, the algorithm still converges
        # slowly although the start vectors are very close the
        # final solution.
        bs = np.zeros((nocc,nvirt,lmax))
        for n in range(0,lmax):
            bs[:,:,n] = omega_sq_inv * XpYguess[n,:,:]
            bs[:,:,n] /= norm(bs[:,:,n])

    l = lmax  # l is the current number of expansion vectors
    k = nstates
    for it in range(0, maxiter):
        time1 = time.time()
        rbs = Rv(bs)
        Hb = np.tensordot(bs,rbs, axes=([0,1],[0,1]))
        w2, Tb = la.eigh(Hb)
        # transform to the canonical basis Tb -> T
        T = np.tensordot(bs, Tb, axes=(2,0))
        if selector != None:
            good, bad, labels = selector.select_vectors(T)
            # place states with the selected properties first
            w2, T, labels = _reorder_vectors(w2, T, good, bad, labels, nstates)
            # number of expansion vectors that have the right properties
            ngood = len(good)
        else:
            ngood = len(w2)
            labels = ["" for n in range(0,ngood)]
        assert ngood > 0, "Not enough states with from the right irrep. Try increasing diag_ifact to include more initial guess vectors."
        # sort vectors so that those with low Lambda2 values come first
        if L2_thresh > 0.0:
            w2, T = reorder_vectors_Lambda2(Oia, w2, T, L2_thresh=L2_thresh)
        #
        w = np.sqrt(w2)
        # residual vectors
        W = Rv(T) - w2*T
        norms = np.zeros(k)
        for n in range(0,k):
            norms[n] = norm(W[:,:,n])
        time2 = time.time()
        if verbose > 0:
            writeIteration(it, norms, labels, w, l)
        # check for convergence
        if (norms < conv).all():
            if verbose > 0:
                print "All %d roots CONVERGED (residual norms < %s) after %d iterations" \
                    % (k, conv, it+1)
            break
        # enlarge dimension of subspace by dk vectors
        # At most k new expansion vectors are added
        dkmax = min(kmax-l,k) # if the full space is reached, stop enlarging
        # count number of non-converged vectors
        # residual vectors that are zero cannot be used as new expansion vectors
        eps = 0.01*conv #1.0e-16
        nc = np.sum(norms > eps)
        dk = min(dkmax,nc)
        #
        Qs = np.zeros((nocc, nvirt, dk))
        nb = 0   # enumerate vectors that are added 
        # select new expansion vectors among the residual vectors
        for n in range(0, dkmax):
            wD = w[n] - omega
            # remove singularities from 1/wD which lead to linear dependencies
            indx = abs(wD) < 1.0e-6
            wD[indx] = 1.0e-6 * omega[indx]
            if norms[n] > eps:
                Qs[:,:,nb] = 1.0/wD * W[:,:,n]
                nb += 1
        assert nb==dk
        # new expansion vectors are bs + Qs
        bs_new = np.zeros((nocc,nvirt, l+dk))
        bs_new[:,:,:l] = bs
        bs_new[:,:,l:] = Qs
        # orthonormalize all expansion vectors among themselves
        if verbose > 0:
            print "orthogonalize new expansion vectors",
        #bs = gram_schmidt_2d(bs_new)
        # Orthogonalize using QR decomposition, which effectively performs a Gram-Schmidt orthogonalization
        nvec = l+dk
        bs_flat = np.reshape(bs_new, (nocc*nvirt, nvec))
        Q,R = sla.qr(bs_flat, mode='economic')
        bs = np.reshape(Q, (nocc,nvirt,nvec))
        if verbose > 0:
            print "... done"
        #
        l = bs.shape[-1]   # l+dk unless linear dependent expansion vectors were removed
    else:
        raise ExcitedStatesNotConverged("Davidson-like diagonalization for Hermitian eigenvalue problem did not converge after %s iterations (threshold: %e) !!! It might help to include more states or increase the number of iterations using the --diag_maxiter=100 option" % (it+1, conv))

    #
    Omega = w[:k]
    XpY = np.zeros((nocc,nvirt,k))
    XmY = np.zeros((nocc,nvirt,k))
    C = np.zeros((nocc,nvirt,k))
    for n in range(0,k):
        # X+Y = 1/sqrt(Omega)*(A-B)^(1/2).T
        XpY[:,:,n] = omega_sq/np.sqrt(Omega[n]) * T[:,:,n]
        # X-Y = sqrt(Omega)*(A-B)^(-1).(X+Y)
        XmY[:,:,n] = np.sqrt(Omega[n]) * omega_sq_inv * T[:,:,n]
        # C = (A-B)^(-1/2).(X+Y) * sqrt(Omega)
        C[:,:,n] = T[:,:,n]

    # XmY, XpY and C have shape (nocc,nvirt, nstates)
    # bring the last axis to the front
    XmY = np.rollaxis(XmY, 2)
    XpY = np.rollaxis(XpY, 2)
    C = np.rollaxis(C, 2)

    return Omega, C, XmY, XpY


# matrix products for iterative solvers

def _ApBv(g0, g0_lr, q_oo, q_ov, q_vv, omega, vs, lc):
    lmax = vs.shape[2]
    us = np.zeros(vs.shape)
    for l in range(0, lmax):
        time1 = time.time()
        v = vs[:,:,l]
        # matrix product u_ia = sum_jb (A+B)_(ia,jb) v_jb
        # 1st term in (A+B).v  - KS orbital energy differences
        u = omega*v
        # 2nd term - Coulomb
        tmp = np.tensordot(q_ov, v, axes=([1,2],[0,1]))
        tmp2 = np.dot(g0, tmp)
        u += 4*np.tensordot(q_ov, tmp2, axes=(0,0))
        if lc == 1:
            # 3rd term - Exchange
            tmp = np.tensordot(q_vv, v, axes=(2,1))
            tmp2 = np.tensordot(g0_lr, tmp, axes=(1,0))
            u -= np.tensordot(q_oo, tmp2, axes=([0,2],[0,2]))
            # 4th term - Exchange
            tmp = np.tensordot(q_ov, v, axes=(1,0))
            tmp2 = np.tensordot(g0_lr, tmp, axes=(1,0))
            u -= np.tensordot(q_ov, tmp2, axes=([0,2],[0,2]))
        us[:,:,l] = u
        time2 = time.time()
        dtime=time2-time1
        if dtime > 2.5:  # print if it takes more than 2.5 seconds per expansion vector
            print "(A+B).v for expansion vectors nr. %4d   time = %10.5f seconds" % (l, dtime)
    return us

def _AmBv(g0, g0_lr, q_oo, q_ov, q_vv, omega, vs, lc):
    lmax = vs.shape[2]
    us = np.zeros(vs.shape)
    for l in range(0, lmax):
        time1 = time.time()
        v = vs[:,:,l]
        # matrix product u_ia = sum_jb (A-B)_(ia,jb) v_jb
        # 1st term, differences in orbital energies
        u = omega*v
        if lc == 1:
            # 2nd term, Exchange
            tmp = np.tensordot(q_ov, v, axes=(1,0))
            tmp2 = np.tensordot(g0_lr, tmp, axes=(1,0))
            u += np.tensordot(q_ov, tmp2, axes=([0,2],[0,2]))
            # 3rd term, Exchange
            tmp = np.tensordot(q_vv, v, axes=(2,1))
            tmp2 = np.tensordot(g0_lr, tmp, axes=(1,0))
            u -= np.tensordot(q_oo, tmp2, axes=([0,2],[0,2]))
        us[:,:,l] = u
        time2 = time.time()
        dtime=time2-time1
        if dtime > 2.5:  # print if it takes more than 2.5 seconds per expansion vector
            print "(A-B).v for expansion vectors nr. %4d   time = %10.5f seconds" % (l, dtime)
    return us

@T.timer
def nonHermitianDavidson(g0, g0_lr,\
                q_oo, q_vv, q_ov, \
                omega, nocc, nvirt, \
                XmYguess, XpYguess, w_guess, \
                selector=None, \
                multiplicity="S", \
                nstates=4, ifact=1, conv=1.0e-14, \
                maxiter=10, verbose=1, lc=1):
    """
    Parameters:
    ===========
    lc: flag to indicate whether the long-range correction should be added or not,
       if g0_lr is zero, this flag saves computation time
    """
    if verbose > 0:
        print "dimension of full vector space: %d" % (nocc*nvirt)
    assert multiplicity == "S"
    # matrix products
    def ApBv(vs):
        """
        # check that python and Fortran product the same result
        import time
        # with Fortran
        ta = time.time()
        bp_f90 = tddftb.tddftb.apbv(g0, g0_lr, q_oo, q_ov, q_vv, omega, vs)
        tb = time.time()
        print "(A+B).v with Fortran took %s seconds" % (tb-ta)
        # with python
        ta = time.time()
        bp_py = _ApBv(g0, g0_lr, q_oo, q_ov, q_vv, omega, vs, lc)
        tb = time.time()
        print "(A+B).v with python  took %s seconds" % (tb-ta)
        err = np.sum(abs(bp_f90-bp_py))
        assert err < 1.0e-10, "err = %s" % err
        """
        if SOLVER_IMPLEMENTATION == "F90":
            # matrix multiplication is coded in Fortran and
            # parallelized with OpenMP. Because the intrisic MATMUL is very slow, this will
            # not be fast than the python code.
            # This part of the code can cause segmentation fault if the stack is too small
            bp = tddftb.tddftb.apbv(g0, g0_lr, q_oo, q_ov, q_vv, omega, vs)
        else:
            # python implementation
            bp = _ApBv(g0, g0_lr, q_oo, q_ov, q_vv, omega, vs, lc)
        return bp
    def AmBv(vs):
        if SOLVER_IMPLEMENTATION == "F90":
            # matrix multiplication is coded in Fortran and
            # parallelized with OpenMP. Because the intrisic MATMUL is very slow, this will
            # not be faster than the python code.
            # This part of the code can cause segmentation fault if the stack is too small
            bm = tddftb.tddftb.ambv(g0, g0_lr, q_oo, q_ov, q_vv, omega, vs)
        else:
            # python implementation
            bm = _AmBv(g0, g0_lr, q_oo, q_ov, q_vv, omega, vs, lc)
        return bm
    def writeIteration(it, norms, labels, w, l):
        dt = time2-time1
        print "Iteration %d: ( %d expansion vectors )   time: %2.6f s" % (it, l, dt)
        nst = len(norms)
        for n in range(0, nst):
            conv_str =     "not converged"
            if norms[n] < conv:
                conv_str = "converged    "
            print "  Eigenvalue %3d      %2.7f    %s   res. norm = %e     %s" % (n+1, w[n], conv_str, norms[n], labels[n])

    kmax = nocc*nvirt  # at most there are nocc*nvirt excited states
    # To achieve fast convergence the solution vectors from a nearby geometry
    # should be used as initial expansion vectors
    if XpYguess is None:
        # initial number of expansion vectors
        lmax = min(ifact*nstates,kmax)  # nocc*nvirt
        bs = initial_expansion_vectors(omega, lmax)
    else:
        if verbose > 0:
            print "start with vectors from previous calculation"
        """
        # use solution vectors from previous iteration
        lmax = nstates
        assert XmYguess.shape[0] == nstates
        assert XpYguess.shape[0] == nstates
        bs = np.zeros((nocc,nvirt,2*lmax))
        # bring axis with state indeces to the back
        bs[:,:,:lmax] = np.rollaxis(XmYguess, 0,2+1)
        bs[:,:,lmax:] = np.rollaxis(XpYguess, 0,2+1)
        # XmY and XpY are not normalized individually
        for n in range(0, 2*lmax):
            bs[:,:,n] /= norm(bs[:,:,n])
        #bs = gram_schmidt_2d(bs)
        # To span the space of both X+Y and X-Y we need twice as many
        # expansion vectors
        lmax *= 2
        """
        # use solution vectors from previous iteration
        assert XmYguess.shape[0] == nstates
        assert XpYguess.shape[0] == nstates

        # For the expansion vectors we use X+Y
        lmaxp = nstates
        # and X-Y if the vectors space has not been exhausted yet
        lmaxm = min(nstates,kmax-nstates)
        lmax = lmaxp+lmaxm
        
        bs = np.zeros((nocc,nvirt,lmax))
        # bring axis with state indeces to the back
        bs[:,:,:lmaxp] = np.rollaxis(XmYguess, 0,2+1)
        bs[:,:,lmaxp:] = np.rollaxis(XpYguess, 0,2+1)[:,:,:lmaxm]
        # XmY and XpY are not normalized individually
        for n in range(0, lmax):
            bs[:,:,n] /= norm(bs[:,:,n])
        #bs = gram_schmidt_2d(bs)

    l = lmax  # l is the current number of expansion vectors
    k = nstates
    for it in range(0, maxiter):
        time1 = time.time()
        if (XpYguess is None) or it > 0:
            # evaluate (A+B).b and (A-B).b
            bp = ApBv(bs)
            bm = AmBv(bs)
            # M^+ = (b_i, (A+B).b_j)
            Mp = np.tensordot(bs, bp, axes=([0,1],[0,1]))
            # M^- = (b_i, (A-B).b_j)
            Mm = np.tensordot(bs, bm, axes=([0,1],[0,1]))
            Mmsq = sla.sqrtm(Mm)
            # Mh is the analog of (A-B)^(1/2).(A+B).(A-B)^(1/2) 
            # in the reduced subspace spanned by the expansion vectors bs
            Mh = np.dot(Mmsq, np.dot(Mp, Mmsq))
            ### check that Mh is hermitian
            err = np.sum(abs(Mh - Mh.conjugate().transpose()))
            #        assert err < 1.0e-10, "Mh is not hermitian! error = %s" % err
            if err > 1.0e-10:
                raise ExcitedStatesError("Mh is not hermitian! error = %s" % err)
            ###
            w2, Tb = la.eigh(Mh)
            if selector != None:
                # transform Tb into the canonical basis
                T = np.tensordot(bs, Tb, axes=(2,0))
                good, bad, labels = selector.select_vectors(T)
                # place states with the selected properties first
                # and remove the others if there are enough good states
                w2, Tb, labels = _reorder_vectors(w2, Tb, good, bad, labels, nstates)
                # number of expansion vectors that have the right properties
                ngood = len(good)
            else:
                ngood = len(w2)
                labels = ["" for n in range(0,ngood)]
            assert ngood > 0, "Not enough states with from the right irrep. Try increasing diag_ifact to include more initial guess vectors."
            w = np.sqrt(w2)
            wsq = np.sqrt(w)
            # approximate right R = (X+Y) and left L = (X-Y) eigenvectors
            # in the basis bs
            # (X+Y) = (A-B)^(1/2).T / sqrt(w)
            Rb = np.dot(Mmsq, Tb) / wsq
            # L = (X-Y) = 1/w * (A+B).(X+Y)
            Lb = np.dot(Mp,Rb) / w
            ### check that (Lb^T, Rb) = 1 is fulfilled
            err = np.sum(abs(np.dot(Lb.transpose(), Rb) - np.eye(len(w))))
            #        assert err < 1.0e-6, "(X+Y) and (X-Y) vectors not orthonormal! error = %s" % err
            if err > 1.0e-3:
                raise ExcitedStatesError("(X+Y) and (X-Y) vectors not orthonormal! error = %s" % err)
            ###
            # transform to the canonical basis Lb -> L, Rb -> R
            # Rn = sum_m^l Rb_(m,n) |b_(i,a,m)>
            #        R = tddftb.tddftb.tensordot20(bs,Rb.real)
            # Ln = sum_m^l Lb_(m,n) |b_(i,a,m)>
            #        L = tddftb.tddftb.tensordot20(bs,Lb.real)
            ## check
            L = np.tensordot(bs, Lb, axes=(2,0))
            R = np.tensordot(bs, Rb, axes=(2,0))
        else:
            # use solution vectors from previous iteration
            R = np.rollaxis(XpYguess, 0,2+1)
            L = np.rollaxis(XmYguess, 0,2+1)
            w = w_guess
            labels = ["" for n in range(0,len(w))]
#        L_py = np.tensordot(bs, Lb, axes=(2,0))
#        R_py = np.tensordot(bs, Rb, axes=(2,0))
#        err = np.sum(abs(R-R_py)) + np.sum(abs(L-L_py))
#        assert err < 1.0e-10, "err = %s" % err
        ### check that (L^T, R) = 1 is fulfilled
#        err = np.sum(abs(np.tensordot(L, R, axes=([0,1],[0,1])) - np.eye(l)))
#        assert err < 1.0e-6, "(X+Y) and (X-Y) vectors not orthonormal! error = %s" % err
        ### 
        # residual vectors
        WL = ApBv(R) - L*w
        WR = AmBv(L) - R*w
        #
        norms = np.zeros(k)
        normsL = np.zeros(k)
        normsR = np.zeros(k)
        for n in range(0, k):
            normsL[n] = norm(WL[:,:,n])
            normsR[n] = norm(WR[:,:,n])
            norms[n] =  normsL[n] + normsR[n]
        time2 = time.time()
        if verbose > 0:
            writeIteration(it, norms, labels, w, l)
        # check for convergence
        if (norms < conv).all() and it > 0:
            if verbose > 0:
                print "All %d roots CONVERGED (residual norms < %s) after %d iterations" \
                    % (k, conv, it+1)
            break
        # enlarge dimension of subspace by dk vectors
        # At most 2*k new expansion vectors are added
        dkmax = min(kmax-l,2*k) # if the full space is reached, stop enlarging
        # count number of non-converged left and right vectors
        # residual vectors that are zero cannot be used as new expansion vectors
        eps = 0.01*conv #1.0e-16
        ncL = np.sum(normsL > eps)
        ncR = np.sum(normsR > eps)
        # Half the new expansion vectors should come from the left residual vectors
        # the other half from the right residual vectors. 
        dkR = min(int(np.ceil(dkmax/2.0)), ncL)
        dkL = min(dkmax-dkR, ncR)
        dk = dkR+dkL
        #
        Qs = np.zeros((nocc, nvirt, dk))
        nb = 0   # enumerate vectors that are added 
        # select new expansion vectors among the non-converged left residual vectors
        for n in range(0, k):
            if nb == dk:
                # got enough new expansion vectors
                break
            wD = w[n] - omega
            # remove singularities from 1/wD which lead to linear dependencies
            indx = abs(wD) < 1.0e-6
            wD[indx] = 1.0e-6 * omega[indx]
            if normsL[n] > eps:
                #print "add L residual vector %d" % n
                Qs[:,:,nb] = 1.0/wD * WL[:,:,n]
                nb += 1
        # select the others among the non-converged right residual vectors
        for n in range(0, k):
            if nb == dk:
                # got enough new expansion vectors
                break
            wD = w[n] - omega
            # remove singularities from 1/wD which lead to linear dependencies
            indx = abs(wD) < 1.0e-6
            wD[indx] = 1.0e-6 * omega[indx]
            if normsR[n] > eps:
                #print "add R residual vector %d" % n
                Qs[:,:,nb] = 1.0/wD * WR[:,:,n]
                nb += 1
        assert nb==dk, "nb = %d vectors were added, but dk = %d vectors should have been added" % (nb,dk)
        # new expansion vectors are bs + Qs
        bs_new = np.zeros((nocc,nvirt, l+dk))
        bs_new[:,:,:l] = bs
        bs_new[:,:,l:] = Qs
        # orthonormalize all expansion vectors among themselves
        if verbose > 0:
            print "orthogonalize new expansion vectors",
        #bs = gram_schmidt_2d(bs_new)
        # Orthogonalize using QR decomposition, which effectively performs a Gram-Schmidt orthogonalization
        nvec = l+dk
        bs_flat = np.reshape(bs_new, (nocc*nvirt, nvec))
        Q,R = sla.qr(bs_flat, mode='economic')
        bs = np.reshape(Q, (nocc,nvirt,nvec))
        if verbose > 0:
            print "... done"
        #
        l = bs.shape[-1]   # l+dk unless linear dependent expansion vectors were removed
    else:
        raise ExcitedStatesNotConverged("Davidson-like diagonalization for non-Hermitian eigenvalue problem did not converge after %s iterations (threshold: %e) !!! It might help to include more states or increase the number of iterations using the --diag_maxiter=100 option." % (it+1, conv))
    
    Omega = w[:k]
    XpY = R[:,:,:k]
    XmY = L[:,:,:k]
    # cis coefficients for 'wavefunction' in the basis bs
    T = np.tensordot(bs, Tb, axes=(2,0))
    C = T[:,:,:k]
    
    # XmY, XpY and C have shape (nocc,nvirt, nstates)
    # bring the last axis to the front
    XmY = np.rollaxis(XmY, 2)
    XpY = np.rollaxis(XpY, 2)
    C = np.rollaxis(C, 2)
    return Omega, C, XmY, XpY

################### LINEAR SOLVERS for Z-vector and CPKS equations #################
    
def Krylov_solver_Zvector(A, Adiag, B, X0, maxiter=1000, conv=1.0e-14, verbose=1):
    """
    solve the matrix equation A.X = B
    iteratively. This function is used for solving the Z-vector equations

    Parameters:
    ===========
    A: linear operator, such that A(X) = A.X
    Adiag: diagonal elements of A-matrix, with dimension (nocc,nvirt)
    B: right hand side of equation, (nocc,nvirt, k)
    X0: initial guess vectors or None
    """
    def writeIteration(it, norms, l):
        dt = time2-time1
        print "Iteration %d: ( %d expansion vectors )   time: %2.6f s" % (it, l, dt)
        k = len(norms)
        for n in range(0, k):
            conv_str =     "not converged"
            if norms[n] < conv:
                conv_str = "converged    "
            print "  Solution vector  %3d      %s   res. norm = %e" % (n, conv_str, norms[n])
    if verbose > 0:
        print "Solve Z-vector equation iteratively"
    # number of vectors
    nocc,nvirt,k = B.shape
    kmax = nocc*nvirt
    l = k
    # bs are the expansion vectors
    # initial guess
    Ainv = 1.0/Adiag # approximate inverse of A by inverse of diagonal
    if type(X0) == type(None):
        # 
        bs = np.zeros((nocc,nvirt,k))
        for n in range(0, k):
            bs[:,:,n] = Ainv * B[:,:,n]
    else:
        # start from solution of previous calculation
        if verbose > 0:
            print "start with solution from previous calculation"
        bs = X0
    #
    for it in range(0, maxiter):
        time1 = time.time()
        # representation of A in the basis of expansion vectors
        Ab = np.tensordot(bs, A(bs), axes=([0,1],[0,1]))
        # RHS in basis of expansion vectors
        Bb = np.tensordot(bs, B, axes=([0,1],[0,1]))
        # solve
        Xb = la.solve(Ab, Bb)
        # transform solution vector back into canonical basis
        X = np.tensordot(bs, Xb, axes=(2,0))
        # residual vectors
        W = A(X) - B
        norms = np.zeros(k)
        for n in range(0, k):
            norms[n] = norm(W[:,:,n])
        time2 = time.time()
        if verbose > 0:
            writeIteration(it, norms, l)
        # check for convergence
        if (norms < conv).all():
            if verbose > 0:
                print "Z-vector equations CONVERGED (residual norms < %s) after %d iterations" \
                    % (conv, it+1)
            break
        # enlarge dimension of subspace by dk vectors
        # At most k new expansion vectors are added
        dkmax = min(kmax-l,k) # if the full space is reached, stop enlarging
        # count number of non-converged vectors
        # residual vectors that are zero cannot be used as new expansion vectors
        eps = 0.01*conv #1.0e-16
        nc = np.sum(norms > eps)
        dk = min(dkmax,nc) 
        #
        Qs = np.zeros((nocc, nvirt, dk))
        nb = 0   # enumerate vectors that are added 
        # select new expansion vectors among the residual vectors
        for n in range(0, dkmax):
            if norms[n] > eps:
                Qs[:,:,nb] = Ainv * W[:,:,n]
                nb += 1
        assert nb==dk
        # new expansion vectors are bs + Qs
        bs_new = np.zeros((nocc,nvirt, l+dk))
        bs_new[:,:,:l] = bs
        bs_new[:,:,l:] = Qs
        # orthonormalize all expansion vectors among themselves
        bs = gram_schmidt_2d(bs_new)
        l = bs.shape[-1]   # l+dk unless linear dependent expansion vectors were removed
    else:
        raise Exception("Z-vector equations: Krylov solver for linear system did not converge after %s iterations (threshold: %e) !!!" % (it+1, conv))
    return X



###################### for Coupled-Perturbed Kohn-Sham equations #########
#
# `gram_schmidt_1d` and `Krylov_solver_CPKS` differ from `gram_schmidt_2d` and `Krylov_solver_Zvector`
# only in the dimension of the vectors: In the first case the vectors are truely 1d, while in the second case
# the vectors are 2d matrices that should be interpreted as vectors.
#

def gram_schmidt_1d(W):
    """
    orthogonalize a set of vectors vs using the Gram-Schmidt procedure.
    W holds all the vectors and has the shape (nocc,nvirt,nvectors)
    A single vector w has the shape (nocc*(nocc-1))

    Parameters:
    ===========
    W: set of linear independent vectors
    
    Returns:
    ========
    V: set of orthonormalized vectors
    """
    def olap(W):
        dim = W.shape[-1]
        # vectors W are not normalized
        norms = np.zeros(dim)
        for i in range(0, dim):
            norms[i] = np.sqrt(np.dot(W[:,i],W[:,i]))
        S = np.zeros((dim,dim))
        for i in range(0, dim):
            for j in range(0, dim):
                S[i,j] = np.dot(W[:,i],W[:,j])/(norms[i]*norms[j])
        return S
    ## check that expansion vectors are linearly independent
    S = olap(W)
    detS = la.det(S)
    if abs(detS) == 0.0:
        print "Trying to remove linear dependent expansion vectors ...",
        indx = linear_independent_cols(S, 1.0e-18)
        print "done"
        W = W[:,indx]
        S = olap(W)
#    assert la.det(S) > 0.0, "Gram-Schmidt orthogonalizer got vectors which are not linearly independent!" 
    ##
    dim = W.shape[-1]
    # V holds the orthonormalized vectors
    V = np.zeros(W.shape)
    V[:,0] = W[:,0]/la.norm(W[:,0])
    for n in range(1, dim):
        V[:,n] = np.copy(W[:,n])
        for m in range(0,n):
#            V[:,n] -= dot(V[:,m],W[:,n])*V[:,m] # numerically unstable classical Gram-Schmidt 
            V[:,n] -= np.dot(V[:,m],V[:,n])*V[:,m] # numerically stable modified Gram-Schmidt
        V[:,n] /= la.norm(V[:,n])
    # check that all vectors are orthonormal
    S = np.zeros((dim,dim))
    for i in range(0, dim):
        for j in range(0, dim):
            S[i,j] = np.dot(V[:,i],V[:,j])
    err = np.sum(abs(S - np.eye(dim)))
#    assert err < 1.0e-10, "Gram-Schmidt orthogonalization failed, error = %s" % err
    if err > 1.0e-10:
        print "WARNING: Gram-Schmidt orthogonalization failed, error = %s" % err
    return V


def Krylov_solver_CPKS(A, Adiag, B, X0, maxiter=1000, conv=1.0e-12, verbose=1):
    """
    solve the matrix equation A.X = B
    iteratively. This function is used for solving the CPKS equations. 

    Parameters:
    ===========
    A: linear operator, such that A(X) = A.X
    Adiag: diagonal elements of A-matrix, with dimension (nocc*(nocc-1))
    B: right hand side of equation, (nocc*(nocc-1), k)
    X0: initial guess vectors or None
    """
    def writeIteration(it, norms, l):
        dt = time2-time1
        print "Iteration %d: ( %d expansion vectors )   time: %2.6f s" % (it, l, dt)
        k = len(norms)
        for n in range(0, k):
            conv_str =     "not converged"
            if norms[n] < conv:
                conv_str = "converged    "
            print "  Solution vector  %3d      %s   res. norm = %e" % (n, conv_str, norms[n])
    if verbose > 0:
        print "Solve CPKS equations iteratively"
    # number of vectors
    dim,k = B.shape
    kmax = dim
    l = k
    # bs are the expansion vectors
    # initial guess
    Ainv = 1.0/Adiag # approximate inverse of A by inverse of diagonal
    if type(X0) == type(None):
        # 
        bs = np.zeros((dim,k))
        for n in range(0, k):
            bs[:,n] = Ainv * B[:,n]
    else:
        # start from solution of previous calculation
        if verbose > 0:
            print "start with solution from previous calculation"
        bs = X0
    #
    for it in range(0, maxiter):
        time1 = time.time()
        # representation of A in the basis of expansion vectors
        Ab = np.tensordot(bs, A(bs), axes=([0],[0]))
        # RHS in basis of expansion vectors
        Bb = np.tensordot(bs, B, axes=([0],[0]))
        # solve
        Xb = la.solve(Ab, Bb)
        # transform solution vector back into canonical basis
        X = np.tensordot(bs, Xb, axes=(1,0))
        # residual vectors
        W = A(X) - B
        norms = np.zeros(k)
        for n in range(0, k):
            norms[n] = la.norm(W[:,n])
        time2 = time.time()
        if verbose > 0:
            writeIteration(it, norms, l)
        # check for convergence
        if (norms < conv).all():
            if verbose > 0:
                print "CPKS equations CONVERGED (residual norms < %s) after %d iterations" \
                    % (conv, it+1)
            break
        # enlarge dimension of subspace by dk vectors
        # At most k new expansion vectors are added
        dkmax = min(kmax-l,k) # if the full space is reached, stop enlarging
        # count number of non-converged vectors
        # residual vectors that are zero cannot be used as new expansion vectors
        eps = 0.01*conv #1.0e-16
        nc = np.sum(norms > eps)
        dk = min(dkmax,nc) 
        #
        Qs = np.zeros((dim, dk))
        nb = 0   # enumerate vectors that are added 
        # select new expansion vectors among the residual vectors
        for n in range(0, dkmax):
            if norms[n] > eps:
                Qs[:,nb] = Ainv * W[:,n]
                nb += 1
        assert nb==dk
        # new expansion vectors are bs + Qs
        bs_new = np.zeros((dim, l+dk))
        bs_new[:,:l] = bs
        bs_new[:,l:] = Qs
        # orthonormalize all expansion vectors among themselves
        bs = gram_schmidt_1d(bs_new)
        l = bs.shape[-1]   # l+dk unless linear dependent expansion vectors were removed
    else:
        raise Exception("CPKS equations: Krylov solver for linear system did not converge after %s iterations (threshold: %e) !!!" % (it+1, conv))
    return X
