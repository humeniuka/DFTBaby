#!/usr/bin/env python
"""
The coupled-perturbed Kohn-Sham equations (CPKS) are solved for the gradients of the
MO coefficients. These gradients are not required for total energy gradients which
are computed using Furche's auxiliary functional method. They are only needed for 
gradients of properties such as the Mulliken charges.

The gradients of the MO coefficients are not necessarily continuous if there are degenerate
orbitals. Also, since the phases of eigenvectors are arbitrary, the numerical gradients of the
MO coefficients will be wrong unless the phases are matched between neighbouring geometries.
"""

import numpy as np
import scipy.linalg as la

from DFTB import Solver
from DFTB.Timer import GlobalTimer as T

@T.timer
def solve_cpks_direct(Grad, verbose=0):
    """
    The CPKS equations

         A.C = B

    are solved for C directly. The C-matrix gives the coordinates of the MO gradients
    in the basis of the MOs:

         d(orbs[a,i])/dR(k) = sum_j orbs[a,j] C[k,j,i]

    Parameters
    ----------
    Grad: instance of DFTB.ExcGradients.Gradients

    Returns
    -------
    grad_en: gradients of orbital energies, numpy array with shape (3*Nat,Norb)
           grad_en[:,i] is the gradient of energy of the i-th orbital
    grad_orbs: gradients of MO coefficients, numpy array with shape (3*Nat,Norb,Norb)
           grad_orbs[:,a,i] is the gradient of the MO coefficient orbs[a,i]
    
    """
    if verbose > 0:
        print "Solve CPKS equations directly for gradients of MO coefficients"
    # 
    # Kohn-Sham orbital energies
    en = Grad.dftb.getKSEnergies()
    # Kohn-Sham MO coefficients
    orbs = Grad.dftb.getKSCoefficients()
    S = Grad.dftb.S
    # occupation numbers, 2 - double occupied, 0 - virtual orbital
    f = Grad.dftb.getOccupation()
    # density matrix depends only on occupied orbitals
    occ_indx = np.where(f > 0.0)[0]
    # explicit gradients
    #   gradS = dS_mn / dR  with density matrix P held constant
    #   gradH = dH_mn / dR  with density matrix P held constant
    gradS = Grad.gradS
    gradH = Grad.gradH

    # nnuc - number of nuclear coordinates, length of gradient vector
    # norb - number of AOs
    nnuc, norb, norb = gradH.shape

    # Dimension of the A-matrix, since the diagonal elements C_ii are spared out
    # there are only norb*(norb-1) unknowns.
    dim = norb*(norb-1)
    # left hand side
    if verbose > 0:
        print "  build A-matrix with dimension (%d,%d)" % (dim,dim)
    A = np.zeros( (dim,dim) )

    kl = 0
    for k in range(0, norb):
        for l in range(0, norb):
            if k == l:
                continue
            # diagonal elements  A_(ij,kl) = delta_ik delta_jl (en[j] - en[i])
            A[kl,kl] = en[l]-en[k]

            if l in occ_indx:
                V_kl = np.outer(orbs[:,k], orbs[:,l])
                h_ijkl = Grad.h_func(V_kl)

                V_lk = np.outer(orbs[:,l], orbs[:,k])
                h_ijlk = Grad.h_func(V_lk)

                hh = h_ijkl + h_ijlk
                
                ij = 0
                for i in range(0, norb):
                    for j in range(0, norb):
                        if i == j:
                            continue
                        # off-diagonal elements    A_(ij,kl) -=  2 delta(l in occ) (h_ijkl + h_ijlk)
                        A[ij,kl] -= 2 * hh[i,j]
                        
                        ij += 1
                assert ij == dim
            kl += 1
    assert kl == dim
    
    if verbose > 0:
        print "  build right hand side B"
    # right hand side
    B = np.zeros( (nnuc,dim) )
    
    for c in range(0, nnuc):
        ij = 0
        for i in range(0, norb):
            for j in range(0, norb):
                if i == j:
                    continue
                
                B[c,ij]  = np.dot(orbs[:,i], np.dot(gradH[c,:,:] - en[j]*gradS[c,:,:], orbs[:,j]))

                ij += 1
        assert ij == dim
        
    for k in occ_indx:
        V_kk = np.outer(orbs[:,k], orbs[:,k])
        h_ijkk = Grad.h_func(V_kk)
        for c in range(0, nnuc):
            Ckk = -0.5 * np.dot(orbs[:,k], np.dot(gradS[c,:,:], orbs[:,k]))
            ij = 0
            for i in range(0, norb):
                for j in range(0, norb):
                    if i == j:
                        continue
                    
                    B[c,ij] += 4 * h_ijkk[i,j] * Ckk
                    
                    ij += 1
            assert ij == dim

    #
    C = np.zeros((nnuc, norb,norb))
    if verbose > 0:
        print "  solve A.c = B for each nuclear coordinate"
    # solve for all nuclear coordinates at once, if numpy is parallelized this may
    # run faster than solving for each coordinate individually.

    # a) solve equations A.c = B exactly, this may give problems if A is singular
    #Cvecs = la.solve(A, B.transpose()).transpose()
    
    # b) solve equations A.c = B in a least square sense
    x, residuals, rank, svals = la.lstsq(A, B.transpose())
    Cvecs = x.transpose()

    for c in range(0, nnuc):
        # solve A.C = B for each nuclear coordinate individually, this is slower
        #Cvec = la.solve(A, B[c,:])
        
        Cvec = Cvecs[c,:]

        ij = 0
        for i in range(0, norb):
            for j in range(0, norb):
                if i == j:
                    # diagonal entries
                    C[c,i,i] = - 0.5 * np.dot(orbs[:,i], np.dot(gradS[c,:,:], orbs[:,i]))
                    continue
                C[c,i,j] = Cvec[ij]
                ij += 1
        assert ij == dim

    if verbose > 0:
        print "  assemble gradients of MO coefficients"
    # compute gradients of MO coefficients d(orbs)/dp
    grad_orbs = np.zeros((nnuc,norb,norb))
    for c in range(0, nnuc):
        grad_orbs[c,:,:] = np.dot(orbs, C[c,:,:])

    if verbose > 0:
        print "  assemble gradients of MO energies"
    # compute gradients of energies
    grad_en = np.zeros((nnuc,norb))
    for c in range(0, nnuc):
        # compute diagonal elements B_ii
        for i in range(0, norb):
            grad_en[c,i] = np.dot(orbs[:,i], np.dot(gradH[c,:,:] - en[i]*gradS[c,:,:], orbs[:,i]))
        for k in occ_indx:
            V_kk = np.outer(orbs[:,k], orbs[:,k])
            h_iikk = np.diag( Grad.h_func(V_kk) )
            grad_en[c,:] -= 2 * h_iikk * np.dot(orbs[:,k], np.dot(gradS[c,:,:], orbs[:,k]))
    # compute sum_(k!=l) A_(ii,kl) C_kl
    for k in range(0, norb):
        for l in range(0, norb):
            if k == l:
                continue
            if l in occ_indx:

                V_kl = np.outer(orbs[:,k], orbs[:,l])
                h_iikl = np.diag( Grad.h_func(V_kl) )

                V_lk = np.outer(orbs[:,l], orbs[:,k])
                h_iilk = np.diag( Grad.h_func(V_lk) )

                A_iikl = -2*(h_iikl + h_iilk)

                grad_en -= np.outer(C[:,k,l], A_iikl)
            
    return grad_en, grad_orbs

@T.timer
def solve_cpks_iterative(Grad, verbose=0):
    """
    The CPKS equations

         A.C = B

    are solved for C iteratively using the Krylov subspace method. 
    The C-matrix gives the coordinates of the MO gradients in the basis of the MOs:

         d(orbs[a,i])/dR(k) = sum_j orbs[a,j] C[k,j,i]

    Parameters
    ----------
    Grad: instance of DFTB.ExcGradients.Gradients

    Returns
    -------
    grad_en: gradients of orbital energies, numpy array with shape (3*Nat,Norb)
           grad_en[:,i] is the gradient of energy of the i-th orbital
    grad_orbs: gradients of MO coefficients, numpy array with shape (3*Nat,Norb,Norb)
           grad_orbs[:,a,i] is the gradient of the MO coefficient orbs[a,i]
    
    """
    if verbose > 0:
        print "Solve CPKS equations iteratively for gradients of MO coefficients"
    # 
    # Kohn-Sham orbital energies
    en = Grad.dftb.getKSEnergies()
    # Kohn-Sham MO coefficients
    orbs = Grad.dftb.getKSCoefficients()
    orbsT = orbs.transpose()
    S = Grad.dftb.S
    # occupation numbers, 2 - double occupied, 0 - virtual orbital
    f = Grad.dftb.getOccupation()
    # density matrix depends only on occupied orbitals
    occ_indx = np.where(f > 0.0)[0]
    # explicit gradients
    #   gradS = dS_mn / dR  with density matrix P held constant
    #   gradH = dH_mn / dR  with density matrix P held constant
    gradS = Grad.gradS
    gradH = Grad.gradH

    # nnuc - number of nuclear coordinates, length of gradient vector
    # norb - number of AOs
    nnuc, norb, norb = gradH.shape

    # Dimension of the A-matrix, since the diagonal elements C_ii are spared out
    # there are only norb*(norb-1) unknowns.
    dim = norb*(norb-1)

    # compute omega_ij = en[j]-en[i]
    # and diagonal part of A-matrix: A_ij,ij = en[j]-en[i] - 2 * delta(j in occ) (h_ijij + h_ijji)
    omega = np.zeros(dim)
    Adiag = np.zeros(dim)
    ij = 0
    for i in range(0, norb):
        for j in range(0, norb):
            if i == j:
                continue
            omega[ij] = en[j]-en[i]

            Adiag[ij] = en[j]-en[i]
            if j in occ_indx:
                V_ij = np.outer(orbs[:,i],orbs[:,j])
                V_ji = np.outer(orbs[:,j],orbs[:,i])
                h_ijij = Grad.h_func(V_ij)[i,j]
                h_ijji = Grad.h_func(V_ji)[i,j]
                Adiag[ij] -= 2 * (h_ijij + h_ijji)
            
            ij += 1
    assert ij == dim

    
    # linear operator A.C
    def Aop(Cs):
        # number of expansion vectors
        nmax = Cs.shape[-1]
        # AC = A.C
        ACs = np.zeros(Cs.shape)
        for n in range(0, nmax):
            # flattened C
            Cvec = Cs[:,n]
            # reshape 1D vector Cvec[kl] into 2D matrix C[k,l]   with k!=l
            kl = 0
            C = np.zeros((norb,norb))
            for k in range(0, norb):
                for l in range(0, norb):
                    if k == l:
                        continue
                    
                    C[k,l] = Cvec[kl]
                    
                    kl += 1
            # 
            V_kl = np.zeros((norb,norb))
            V_lk = np.zeros((norb,norb))

            for a in range(0, norb):
                for b in range(0, norb):
                    
                    V_kl[a,b] = np.dot(orbs[a,:], np.dot(C[:,occ_indx], orbs[b,occ_indx]))
                    V_lk[a,b] = np.dot(orbs[b,:], np.dot(C[:,occ_indx], orbs[a,occ_indx]))

            h_kl = Grad.h_func(V_kl)
            h_lk = Grad.h_func(V_lk)
            # flatten hh[ij] = h_kl[i,j] + h_lk[i,j]
            hh = np.zeros(dim)
            ij = 0
            for i in range(0, norb):
                for j in range(0, norb):
                    if i == j:
                        continue

                    hh[ij] = h_kl[i,j] + h_lk[i,j]

                    ij += 1
                       
            AC = omega*Cvec - 2*hh

            ACs[:,n] = AC

        return ACs
        
    # right hand side for nuclear coordinate with index c
    def Bvec(c):
        Bc = np.zeros(dim)

        ij = 0
        for i in range(0, norb):
            for j in range(0, norb):
                if i == j:
                    continue
                
                Bc[ij]  = np.dot(orbs[:,i], np.dot(gradH[c,:,:] - en[j]*gradS[c,:,:], orbs[:,j]))
                
                for k in occ_indx:
                    V_kk = np.outer(orbs[:,k], orbs[:,k])
                    h_ijkk = Grad.h_func(V_kk)[i,j]
                    Bc[ij] -= 2 * h_ijkk * np.dot(orbs[:,k], np.dot(gradS[c,:,:], orbs[:,k]))
                    
                ij += 1
        assert ij == dim
        return Bc
    
    #
    C = np.zeros((nnuc, norb,norb))
    # solve A.C = B for each nuclear coordinate
    for c in range(0, nnuc):
        if verbose > 0:
            print "(%d of %d) CPKS for nuclear coordinate %d" % (c+1, nnuc, c+1)
        Bc = np.reshape(Bvec(c), (dim,1))
        Cvec = Solver.Krylov_solver_CPKS(Aop, Adiag, Bc, None, verbose=verbose)

        ij = 0
        for i in range(0, norb):
            for j in range(0, norb):
                if i == j:
                    # diagonal entries
                    C[c,i,i] = - 0.5 * np.dot(orbs[:,i], np.dot(gradS[c,:,:], orbs[:,i]))
                    continue
                C[c,i,j] = Cvec[ij]
                ij += 1
        assert ij == dim
        
    # compute gradients of MO coefficients d(orbs)/dp
    grad_orbs = np.zeros((nnuc,norb,norb))
    for c in range(0, nnuc):
        grad_orbs[c,:,:] = np.dot(orbs, C[c,:,:])

    # compute gradients of energies
    grad_en = np.zeros((nnuc,norb))
    for c in range(0, nnuc):
        # compute diagonal elements B_ii
        for i in range(0, norb):
            grad_en[c,i] = np.dot(orbs[:,i], np.dot(gradH[c,:,:] - en[i]*gradS[c,:,:], orbs[:,i]))
        for k in occ_indx:
            V_kk = np.outer(orbs[:,k], orbs[:,k])
            h_iikk = np.diag( Grad.h_func(V_kk) )
            grad_en[c,:] -= 2 * h_iikk * np.dot(orbs[:,k], np.dot(gradS[c,:,:], orbs[:,k]))
    # compute sum_(k!=l) A_(ii,kl) C_kl
    for k in range(0, norb):
        for l in range(0, norb):
            if k == l:
                continue
            if l in occ_indx:

                V_kl = np.outer(orbs[:,k], orbs[:,l])
                h_iikl = np.diag( Grad.h_func(V_kl) )

                V_lk = np.outer(orbs[:,l], orbs[:,k])
                h_iilk = np.diag( Grad.h_func(V_lk) )

                A_iikl = -2*(h_iikl + h_iilk)

                grad_en -= np.outer(C[:,k,l], A_iikl)
    
    return grad_en, grad_orbs


############ TESTS  #########################################################
def dftb_numerical_mo_gradients(dftb2, atomlist, dp=0.0000001):
    """
    compute gradients of orbital energies and MO coefficients numerically by finite differences

    Optional:
    ---------
    dp: step size for numerical differentiation
    """
    from DFTB import XYZ

    dftb2.setGeometry(atomlist)

    dftb2.getEnergy()
    orbe0 = dftb2.getKSEnergies()
    orbs0 = dftb2.getKSCoefficients()
    #
    atomlist0 = dftb2.getGeometry()
    vec = XYZ.atomlist2vector(atomlist0)
    Nat = len(atomlist)
    Norb,Norb = orbs0.shape
    dorbe = np.zeros((3*Nat, Norb))
    dorbs = np.zeros((3*Nat, Norb, Norb))
    for i in range(0, 3*Nat):
        vec_dp = np.copy(vec)
        vec_dp[i] += dp
        atomlist = XYZ.vector2atomlist(vec_dp, atomlist0)
        dftb2.setGeometry(atomlist)

        dftb2.getEnergy()
        orbe1 = dftb2.getKSEnergies()
        orbs1 = dftb2.getKSCoefficients()
    
        dorbe[i,:] = (orbe1-orbe0)/dp
        dorbs[i,:,:] = (orbs1-orbs0)/dp

    return dorbe, dorbs

def dftb_numerical_charge_gradients(dftb2, atomlist, dp=0.0000001):
    """
    compute gradients of Mulliken charges numerically

    Optional:
    ---------
    dp: step size for numerical differentiation
    """
    from DFTB import XYZ
    
    dftb2.setGeometry(atomlist)

    dftb2.getEnergy(density_mixer=None, scf_conv=1.0e-14, start_from_previous=0)
    charges0 = dftb2.getPartialCharges()
    #
    atomlist0 = dftb2.getGeometry()
    vec = XYZ.atomlist2vector(atomlist0)
    Nat = len(atomlist)
    # charge gradient
    dQdp = np.zeros((3*Nat, Nat))
    for i in range(0, 3*Nat):
        # vec + dp
        vec_dp = np.copy(vec)
        vec_dp[i] += dp
        atomlist = XYZ.vector2atomlist(vec_dp, atomlist0)
        dftb2.setGeometry(atomlist)

        dftb2.getEnergy(density_mixer=None, scf_conv=1.0e-14, start_from_previous=0)

        charges1 = dftb2.getPartialCharges()

        dQdp[i,:] = (charges1-charges0)/dp

    return dQdp
    

def test_dftb_eigenvector_derivative():
    """compare analytical and numerical gradients of MO coefficients and orbital energies"""
    from DFTB.XYZ import read_xyz, extract_keywords_xyz
    from DFTB.DFTB2 import DFTB2
    from DFTB.Analyse.Cube import CubeExporterEx
    from DFTB.Molden import MoldenExporter
    from DFTB import utils
    from DFTB.LR_TDDFTB import LR_TDDFTB
    from DFTB.ExcGradients import Gradients
    
    import sys
    import os

    usage  = "Usage: %s <xyz-file>\n" % sys.argv[0]
    usage += "   --help option will give more information\n"

    parser = utils.OptionParserFuncWrapper([\
       DFTB2.__init__, DFTB2.runSCC, \
       LR_TDDFTB.getEnergies, LR_TDDFTB.saveAbsorptionSpectrum, LR_TDDFTB.analyseParticleHole, \
       LR_TDDFTB.graphical_analysis, \
       CubeExporterEx.exportCubes, MoldenExporter.export, \
       Gradients.getGradients], \
                usage)

    (options,args) = parser.parse_args(DFTB2.__init__)

    if len(args) < 1:
        print usage
        exit(-1)

    xyz_file = args[0]
    atomlist = read_xyz(xyz_file)[0]
    kwds = extract_keywords_xyz(xyz_file)

    tddftb = LR_TDDFTB(atomlist, **options)

    (options,args) = parser.parse_args(tddftb.getEnergies)
    (scf_options,args) = parser.parse_args(tddftb.dftb2.runSCC)
    options.update(scf_options)
    tddftb.setGeometry(atomlist, charge=kwds.get("charge", 0.0))
    tddftb.getEnergies(**options)

    
    grad = Gradients(tddftb)
    grad.gradient(I=0, save_intermediates_CPKS=1)

    dE_an, dX_an = grad.getMOgradients()

    dftb2 = tddftb.dftb2
    E = dftb2.getKSEnergies()
    X = dftb2.getKSCoefficients()
    dE_num, dX_num = dftb_numerical_mo_gradients(dftb2, atomlist)
    print "eigenvalues"
    print E
    print "numerical eigenvalue gradients"
    print dE_num
    print "analytical eigenvalue gradients"
    print dE_an
#    print "numerical eigenvector gradients"
#    print dX_num
#    print "analytical eigenvector gradients"
#    print dX_an
    err_dE = la.norm(dE_num-dE_an)
    err_dX = la.norm(dX_num-dX_an)
    assert err_dE < 1.0e-3, "err(dE) = %s" % err_dE
    # Numerical gradients of MO coefficients are wrong because the phases (+1 or -1) of
    # the eigenvectors are arbitrary. The phases should be aligned between neighbouring
    # geometries. 
    #assert err_dX < 1.0e-3, "err(dX) = %s" % err_dX

def test_dftb_charge_derivative():
    """compare analytical and numerical gradients of Mulliken charges"""
    from DFTB.XYZ import read_xyz, extract_keywords_xyz
    from DFTB.DFTB2 import DFTB2
    from DFTB.Analyse.Cube import CubeExporterEx
    from DFTB.Molden import MoldenExporter
    from DFTB import utils
    from DFTB.LR_TDDFTB import LR_TDDFTB
    from DFTB.ExcGradients import Gradients
    
    import sys
    import os

    usage  = "Usage: %s <xyz-file>\n" % sys.argv[0]
    usage += "   --help option will give more information\n"

    parser = utils.OptionParserFuncWrapper([\
       DFTB2.__init__, DFTB2.runSCC, \
       LR_TDDFTB.getEnergies, LR_TDDFTB.saveAbsorptionSpectrum, LR_TDDFTB.analyseParticleHole, \
       LR_TDDFTB.graphical_analysis, \
       CubeExporterEx.exportCubes, MoldenExporter.export, \
       Gradients.getGradients], \
                usage)

    (options,args) = parser.parse_args(DFTB2.__init__)

    if len(args) < 1:
        print usage
        exit(-1)

    xyz_file = args[0]
    atomlist = read_xyz(xyz_file)[0]
    kwds = extract_keywords_xyz(xyz_file)

    tddftb = LR_TDDFTB(atomlist, **options)

    (options,args) = parser.parse_args(tddftb.getEnergies)
    (scf_options,args) = parser.parse_args(tddftb.dftb2.runSCC)
    options.update(scf_options)
    tddftb.setGeometry(atomlist, charge=kwds.get("charge", 0.0))
    tddftb.getEnergies(**options)

    
    grad = Gradients(tddftb)
    grad.gradient(I=0, save_intermediates_CPKS=1)

    dQdp_ana = grad.getChargeGradients()
    
    dftb2 = tddftb.dftb2
    dQdp_num = dftb_numerical_charge_gradients(dftb2, atomlist)
    print "partial Mulliken charges"
    print dftb2.getPartialCharges()
    print "numerical charge gradients"
    print dQdp_num
    print "analytical charge gradients"
    print dQdp_ana
    print "difference"
    print dQdp_num-dQdp_ana
    
    err_dQ = la.norm(dQdp_num-dQdp_ana)
    print "err(dQdp) = %s" % err_dQ
    #assert err_dQ < 1.0e-4, "err(dQdp) = %s" % err_dQ
    
    # show timings
    print T


if __name__ == "__main__":
#    test_eigenvector_derivative()
    test_dftb_charge_derivative()
#    test_dftb_eigenvector_derivative()
