"""
This module applies a simple correction to the diagonal elements of the TD-DFTB coupling
matrix which removes spurious low-energy charge transfer states.

see "Assessment of a simple correction for the long-range charge-transfer problem in TD-DFT"
by Neugebaueer, J.  Gritsenko, O. and Baerends, E.
J. Chem. Phys. 124, 214102 (2006)
"""
import numpy as np

from DFTB import AtomicData

class LongRangeCTCorrection_Neugebauer:
    def __init__(self, tddftb):
        self.tddftb = tddftb
        gaussOmega = tddftb.dftb2._construct_gaussian_overlap()
        occ, virt = tddftb.getActiveOrbitals()
        nocc, nvirt = len(occ), len(virt)
        qoo, qvv, qov = tddftb.getTransitionCharges()
        # compute orbitals density overlap
        # Sia = integral phi_i^2(r) phi_a^2(r) d^3r
        #     ~ sum_A sum_B q_A^(ii) gaussOmega_AB q_B^(aa)
        # with the usual DFTB approximations
        Sia = np.zeros((nocc,nvirt))
        # diagonal elements of qoo and qvv
        Nat = qoo.shape[0]
        qi = np.zeros((Nat, nocc)) # q_A^(ii)
        qa = np.zeros((Nat, nvirt)) # q_B^(aa)
        for i in range(0, nocc):
            qi[:,i] = qoo[:,i,i]
        for a in range(0, nvirt):
            qa[:,a] = qvv[:,a,a]
        # 
        self.Sia = np.dot(qi.transpose(), np.dot(gaussOmega, qa))
        assert self.Sia.shape == (nocc, nvirt)
        #
        # compute average distance between the orbital densities
        # using the Mulliken approximation
        # (Xia,Yia,Zia) = integral phi_i(r) r phi_i(r) d^3r 
        #               ~ sum_A [q_A^(ii) - q_A^(aa)] R_A  (R_A is a vector)
        atomlist = tddftb.dftb2.getGeometry()
        Xia = np.zeros((nocc,nvirt))
        Yia = np.zeros((nocc,nvirt))
        Zia = np.zeros((nocc,nvirt))
        for A,(ZA,posA) in enumerate(atomlist):
            xA,yA,zA = posA
            dqia = np.outer(qi[A,:], np.ones(nvirt)) - np.outer(np.ones(nocc), qa[A,:])
            Xia += dqia * xA
            Yia += dqia * yA
            Zia += dqia * zA
        self.Ria = np.sqrt(Xia**2 + Yia**2 + Zia**2)
        # approximation for the asymptotic limit of the correction
        # Delta ~ -orbe_a
        orbe = tddftb.dftb2.getKSEnergies()
        orbe_virt = orbe[virt]
        self.Dia = -np.outer(np.ones(nocc), orbe_virt)

    def diagonal_shifts(self, omega, Sc=0.0001, Rc=0.1, verbose=1):
        assert omega.shape == self.Ria.shape
        # compute the diagonal correction Tia
        Sia,Ria,Dia = self.Sia, self.Ria, self.Dia
        # The term exp(-(Rc/Ria)**2) is new, it avoids the correction
        # for delocalized orbitals
        Tia = np.exp(-(Sia/Sc)**2) \
            * np.exp(-(Rc/Ria)**2) \
            * (Dia - 1.0/Ria + (Dia - 1.0/Ria)**2/(2*omega))
        # correction is only added if Tia > 0
        Tia = np.clip(Tia,0.0, 1.0e10)
        if verbose > 0:
            occ, virt = self.tddftb.getActiveOrbitals()
            # show which excitations seem to have CT character
            thresh =1.0e-3
            # sort by KS energy differences
            indx = omega.argsort(axis=None, kind='mergesort')
            indi, inda = np.unravel_index(indx, omega.shape)
            print "Excitations i->a with Charge Transfer Character (Tia > %e)" % thresh
            print "======================================================================================================="
            print "  i  ->   a         Sia              Ria              Tia              w_ia / eV        w_ia+Tia / eV "
            print "-------------------------------------------------------------------------------------------------------"
            for ia in range(0, len(indi)):
                i,a = indi[ia], inda[ia]
                if Tia[i,a] > thresh:
                    print " %3d -> %3d    %13.8f    %13.8f    %13.8f    %13.8f    %13.8f" \
                % (occ[i],virt[a], Sia[i,a], Ria[i,a], Tia[i,a], \
                       omega[i,a]*AtomicData.hartree_to_eV, \
                       (omega[i,a]+Tia[i,a])*AtomicData.hartree_to_eV)
                # only show lowest energy excitations
                if omega[i,a] > (40.0/AtomicData.hartree_to_eV):
                    # stop at 40 eV
                    break
            print ""
        return Tia

