"""
gradients of excitation energies (using Furche's variational method)
and gradients of MO coefficients (iterative method similar to CPHF equations)
"""
import numpy as np
import numpy.linalg as la

from DFTB import XYZ, utils
from DFTB.Timer import GlobalTimer as T
from DFTB.Solver import Krylov_solver_Zvector, _ApBv # for solving Z-vector equations iteratively
from DFTB import CPKS
from DFTB import GammaApproximation

from DFTB.SKMatrixElements import SLAKO_IMPLEMENTATION

from DFTB.DiskMemory import GlobalDiskMemory as Mem

@T.timer
def gradients_H0andS(atomlist, valorbs, SKT, Mproximity):
    """
    # check that F90 and Python implementation give the same results
    # parallelized fortran code
    import time
    ta = time.time()
    gradS_f90, gradH0_f90 = _gradients_H0andS_f90(atomlist, valorbs, SKT, Mproximity)
    tb = time.time()
    time_f90 = tb-ta
    # serial code
    ta = time.time()
    gradS_py, gradH0_py = _gradients_H0andS(atomlist, valorbs, SKT, Mproximity)
    tb = time.time()
    time_py = tb-ta
    # compare
#    from DFTB import utils
#    orb_labels = ["%d" % i for i in range(0, gradS_f90.shape[1])]
#    Nat = len(atomlist)
#    for i in range(0, 3*Nat):
#        print "dS/dx%d F90" % i
#        print utils.annotated_matrix(gradS_f90[i,:,:], orb_labels, orb_labels)
#        print "dS/dx%d py" % i
#        print utils.annotated_matrix(gradS_py[i,:,:], orb_labels, orb_labels)
    errS = (abs(gradS_py-gradS_f90)).max()
    assert errS < 1.0e-3, "python and Fortran implementations disagree for grad S, error = %s" % errS
#    for i in range(0, 3*Nat):
#        print "dH/dx%d F90" % i
#        print utils.annotated_matrix(gradH0_f90[i,:,:], orb_labels, orb_labels)
#        print "dH/dx%d py" % i
#        print utils.annotated_matrix(gradH0_py[i,:,:], orb_labels, orb_labels)
    errH0 = (abs(gradH0_py-gradH0_f90)).max()
    assert errH0 < 1.0e-3, "python and Fortran implementations disagree for grad H0, error = %s" % errH0

    print "Timing:"
    print "Fortran: %s seconds" % time_f90
    print "Python : %s seconds" % time_py
    exit(-1)
    """
    # 
    if SLAKO_IMPLEMENTATION == "F90":
        # Fortran code
        gradS, gradH0 = _gradients_H0andS_f90(atomlist, valorbs, SKT, Mproximity)
    else:
        # python code
        gradS, gradH0 = _gradients_H0andS(atomlist, valorbs, SKT, Mproximity)

    gradS = np.asarray(gradS)
    gradH0 = np.asarray(gradH0)
    
    return gradS, gradH0

def _gradients_H0andS(atomlist, valorbs, SKT, Mproximity):
    """
    gradients of overlap matrix S and 0-order hamiltonian matrix H0
    using Slater-Koster Rules

    Parameters:
    ===========
    atomlist: list of tuples (Zi,[xi,yi,zi]) of atom types and positions
    valorbs: list of valence orbitals with quantum numbers (ni,li,mi)
    SKT: Slater Koster table
    Mproximity: M[i,j] == 1, if the atoms i and j are close enough 
      so that the gradients for matrix elements
      between orbitals on i and j should be computed

    
    Returns:
    ========
    grad S : gradS[i,m,n] is the derivative of S_mn with respect to the i-th coordinate
    grad H0: gradH0[i,m,n] is the derivative of H0_mn with respect to the i-th coordinate
    """
    Nat = len(atomlist)
    # count valence orbitals
    Norb = 0
    for i,(Zi,posi) in enumerate(atomlist):
        Norb += len(valorbs[Zi])
    gradS = np.zeros((3*Nat, Norb, Norb))
    gradH0 = np.zeros((3*Nat, Norb, Norb))
    # iterate over atoms
    mu = 0
    for i,(Zi,posi) in enumerate(atomlist):
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbs[Zi]:
            # iterate over atoms
            nu = 0
            for j,(Zj,posj) in enumerate(atomlist):
                # iterate over orbitals on center j
                for (nj,lj,mj) in valorbs[Zj]:
                    if mu == nu:
                        pass
                    elif Mproximity[i,j] != 1:
                        pass
                    else:
                        if Zi <= Zj:
                            # the first atom given to getHamiltonian() or getOverlap()
                            # has to be always the one with lower atomic number
                            if i == j:
                                # different orbitals on the same atom
                                # no contribution
                                assert mu != nu
                                s_deriv = np.zeros(3)
                                h0_deriv = np.zeros(3)
                            else: 
                                # the hardcoded Slater-Koster rules compute the gradient 
                                # with respect to r = posj - posi
                                # but we want the gradient with respect to posi, so an additional
                                # minus sign is introduced
                                s_deriv  = -np.array(SKT[(Zi,Zj)].getOverlap(li,mi,posi, lj,mj,posj, deriv=1))
                                h0_deriv = -np.array(SKT[(Zi,Zj)].getHamiltonian0(li,mi,posi, lj,mj,posj, deriv=1))
                        else:
                            # swap atoms if Zj > Zi, since posi and posj are swapped, the gradient
                            # with respect to r = posi - posj equals the gradient with respect to
                            # posi, so no additional minus sign is needed.
                            s_deriv  = np.array(SKT[(Zj,Zi)].getOverlap(lj,mj,posj, li,mi,posi, deriv=1))
                            h0_deriv = np.array(SKT[(Zj,Zi)].getHamiltonian0(lj,mj,posj, li,mi,posi, deriv=1))
                            #
                        gradS[3*i:3*(i+1),mu,nu] += s_deriv
                        gradH0[3*i:3*(i+1),mu,nu] += h0_deriv
                        # S and H0 are hermitian/symmetric
                        gradS[3*i:3*(i+1),nu,mu] += s_deriv
                        gradH0[3*i:3*(i+1),nu,mu] += h0_deriv

                    nu += 1
            mu += 1
    return gradS, gradH0

# FASTER SLATER-KOSTER RULES FOR GRADIENTS WITH FORTRAN
from DFTB import XYZ
from DFTB.extensions import slako, grad
from DFTB.SlaterKoster.SKIntegrals import combine_slako_tables_f90
from DFTB.SKMatrixElements import atomlist2orbitals

def _gradients_H0andS_f90(atomlist, valorbs, SKT, Mproximity):
    """
    This function calls an external Fortran function that computes the matrix elements.
    Since hashes and tuples cannot be used easily in Fortran, the arguments are brought
    into a form that can be fed into the Fortran function.
    """
    atom_type_dic, spline_deg, tab_filled_SH0, tab_filled_D, \
        (S_knots, S_coefs, H_knots, H_coefs, D_knots, D_coefs) = \
                            combine_slako_tables_f90(SKT)
    
    # number of atoms
    Nat = len(atomlist)
    # orbitals
    atom_indeces,atom_types,ls,ms = atomlist2orbitals(atomlist, valorbs, atom_type_dic)
    # count valence orbitals
    Norb = len(ls)
    # distances and direction
    pos = XYZ.atomlist2vector(atomlist)
    pos = np.reshape(pos,(Nat,3)).transpose()  # pos(:,i) is the 3d position of atom i
    r,x,y,z = slako.slako.directional_cosines(pos,pos)
    # call Fortran function
    gradS,gradH0 = slako.slako.gradients_h0ands(atom_indeces,atom_types,ls,ms,
                                           r,x,y,z,
                                           Mproximity,
                                           spline_deg, tab_filled_SH0,
                                           S_knots, S_coefs, H_knots, H_coefs)
    return gradS,gradH0

#

def gradient_Vrep(atomlist, distances, directions, VREP):
    """
    compute the gradient of the repulsive potential

    Parameters:
    ===========
    atomlist: list of tuples (Zi, [xi,yi,zi]) for each atom
    distances: matrix with distances between atoms, distance[i,j]
      is the distance between atoms i and j
    directions: directions[i,j,:] is the unit vector pointing from
      atom j to atom i
    VREP: dictionary, VREP[(Zi,Zj)] has to be an instance of RepulsivePotential
      for the atom pair Zi-Zj
    """
    Nat = len(atomlist)
    grad = np.zeros(3*Nat)
    for i,(Zi,posi) in enumerate(atomlist):
        for j,(Zj,posj) in enumerate(atomlist):
            if i == j:
                continue
            if Zi > Zj:
                Z1 = Zj
                Z2 = Zi
            else:
                Z1 = Zi
                Z2 = Zj
            Rij = distances[i,j]
            eij = directions[i,j,:]
            Vij_deriv = VREP[(Z1,Z2)].getVrep(Rij, deriv=1)
            grad[i*3:i*3+3] += Vij_deriv * eij
    return grad

class Gradients:
    def __init__(self, tddftb):
        self.tddftb = tddftb
        self.dftb = tddftb.dftb2
        # save last calculations
        self.Zia = None
        self.I = None
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
        # 
        #assert nocc+nvirt == len(orbe)
        #"Analytical gradients only work if all single excitations are included, the active space cannot be restricted using the 'nr_active_occ' or 'nr_active_virt' options! The RPA eigenvalue problem should be solved for the lowest states using the Davidson-like algorithm."
        return orbe_occ, orbe_virt, orbs_occ, orbs_virt, nocc, nvirt
    def gradient(self, I=0, check_zvector=0, save_intermediates_CPKS=0):
        ###########################################
        has_active_space = self.tddftb.hasActiveSpace()
        if has_active_space == True:
            # UGLY HACK for limited active space
            # retrieve limited active space
            active_occ, active_virt = self.tddftb.getActiveOrbitals()
            # set active space to the full space 
            self.tddftb.setActiveOrbitals(nr_active_occ=None, nr_active_virt=None)
            # save Mulliken transition charges for active orbitals
            qoo_act, qvv_act, qov_act = self.tddftb.qtrans_oo, self.tddftb.qtrans_vv, self.tddftb.qtrans_ov
            # recalculate transition charges for ALL orbitals
            self.tddftb._MullikenTransitionCharges()
            # add 0s in the excitation vectors X+Y and X-Y for the inactive transitions
            occ, virt = self.tddftb.getActiveOrbitals()
            nocc, nvirt = len(occ), len(virt)
            XmYact, XpYact = self.tddftb.getXY()  # transition amplitudes between active orbitals
            nstates = XmYact.shape[0]
            XmY = np.zeros((nstates, nocc,nvirt))
            XpY = np.zeros((nstates, nocc,nvirt))
            for n in range(0, nstates):
                for i,occ in enumerate(active_occ):
                    for a,virt in enumerate(active_virt):
                        XmY[n,occ,virt-nocc] = XmYact[n,i,a]
                        XpY[n,occ,virt-nocc] = XpYact[n,i,a]        
            self.tddftb.XmY = XmY
            self.tddftb.XpY = XpY        
        # END OF HACK ############################
#        assert self.tddftb.ct_correction == 0, "Gradients for CT correction not implemented yet!"
        if self.tddftb.ct_correction != 0:
            print "WARNING: Gradient of CT correction not implemented yet!"
            print "         The contribution from the CT correction is neglected!"
        #
        if self.dftb.long_range_correction == 0:
            # use special formula
            gradVrep, gradE0, gradExc = self.gradient_nolc(I=I, check_zvector=check_zvector,
                                                           save_intermediates_CPKS=save_intermediates_CPKS)
        else:
            gradVrep, gradE0, gradExc = self.gradient_lc(I=I, check_zvector=check_zvector,
                                                         save_intermediates_CPKS=save_intermediates_CPKS)
        # REVERT UGLY HACK #######################
        if has_active_space == True:
            # go back to limited active space
            self.tddftb.setActiveOrbitals(len(active_occ), len(active_virt))
            # restore Mulliken transition charges between active orbitals
            self.tddftb.qtrans_oo, self.tddftb.qtrans_vv, self.tddftb.qtrans_ov = qoo_act, qvv_act, qov_act
            self.tddftb.XmY = XmYact
            self.tddftb.XpY = XpYact
        # END OF HACK ############################
        return gradVrep, gradE0, gradExc
    @T.timer
    @Mem.memmap
    def gradient_lc(self, I=0, check_zvector=0, save_intermediates_CPKS=0):
        """
        solve for the Lagrange multipliers Z and W
        and build gradients

        Parameters:
        ===========
        I: index of state, 0 = ground state
        check_zvector: 1 => compare iterative solution of Z-vector equations with exact one
                       0 => skip check
        save_intermediates_CPKS: save intermediate variables such as `gradS`, `gradH` that
               are required for the calculation of MO gradients using the coupled-perturbed
               Kohn-Sham equations (CPKS)

        Returns:
        ========
        gradVrep: gradient of repulsive potential
        gradE0: gradient of ground state energy
        gradExc: gradient of I-th excitation energy

        The total gradient is the sum of the 3 contributions
        """
        if self.dftb.verbose > 0:
            print "Compute gradient for state %d" % I
        # Gamma matrices
        atomlist = self.dftb.getGeometry()
        valorbs = self.dftb.valorbs
        distances, directions = self.dftb.distances, self.dftb.directions
        if self.dftb.verbose > 0:
            print "  gradients of gamma matrices"
        g0, g1, g0_AO, g1_AO = \
            self.dftb.gm.gamma_AOwise(atomlist, valorbs, distances, directions)
        g0lr, g1lr, g0lr_AO, g1lr_AO = \
            self.dftb.gm_lc.gamma_AOwise(atomlist, valorbs, distances, directions)
        
        # screening by solvation cavity
        if hasattr(self.dftb, "solvent_cavity") and self.dftb.solvent_cavity.implicit_solvent == 1:
            # We assume that a COSMO calculation has already been performed when the energy was calculated,
            # so that the surface points are available and only the gradient of gamma_solvent needs to
            # be calculated.
            gamma_solvent, grad_gamma_solvent = self.dftb.solvent_cavity.constructCOSMOgradient()
            # enlarge gamma matrices, so that the indeces refer to AOs instead of atoms
            g0_solvent_AO, g1_solvent_AO = GammaApproximation.atomwise2AOwise(gamma_solvent, grad_gamma_solvent, atomlist, valorbs)
            # screening changes the gamma-matrices
            g0_AO += g0_solvent_AO
            g1_AO += g1_solvent_AO
            
        # KS orbitals and energies
        orbe_occ, orbe_virt, orbs_occ, orbs_virt, nocc,nvirt = self.activeOrbitals()
        ei = np.diag(orbe_occ)
        ea = np.diag(orbe_virt)

        Nat, Nat = g0.shape

        if self.dftb.verbose > 0:
            print "  gradients of S and H0 from Slater-Koster rules"
        gradS, gradH0 = gradients_H0andS(atomlist, valorbs,
                                         self.dftb.SKT, self.dftb.Mproximity)
        S = self.dftb.S
        # linear operators
        Norb, Norb = g0_AO.shape
        # python implementation
        def _F(v):
            vp = v+v.transpose()
            Sv = np.sum(S*vp, axis=1)
            gSv = np.dot(g0_AO, Sv)
            gdSv = np.zeros((3*Nat,Norb))
            dgSv = np.zeros((3*Nat,Norb))
            for nc in range(0, 3*Nat):
                gdSv[nc,:] = np.dot(g0_AO,np.sum(gradS[nc,:,:]*vp, axis=1))
                dgSv[nc,:] = np.dot(g1_AO[nc,:,:], Sv)
            f = np.zeros((3*Nat,Norb,Norb))
            # a and b are AOs
            for a in range(0, Norb):
                for b in range(0, Norb):
                    f[:,a,b] = \
                        gradS[:,a,b]*(gSv[a] + gSv[b]) \
                       +S[a,b]*(dgSv[:,a]+dgSv[:,b]) \
                       +S[a,b]*(gdSv[:,a]+gdSv[:,b])
            f *= 0.25
            return f
        def _Flr(v):
            Sv = np.dot(S, v)
            SvT = np.dot(S, v.transpose())
            gv = g0lr_AO * v
            dSv24 = np.tensordot(gradS, v, axes=(2,1))
            dSv23 = np.tensordot(gradS, v, axes=(2,0))
            dgv = np.zeros((3*Nat,Norb,Norb))
            for nc in range(0,3*Nat):
                dgv[nc,:,:] = g1lr_AO[nc,:,:]*v
            #
            flr = np.zeros((3*Nat,Norb,Norb))
            # 1st term in eqn. for F_ab^lr
            tmp1 = np.tensordot(gradS, Sv, axes=(2,1))
            for nc in range(0, 3*Nat):
                flr[nc,:,:] += g0lr_AO*tmp1[nc,:,:]
            # 2nd term
            tmp2 = np.zeros((3*Nat,Norb,Norb))
            for nc in range(0, 3*Nat):
                tmp2[nc,:,:] = dSv24[nc,:,:]*g0lr_AO
            flr += np.tensordot(tmp2, S, axes=(2,1))
            # 3rd term
            flr += np.tensordot(gradS, Sv*g0lr_AO, axes=(2,1))
            # 4th term
            flr += np.tensordot(gradS, np.dot(S, gv), axes=(2,1))
            # 5th term
            tmp5 = np.swapaxes(np.tensordot(S, dSv23, axes=(1,2)), 0,1)
            for nc in range(0, 3*Nat):
                flr[nc,:,:] += g0lr_AO*tmp5[nc,:,:]
            # 6th term
            flr += np.swapaxes(np.tensordot(SvT*g0lr_AO, gradS, axes=(1,2)), 0,1)
            # 7th term
            tmp7 = np.zeros((3*Nat,Norb,Norb))
            for nc in range(0, 3*Nat):
                tmp7[nc,:,:] = dSv23[nc,:,:]*g0lr_AO
            flr += np.swapaxes(np.tensordot(S, tmp7, axes=(1,2)), 0,1)
            # 8th term
            tmp8 = np.tensordot(gradS, gv, axes=(2,0))
            flr += np.swapaxes(np.tensordot(S, tmp8, axes=(1,2)), 0,1)
            # 9th term
            tmp9 = np.dot(S, Sv.transpose())
            for nc in range(0, 3*Nat):
                flr[nc,:,:] += g1lr_AO[nc,:,:] * tmp9
            # 10th term
            tmp10 = np.zeros((3*Nat,Norb,Norb))
            for nc in range(0, 3*Nat):
                tmp10[nc,:,:] = SvT * g1lr_AO[nc,:,:]
            flr += np.tensordot(tmp10, S, axes=(2,1))
            # 11th term
            tmp11 = np.zeros((3*Nat,Norb,Norb))
            for nc in range(0, 3*Nat):
                tmp11[nc,:,:] = Sv * g1lr_AO[nc,:,:]
            flr += np.swapaxes(np.tensordot(S, tmp11, axes=(1,2)), 0,1)
            # 12th term
            tmp12 = np.tensordot(dgv, S, axes=(1,1))
            flr += np.swapaxes(np.tensordot(S, tmp12, axes=(1,1)), 0,1)
            ###
            flr *= 0.25
            return flr
        # slightly fast Fortran implementation of F(v) and Flr(v)
        def _F_f90(v):
            # Fortran implementation requires the axis with dimension 3*Nat to be the
            # last one
            f = grad.grad.fop(S,np.rollaxis(gradS,0,3), g0_AO,np.rollaxis(g1_AO,0,3), v)
            return np.rollaxis(f,2)
        def _Flr_f90(v):
            # Fortran implementation requires the axis with dimension 3*Nat to be the
            # last one
            flr = grad.grad.flrop(S,np.rollaxis(gradS,0,3), g0lr_AO,np.rollaxis(g1lr_AO,0,3), v)
            return np.rollaxis(flr,2)
        def F(v):
            """
            ### TEST grad.Fop extension
            import time
            ta = time.time()
            Fv_py = _F(v)
            tb = time.time()
            print "Time for F(v) with python:  %s seconds" % (tb-ta)
            ta = time.time()
            Fv_f90 = _F_f90(v)
            tb = time.time()
            print "Time for F(v) with Fortran: %s seconds" % (tb-ta)
            err = np.sum(abs(Fv_py-Fv_f90))
            assert err < 1.0e-10, "err = %s" % err
            ###
            """
            if SLAKO_IMPLEMENTATION == "F90":
                f = _F_f90(v)
            else:
                f = _F(v)
            return f
        def Flr(v):
            """
            ### TEST grad.Flrop extension
            import time
            ta = time.time()
            Flrv_py = _Flr(dD)
            tb = time.time()
            print "Time for Flr(v) with python:  %s seconds" % (tb-ta)
            ta = time.time()
            Flrv_f90 = _Flr_f90(dD)
            tb = time.time()
            print "Time for Flr(v) with Fortran: %s seconds" % (tb-ta)
            err = np.sum(abs(Flrv_py-Flrv_f90))
            assert err < 1.0e-10, "err = %s" % err
            ###
            """
            if SLAKO_IMPLEMENTATION == "F90":
                f = _Flr_f90(v)
            else:
                f = _Flr(v)
            return f
        
        if self.dftb.verbose > 0:
            print "  compute F(D-D0)"
        # density matrix
        D = 2*np.dot(orbs_occ, orbs_occ.transpose())
        # reference density matrix
        D0 = self.dftb.density_matrix_ref() 
        # 
        dD = D-D0
        FDmD0 = F(dD)
        if self.dftb.verbose > 0:
            print "  compute Flr(D-D0)"
        if self.dftb.lc_implementation == "old":
            # In JCP 143, 134120 (2015) the lc correction is applied to the full density matrix
            FlrDmD0 = Flr(D)
        else:
            FlrDmD0 = Flr(dD)            
            
        # GRADIENT OF GROUND STATE ENERGY
        if self.dftb.verbose > 0:
            print "  gradient of ground state energy"
        # energy weighted density matrix
        Den = 2*np.dot(orbs_occ, np.dot(ei,orbs_occ.transpose()))

        gradE0 = np.tensordot(gradH0, D, axes=([1,2],[0,1]))
        gradE0 += 0.5*np.tensordot(FDmD0, dD, axes=([1,2],[0,1]))
        if self.dftb.lc_implementation == "old":
            gradE0 -= 0.25*np.tensordot(FlrDmD0, D, axes=([1,2],[0,1]))
        else:
            gradE0 -= 0.25*np.tensordot(FlrDmD0, dD, axes=([1,2],[0,1]))
        gradE0 -= np.tensordot(gradS, Den, axes=([1,2],[0,1]))
        #
        # GRADIENT OF REPULSIVE POTENTIAL
        if self.dftb.verbose > 0:
            print "  gradient of repulsive potential"
        gradVrep = gradient_Vrep(self.dftb.getGeometry(), \
                                 distances, directions, self.dftb.VREP)
        # GRADIENT FROM DISPERSION CORRECTION
        if self.dftb.verbose > 0:
            print "  gradient of dispersion correction"
        if hasattr(self.dftb, "dispersion"):
            gradVrep += self.dftb.dispersion.getGradient(self.dftb.getGeometry(), \
                    distances, directions)

        if save_intermediates_CPKS == 1:
            # these quantities are needed for MO gradients
            self.g0_AO = g0_AO
            self.g0lr_AO = g0lr_AO
            self.gradS = gradS
            # dH/dR
            self.gradH = gradH0 + FDmD0 - 0.5 * FlrDmD0

        #
        if I == 0:
            # only ground state gradient needed
            # set gradient of excitation energy to 0
            gradExc = 0*gradE0
            return gradVrep, gradE0, gradExc
        # GRADIENT OF EXCITATION ENERGY
        if self.dftb.verbose > 0:
            print "  gradient of excitation energy"

        XmY, XpY = self.tddftb.getXY()
        Omega = self.tddftb.getExcEnergies()

        # select state I, use I-1 since we are counting from 0
        XmY = XmY[I-1,:,:]
        XpY = XpY[I-1,:,:]
        assert abs(np.sum(XmY*XpY)-1.0) < 1.0e-10
        Omega = Omega[I-1]
        
        qoo,qvv,qov = self.tddftb.getTransitionCharges()

        # vectors U, V and T
        Uab = np.tensordot(XpY, XmY, axes=(0,0)) \
             +np.tensordot(XmY, XpY, axes=(0,0))
        Uij = np.tensordot(XpY, XmY, axes=(1,1)) \
             +np.tensordot(XmY, XpY, axes=(1,1))
        Vab = np.tensordot(np.dot(ei,XpY),XpY, axes=(0,0)) \
             +np.tensordot(np.dot(ei,XmY),XmY, axes=(0,0))
        Vij = np.tensordot(np.dot(XpY,ea),XpY, axes=(1,1)) \
             +np.tensordot(np.dot(XmY,ea),XmY, axes=(1,1))
        Tab = np.tensordot(XpY, XpY, axes=(0,0)) \
             +np.tensordot(XmY, XmY, axes=(0,0))
        Tab *= 0.5
        Tij = np.tensordot(XpY, XpY, axes=(1,1)) \
             +np.tensordot(XmY, XmY, axes=(1,1))
        Tij *= 0.5

        def Hplus(q_pq,q_rs,q_pr,q_qs, q_ps, q_qr, v_rs):
            # term 1 in definition of H^+
            tmp = np.tensordot(q_rs,v_rs, axes=([1,2],[0,1]))
            tmp2 = np.dot(g0,tmp)
            hplus_pq = 4*np.tensordot(q_pq, tmp2, axes=(0,0))
            # term 2 
            tmp = np.tensordot(q_qs, v_rs, axes=(2,1))
            tmp2 = np.tensordot(g0lr, tmp, axes=(1,0))
            hplus_pq -= np.tensordot(q_pr, tmp2, axes=([0,2],[0,2]))
            # term 3
            tmp = np.tensordot(q_qr, v_rs, axes=(2,0))
            tmp2 = np.tensordot(g0lr, tmp, axes=(1,0))
            hplus_pq -= np.tensordot(q_ps, tmp2, axes=([0,2],[0,2]))
            #
            return hplus_pq
        def Hminus(q_ps, q_qr, q_pr, q_qs, v_rs):
            # term 1 in definition of H^-
            tmp = np.tensordot(q_qr, v_rs, axes=(2,0))
            tmp2 = np.tensordot(g0lr, tmp, axes=(1,0))
            hminus_pq = np.tensordot(q_ps, tmp2, axes=([0,2],[0,2]))
            # term 2
            tmp = np.tensordot(q_qs, v_rs, axes=(2,1))
            tmp2 = np.tensordot(g0lr, tmp, axes=(1,0))
            hminus_pq -= np.tensordot(q_pr, tmp2, axes=([0,2],[0,2]))
            #
            return hminus_pq
        # H^+_ij[T_ab]
        HpijTab = Hplus(qoo, qvv, qov, qov, qov, qov, Tab)
        # H^+_ij[T_ij]
        HpijTij = Hplus(qoo, qoo, qoo, qoo, qoo, qoo, Tij)
        Gij = HpijTab - HpijTij
        
        # build Q
        qvo = np.swapaxes(qov, 1,2)
        # Qij
        Qij = Omega * Uij - Vij + Gij

        # Qia
        Qia =  np.tensordot(XpY, Hplus(qvv,qov,qvo,qvv,qvv,qvo,XpY), \
                               axes=(1,1))
        Qia += np.tensordot(XmY, Hminus(qvv,qvo,qvo,qvv, XmY), \
                               axes=(1,1))
        Qia += Hplus(qov,qvv,qov,qvv,qov,qvv, Tab)
        Qia -= Hplus(qov,qoo,qoo,qvo,qoo,qvo, Tij)
        # Qai
        Qai = np.tensordot(XpY, Hplus(qoo,qov,qoo,qov,qov,qoo, XpY), \
                               axes=(0,0)) \
             +np.tensordot(XmY, Hminus(qov,qoo,qoo,qov, XmY), \
                               axes=(0,0))
        # Qab
        Qab = Omega * Uab + Vab
        
        # right hand side
        Ria = Qai.transpose()-Qia
        ##########################################
        # SOLVE Z-VECTOR EQUATION ITERATIVELY

        # build omega
        omega = np.outer(np.ones(orbe_occ.shape), orbe_virt) \
               -np.outer(orbe_occ,np.ones(orbe_virt.shape))
        def ApBv(vs):
            return _ApBv(g0, g0lr, qoo, qov, qvv, omega, vs, 1)
        if self.I == I:
            # use Z-vector from last calculation if gradient state didn't change
            Zia_guess = np.reshape(self.Zia, (nocc,nvirt,1))
        else:
            Zia_guess = None
        Zia = Krylov_solver_Zvector(ApBv, omega, np.reshape(Ria, (nocc,nvirt,1)),
                                    Zia_guess, verbose=self.dftb.verbose)
        Zia = np.reshape(Zia,(nocc,nvirt))
        # save Z-vector
        self.Zia = Zia
        # save gradient state
        self.I = I

        if check_zvector == 1:
            # compare with full solution
            ##### SOLVE Z-VECTOR EQUATION by constructing full TD-DFTB matrix ############
            #
            gq_ov = np.tensordot(g0, qov, axes=(1,0))
            gq_lr_oo = np.tensordot(g0lr, qoo, axes=(1,0))
            gq_lr_ov = np.tensordot(g0lr, qov, axes=(1,0))
            gq_lr_vv = np.tensordot(g0lr, qvv, axes=(1,0))
            # build (A+B)_(ia,jb)
            ApB = np.reshape(np.diag(np.reshape(omega, nocc*nvirt)), (nocc,nvirt,nocc,nvirt))
            ApB += 4*np.tensordot(qov, gq_ov, axes=(0,0))
            tmp = np.tensordot(qoo, gq_lr_vv, axes=(0,0))
            ApB -= np.swapaxes(tmp, 1, 2)
            tmp = np.tensordot(qov, gq_lr_ov, axes=(0,0))
            ApB -= np.swapaxes(tmp, 1, 3)

            ApB = np.reshape(ApB, (nocc*nvirt,nocc*nvirt))
            #
            err = np.sum(abs(np.dot(ApB, np.reshape(XpY,nocc*nvirt)) - Omega*np.reshape(XmY, nocc*nvirt)))
            assert err < 1.0e-5, "err = %s" % err
            #
            Ria_flat = np.reshape(Ria, nocc*nvirt)
            # solve for Z
            Z = la.solve(ApB, Ria_flat)
            Zia_full = np.reshape(Z, (nocc,nvirt))

            # compare with iterative solution
            err = np.sum(abs(Zia-Zia_full))
            assert err < 1.0e-10, "iterative solution of Z-vector equation failed, error = %s" % err
        ##########################################
        # build W
        Wij = Qij + Hplus(qoo,qov,qoo,qov,qov,qoo , Zia)
        Wij[np.diag_indices_from(Wij)] /= 2.0
        Wia = Qai.transpose() + np.dot(ei,Zia)
        Wai = Wia.transpose()
        Wab = Qab
        Wab[np.diag_indices_from(Wab)] /= 2.0
        W = np.bmat([[Wij, Wia], [Wai, Wab]])
        W = np.array(W)

        # ASSEMBLE GRADIENTS
        ##############################
        if self.dftb.verbose > 0:
            print "  assemble total gradient from the individual parts"

        # dH/dR
        gradH = gradH0 + FDmD0 - 0.5 * FlrDmD0

        # transform vectors to AO basis
        Too = np.dot(orbs_occ,np.dot(Tij,orbs_occ.transpose()))
        Tvv = np.dot(orbs_virt,np.dot(Tab,orbs_virt.transpose()))
        Zao = np.dot(orbs_occ,np.dot(Zia,orbs_virt.transpose()))
        orbs = np.hstack((orbs_occ,orbs_virt))
        Wtriu = np.triu(W)
        Wao = np.dot(orbs,np.dot(Wtriu, orbs.transpose()))
 
        XpYao = np.dot(orbs_occ,np.dot(XpY,orbs_virt.transpose()))
        XmYao = np.dot(orbs_occ,np.dot(XmY,orbs_virt.transpose()))

        gradExc = np.zeros(3*Nat)
        f = F(XpYao)
        flr_p = Flr(XpYao+XpYao.transpose())
        flr_m = -Flr(XmYao-XmYao.transpose()) # minus sign because Flr[v] actually calculates Flr[v^T]
        gradExc  = np.tensordot(gradH, Tvv-Too+Zao, axes=([1,2],[0,1]))
        gradExc -= np.tensordot(gradS, Wao, axes=([1,2],[0,1]))
        gradExc += 2*np.tensordot(XpYao, f, axes=([0,1],[1,2]))
        gradExc -= 0.5*np.tensordot(XpYao, flr_p, axes=([0,1],[1,2]))
        gradExc -= 0.5*np.tensordot(XmYao, flr_m, axes=([0,1],[1,2]))

        return gradVrep, gradE0, gradExc
    @T.timer
    def gradient_nolc(self, I=0, check_zvector=0, save_intermediates_CPKS=0):
        """
        solve for the Lagrange multipliers Z and W
        and build gradients for the special case that g0_lr = 0

        Parameters:
        ===========
        I: index of state, 0 = ground state
        check_zvector: 1 => compare iterative solution of Z-vector equations with exact one
                       0 => skip check
        save_intermediates_CPKS: save intermediate variables such as `gradS`, `gradH` that
               are required for the calculation of MO gradients using the coupled-perturbed
               Kohn-Sham equations (CPKS)

        Returns:
        ========
        gradVrep: gradient of repulsive potential
        gradE0: gradient of ground state energy
        gradExc: gradient of I-th excitation energy

        The total gradient is the sum of the 3 contributions
        """
        if self.dftb.verbose > 0:
            print "Compute gradient for state %d" % I
        # Gamma matrices
        atomlist = self.dftb.getGeometry()
        valorbs = self.dftb.valorbs
        if self.dftb.verbose > 0:
            print "  gradient of gamma matrix"
        distances, directions = self.dftb.distances, self.dftb.directions
        g0, g1, g0_AO, g1_AO = \
            self.dftb.gm.gamma_AOwise(atomlist, valorbs, distances, directions)
        # KS orbitals and energies
        orbe_occ, orbe_virt, orbs_occ, orbs_virt, nocc,nvirt = self.activeOrbitals()
        ei = np.diag(orbe_occ)
        ea = np.diag(orbe_virt)

        Nat, Nat = g0.shape
        if self.dftb.verbose > 0:
            print "  gradients of H0 and S from Slater-Koster tables"
        gradS, gradH0 = gradients_H0andS(atomlist, valorbs,
                                         self.dftb.SKT, self.dftb.Mproximity)
        S = self.dftb.S
        # linear operators
        Norb, Norb = g0_AO.shape
        def _F(v):
            vp = v+v.transpose()
            Sv = np.sum(S*vp, axis=1)
            gSv = np.dot(g0_AO, Sv)
            gdSv = np.zeros((3*Nat,Norb))
            dgSv = np.zeros((3*Nat,Norb))
            for nc in range(0, 3*Nat):
                gdSv[nc,:] = np.dot(g0_AO,np.sum(gradS[nc,:,:]*vp, axis=1))
                dgSv[nc,:] = np.dot(g1_AO[nc,:,:], Sv)
            f = np.zeros((3*Nat,Norb,Norb))
            # a and b are AOs
            for a in range(0, Norb):
                for b in range(0, Norb):
                    f[:,a,b] = \
                        gradS[:,a,b]*(gSv[a] + gSv[b]) \
                       +S[a,b]*(dgSv[:,a]+dgSv[:,b]) \
                       +S[a,b]*(gdSv[:,a]+gdSv[:,b])
            f *= 0.25
            return f
                # slightly fast Fortran implementation of F(v) and Flr(v)
        def _F_f90(v):
            # Fortran implementation requires the axis with dimension 3*Nat to be the
            # last one
            f = grad.grad.fop(S,np.rollaxis(gradS,0,3), g0_AO,np.rollaxis(g1_AO,0,3), v)
            return np.rollaxis(f,2)
        def F(v):
            """
            ### TEST grad.Fop extension
            import time
            ta = time.time()
            Fv_py = _F(v)
            tb = time.time()
            print "Time for F(v) with python:  %s seconds" % (tb-ta)
            ta = time.time()
            Fv_f90 = _F_f90(v)
            tb = time.time()
            print "Time for F(v) with Fortran: %s seconds" % (tb-ta)
            err = np.sum(abs(Fv_py-Fv_f90))
            assert err < 1.0e-10, "err = %s" % err
            ###
            """
            if SLAKO_IMPLEMENTATION == "F90":
                f = _F_f90(v)
            else:
                f = _F(v)
            return f
        # density matrix
        D = 2*np.dot(orbs_occ, orbs_occ.transpose())
        # reference density matrix
        D0 = self.dftb.density_matrix_ref() 
        # 
        dD = D-D0
        if self.dftb.verbose > 0:
            print "  computing F(D-D0)"
        FDmD0 = F(dD)

        # GRADIENT OF GROUND STATE ENERGY
        if self.dftb.verbose > 0:
            print "  gradient of ground state energy"

        # energy weighted density matrix
        Den = 2*np.dot(orbs_occ, np.dot(ei,orbs_occ.transpose()))

        gradE0 = np.tensordot(gradH0, D, axes=([1,2],[0,1]))
        gradE0 += 0.5*np.tensordot(FDmD0, dD, axes=([1,2],[0,1]))
        gradE0 -= np.tensordot(gradS, Den, axes=([1,2],[0,1]))
        #
        # GRADIENT OF REPULSIVE POTENTIAL
        gradVrep = gradient_Vrep(self.dftb.getGeometry(), \
                                 distances, directions, self.dftb.VREP)
        # GRADIENT FROM DISPERSION CORRECTION
        if hasattr(self.dftb, "dispersion"):
            gradVrep += self.dftb.dispersion.getGradient(self.dftb.getGeometry(), \
                    distances, directions)
        #

        if save_intermediates_CPKS == 1:
            # these quantities are needed for MO gradients
            self.g0_AO = g0_AO
            self.gradS = gradS
            # dH/dR
            self.gradH = gradH0 + FDmD0

        if I == 0:
            # only ground state gradient needed
            # set gradient of excitation energy to 0
            gradExc = 0*gradE0
            return gradVrep, gradE0, gradExc
        # GRADIENT OF EXCITATION ENERGY
        if self.dftb.verbose > 0:
            print "  gradient of excitation energy"
        XmY, XpY = self.tddftb.getXY()
        Omega = self.tddftb.getExcEnergies()

        # select state I, use I-1 since we are counting from 0
        XmY = XmY[I-1,:,:]
        XpY = XpY[I-1,:,:]
        assert abs(np.sum(XmY*XpY)-1.0) < 1.0e-10
        Omega = Omega[I-1]
        
        qoo,qvv,qov = self.tddftb.getTransitionCharges()

        # vectors U, V and T
        Uab = np.tensordot(XpY, XmY, axes=(0,0)) \
             +np.tensordot(XmY, XpY, axes=(0,0))
        Uij = np.tensordot(XpY, XmY, axes=(1,1)) \
             +np.tensordot(XmY, XpY, axes=(1,1))
        Vab = np.tensordot(np.dot(ei,XpY),XpY, axes=(0,0)) \
             +np.tensordot(np.dot(ei,XmY),XmY, axes=(0,0))
        Vij = np.tensordot(np.dot(XpY,ea),XpY, axes=(1,1)) \
             +np.tensordot(np.dot(XmY,ea),XmY, axes=(1,1))
        Tab = np.tensordot(XpY, XpY, axes=(0,0)) \
             +np.tensordot(XmY, XmY, axes=(0,0))
        Tab *= 0.5
        Tij = np.tensordot(XpY, XpY, axes=(1,1)) \
             +np.tensordot(XmY, XmY, axes=(1,1))
        Tij *= 0.5

        def Hplus(q_pq,q_rs,q_pr,q_qs, q_ps, q_qr, v_rs):
            # term 1 in definition of H^+
            tmp = np.tensordot(q_rs,v_rs, axes=([1,2],[0,1]))
            tmp2 = np.dot(g0,tmp)
            hplus_pq = 4*np.tensordot(q_pq, tmp2, axes=(0,0))
            #
            return hplus_pq
        # H^+_ij[T_ab]
        HpijTab = Hplus(qoo, qvv, qov, qov, qov, qov, Tab)
        # H^+_ij[T_ij]
        HpijTij = Hplus(qoo, qoo, qoo, qoo, qoo, qoo, Tij)
        Gij = HpijTab - HpijTij
        
        # build Q
        qvo = np.swapaxes(qov, 1,2)
        # Qij
        Qij = Omega * Uij - Vij + Gij

        # Qia
        Qia =  np.tensordot(XpY, Hplus(qvv,qov,qvo,qvv,qvv,qvo,XpY), \
                               axes=(1,1))
        Qia += Hplus(qov,qvv,qov,qvv,qov,qvv, Tab)
        Qia -= Hplus(qov,qoo,qoo,qvo,qoo,qvo, Tij)
        # Qai
        Qai = np.tensordot(XpY, Hplus(qoo,qov,qoo,qov,qov,qoo, XpY), \
                               axes=(0,0)) \
        # Qab
        Qab = Omega * Uab + Vab
        
        # right hand side
        Ria = Qai.transpose()-Qia
        ##########################################
        # SOLVE Z-VECTOR EQUATION ITERATIVELY
        # build omega
        omega = np.outer(np.ones(orbe_occ.shape), orbe_virt) \
               -np.outer(orbe_occ,np.ones(orbe_virt.shape))
        def ApBv(vs):
            return _ApBv(g0, None, None, qov, None, omega, vs, 0)
        if self.I == I:
            # use Z-vector from last calculation if gradient state didn't change
            Zia_guess = np.reshape(self.Zia, (nocc,nvirt,1))
        else:
            Zia_guess = None
        Zia = Krylov_solver_Zvector(ApBv, omega, np.reshape(Ria, (nocc,nvirt,1)),
                                    Zia_guess, verbose=self.dftb.verbose)
        Zia = np.reshape(Zia,(nocc,nvirt))
        # save Z-vector
        self.Zia = Zia
        # save gradient state
        self.I = I

        if check_zvector == 1:
            # compare with full solution
            ##### SOLVE Z-VECTOR EQUATION by constructing full TD-DFTB matrix ############
            #
            gq_ov = np.tensordot(g0, qov, axes=(1,0))
            # build (A+B)_(ia,jb)
            ApB = np.reshape(np.diag(np.reshape(omega, nocc*nvirt)), (nocc,nvirt,nocc,nvirt))
            ApB += 4*np.tensordot(qov, gq_ov, axes=(0,0))
            ApB = np.reshape(ApB, (nocc*nvirt,nocc*nvirt))
            #
            err = np.sum(abs(np.dot(ApB, np.reshape(XpY,nocc*nvirt)) - Omega*np.reshape(XmY, nocc*nvirt)))
            assert err < 1.0e-5, "err = %s" % err
            #
            Ria_flat = np.reshape(Ria, nocc*nvirt)
            # solve for Z
            Z = la.solve(ApB, Ria_flat)
            Zia_full = np.reshape(Z, (nocc,nvirt))

            # compare with iterative solution
            err = np.sum(abs(Zia-Zia_full))
            assert err < 1.0e-10, "iterative solution of Z-vector equation failed, error = %s" % err
        ###############################################
        # build W
        Wij = Qij + Hplus(qoo,qov,qoo,qov,qov,qoo , Zia)
        Wij[np.diag_indices_from(Wij)] /= 2.0
        Wia = Qai.transpose() + np.dot(ei,Zia)
        Wai = Wia.transpose()
        Wab = Qab
        Wab[np.diag_indices_from(Wab)] /= 2.0
        W = np.bmat([[Wij, Wia], [Wai, Wab]])
        W = np.array(W)

        # ASSEMBLE GRADIENTS
        ##############################
        # dH/dR
        gradH = gradH0 + FDmD0

        # transform vectors to AO basis
        Too = np.dot(orbs_occ,np.dot(Tij,orbs_occ.transpose()))
        Tvv = np.dot(orbs_virt,np.dot(Tab,orbs_virt.transpose()))
        Zao = np.dot(orbs_occ,np.dot(Zia,orbs_virt.transpose()))
        orbs = np.hstack((orbs_occ,orbs_virt))
        Wtriu = np.triu(W)
        Wao = np.dot(orbs,np.dot(Wtriu, orbs.transpose()))
 
        XpYao = np.dot(orbs_occ,np.dot(XpY,orbs_virt.transpose()))
        XmYao = np.dot(orbs_occ,np.dot(XmY,orbs_virt.transpose()))

        gradExc = np.zeros(3*Nat)
        f = F(XpYao)
        gradExc  = np.tensordot(gradH, Tvv-Too+Zao, axes=([1,2],[0,1]))
        gradExc -= np.tensordot(gradS, Wao, axes=([1,2],[0,1]))
        gradExc += 2*np.tensordot(XpYao, f, axes=([0,1],[1,2]))
        
        return gradVrep, gradE0, gradExc
    def getGradients(self, gradient_state=2, gradient_file=None, gradient_check=0):
        """
        compute gradients of total energy analytically

        see "Analytical Excited State Forces for the TD-DFTB Method"
        by Heringer et al., J. Comput. Chem., (2007), 28, 16. 
        and
            "Adiabatic time-dependent density functional methods for 
             excited state properties"
             by Furche and Ahlrichs, J. Chem. Phys. 117, 7433 (2002)
        and "Excited state geometry optimization by analytical energy
             gradient of long-range corrected time-dependent density
             functional theory" by Chiba, Tsuneda and Hirao
             J. Chem. Phys. 124, 144106 (2006)

        Parameters:
        ===========
        Gradients.gradient_state: Compute the gradient on this state analytically, 0 means ground state, 1 means 1st excited states etc.
        Gradients.gradient_file: save the total gradient (Vrep + E0 + ExcI) to this xyz-file. If this option is not set, the gradient will not be calculated at all!
        Gradients.gradient_check: compute gradient both numerically and analytically and compare (0 -> off, 1 -> perform check). You should increase the SCF-convergence threshold to 1.0e-14 and disable the DIIS mixer.
        """
        if gradient_file == None:
            return
        gradVrep, gradE0, gradExc = self.gradient(I=gradient_state)

        # compute gradients numerically
        atomlist = self.dftb.atomlist
        x0 = XYZ.atomlist2vector(atomlist)

        if gradient_check == 1:
            # repulsive potential
            gradVrepNum = utils.numerical_gradient(lambda x: self.dftb._energy_func(x, self.dftb.atomlist, "Enuc"), x0)
            print "REPULSIVE POTENTIAL"
            print "ANALYTICAL GRADIENT"
            print gradVrep
            print "NUMERICAL GRADIENT"
            print gradVrepNum
            err_vrep = la.norm(gradVrep - gradVrepNum)


            # ground state
            gradE0Num = utils.numerical_gradient(lambda x: self.dftb._energy_func(x, self.dftb.atomlist, "Eelec"), x0)
            print "GROUND STATE"
            print "ANALYTICAL GRADIENT"
            print gradE0

            print "NUMERICAL GRADIENT"
            print gradE0Num
            err_e0 = la.norm(gradE0 - gradE0Num)


            # excited states
            if gradient_state > 0:
                gradExcNum = utils.numerical_gradient(lambda x: self.tddftb.energy_func(x, self.dftb.atomlist, gradient_state-1), x0)
                print "STATE %d" % gradient_state
                print "ANALYTICAL GRADIENT"
                print gradExc

                print "NUMERICAL GRADIENT"
                print gradExcNum

                err_exc = la.norm(gradExc - gradExcNum)
            else:
                err_exc = 0.0

            print "ERROR REPULSIVE: %s" % err_vrep
            print "ERROR E0: %s" % err_e0
            print "ERROR EXCITATION: %s" % err_exc
            
            print "NOTE: To check the gradients you should tighten the SCF convergence setting (using --scf_conv=1.0e-14) and disable the density mixer (using --density_mixer=None). You also have to turn off QM/MM. If the space of excitations has been restricted to those between active occupied and virtual orbitals, the gradients will be wrong unless the state of interest is contained fully in the active space."
            assert err_vrep < 1.0e-4, "err(gradient vrep) = %s" % err_vrep
            assert err_e0 < 1.0e-4, "err(E0) = %s" % err_e0
            assert err_exc < 1.0e-4, "err(E_exc) = %s" % err_exc
        # total gradient
        gradEtot = gradVrep + gradE0 + gradExc
        if self.dftb.qmmm == None:
            atomlist = self.dftb.getGeometry()
        else:
            atomlist = self.dftb.qmmm.getGeometryFull()
            # expand gradient to the full system
            gradEtot = self.dftb.qmmm.getGradientFull(gradEtot)
        # confining cavity
        if self.dftb.cavity != None:
            if self.dftb.qmmm != None:
                atomlist = self.dftb.qmmm.getGeometryFull()
            else:
                atomlist = self.atomlist
            gradEtot += self.dftb.cavity.getGradient(atomlist)

        grad_atomlist = XYZ.vector2atomlist(gradEtot, atomlist)
        # save gradient to file in xyz-format
        XYZ.write_xyz(gradient_file, [grad_atomlist], title="Gradient of total energy of state %d" % gradient_state, units="au")
        print "Gradient of state %d written to file %s" % (gradient_state, gradient_file)

        return gradEtot

    ###### FOR GRADIENTS OF MO-COEFFICIENTS ########
    def _prepare_cpks(self):
        """
        The DFTB Hamiltonian depends explicitly on the nuclear geometries through the overlap matrix S
        and the 0-th order hamiltonian H0 and implicitly through the MO coefficients `orbs`. Calling any
        external parameter or nuclear coordinate `p`, the following quantities are required for computing
        the derivatives of the eigenvectors:

               gradS = dS/dp

               gradH = dH/dp = (dH0/dp + dHcoul/dp + dHx/dp)                explicit dependence
        
        and a matrix function 

               h(V)_ij = sum_(m,n) X_mi X_nj sum_(a,b) dH_mn/dP_ab V_ab

        This function does not return anything, but makes sure the following member variables are defined:

           gradS: dS/dp, numpy array with shape (3*Nat,Nao,Nao)
           gradH: dH/dp, numpy array with shape (3*Nat,Nao,Nao)
           h_func: callable that computes the matrix function h(V)

        """
        # define the matrix function h(V)
        def h_func(V):
            """
            evaluate the matrix function h(V) given the matrix V as an argument,

              h : R^(nao x nao) -> R^(nao x nao)

              h(V)_ij = sum_(m,n) X_mi X_nj sum_(a,b) dH_mn/dP_ab V_ab

            Here m,n,a,b enumerate AOs and i,j enumerate MOs. X_ai are the MO coefficients
            for MO i and P_ab is the density matrix. 
            where X_ai are the MO coefficients for orbital i

            Parameters
            ----------
            V: numpy array with shape (nao,nao)

            Returns
            -------
            h(V): numpy array with shape (nao,nao)
            """
            nao,nao = V.shape
            # MO coefficients
            orbs = self.dftb.getKSCoefficients()            
            # overlap matrix
            S = self.dftb.S
            G = self.g0_AO
            # --- Coulombic part ---
            # sum_ab dH^coul_mn / dP_ab V_ab
            VS = np.sum(V*S, axis=1)    #   VS_a = sum_b V_ab S_ab
            GVS = np.dot(G, VS)          #   GVS_m = sum_a G_ma VS_a
            # sum_ab dHcoul_mn/dP_ab V_ab = 1/2 S_mn sum_a (G_ma + G_na) sum_b V_ab S_ab
            #                             = 1/2 S_mn (GVS_m + GVS_n)
            dHcoulV = 0.5 * S * np.add.outer(GVS,GVS)
            #
            # h(V)_mn = sum_ab dH_mn / dP_ab V_ab
            hV = dHcoulV
            if self.dftb.long_range_correction == 1:
                # --- long-range exchange ---
                G_lr = self.g0lr_AO
                dHxV = np.zeros((nao,nao))
                dHxV += np.dot(G_lr * np.dot(S, V), S)
                dHxV += G_lr * np.dot(np.dot(S, V), S)
                dHxV += np.dot(np.dot(S, V*G_lr), S)
                dHxV += np.dot(S, np.dot(V,S) * G_lr)
                dHxV *= -1.0/8.0
                # add contribution from long-range exchange
                hV += dHxV
            # transform dHV into the basis of MOs
            # h(V)_ij = sum_mn X_mi X_nj h(V)_mn
            hV = np.dot(orbs.transpose(), np.dot(hV, orbs))

            return hV

        self.h_func = h_func
        
    def getMOgradients(self):
        """
        gradients of orbital energies and MO coefficients with respect to nuclear coordinates
        are calculated by solving the CPKS equations.

        Returns
        -------
        grad_en: gradients of orbital energies, numpy array with shape (3*Nat, Nao)
        grad_orbs: gradients of MO coefficients, numpy array with shape (3*Nat, Nao, Nao)
              grad_orbs[:,a,i] is the gradient of the MO coefficient orbs[a,i].
        """
        assert hasattr(self, "g0_AO"), "Calculate the ground state gradient first before calculating the MO gradients!"

        self._prepare_cpks()

        if self.dftb.cpks_solver == "direct":
            # direct solution of CPHF equations
            # `self` is passed as an argument, so that solve_cpks_# has access to all variables
            grad_en, grad_orbs = CPKS.solve_cpks_direct(self, verbose=self.dftb.verbose)
        elif self.dftb.cpks_solver == "iterative":
            # iterative solution using Krylov subspace method
            grad_en, grad_orbs = CPKS.solve_cpks_iterative(self, verbose=self.dftb.verbose)
        else:
            raise ValueError("cpks_solver should be 'direct' or 'iterative', but not '%s'!" % (self.dftb.cpks_solver))
            
        return grad_en, grad_orbs
        
    def getChargeGradients(self):
        """
        gradients of Mulliken charges. The partial charges are compute as

            q_A = sum_(mu in A) sum_nu (P-P0)_(mu,nu) S_(mu,nu)

        Returns
        -------
        dQdp: gradients of charges, numpy array with shape (3*Nat,Nat)
              dQdp[:,A] is the gradient for the charge on atom A.

        """
        grad_en, grad_orbs = self.getMOgradients()        
        # occupation numbers, 2 - doubly occupied, 0 - virtual orbital
        f = self.dftb.getOccupation()
        # density matrix depends only on occupied orbitals
        orbs = self.dftb.getKSCoefficients()
        occ_indx = np.where(f > 0.0)[0]
        CdCdp = np.tensordot(grad_orbs[:,:,occ_indx], orbs[:,occ_indx], axes=(2,1))
        # gradients of density matrix
        gradP = 2 * (CdCdp + np.swapaxes(CdCdp, 1,2))
        # overlap matrix 
        S = self.dftb.S
        # gradients of overlap matrix
        gradS = self.gradS
        
        # number of atoms
        atomlist = self.dftb.getGeometry()
        valorbs = self.dftb.valorbs
        Nat = len(atomlist)
        # charge gradients
        dQdp = np.zeros((3*Nat,Nat))

        # difference density matrix relative to reference density matrix
        PmP0 = self.dftb.P-self.dftb.P0
        # iterate over atoms
        mu = 0
        for i,(Zi,posi) in enumerate(atomlist):
            # iterate over orbitals on center i
            for (ni,li,mi) in valorbs[Zi]:
                # iterate over atoms
                nu = 0
                for j,(Zj,posj) in enumerate(atomlist):
                    # iterate over orbitals on center j
                    for (nj,lj,mj) in valorbs[Zj]:
                        dQdp[:,i] += gradP[:,mu,nu] * S[mu,nu] + PmP0[mu,nu] * gradS[:,mu,nu]
                        nu += 1
                mu += 1

        return dQdp
        
