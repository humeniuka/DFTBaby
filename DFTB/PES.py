"""
This is the interface through which the DFTB module should be
accessed by external programs. As an example DFTB/optimize.py shows how to 
optimize excited states with numpy's optimization routines.
"""
from DFTB.LR_TDDFTB import LR_TDDFTB
from DFTB.ExcGradients import Gradients
from DFTB.ScalarCoupling import ScalarCoupling
from DFTB import XYZ, utils
from DFTB.DFTB2 import DFTB2
from DFTB.ImplicitSolvent import SolventCavity
from DFTB.Timer import GlobalTimer as T
from DFTB import optparse

import numpy as np

class PotentialEnergySurfaces(object):
    def __init__(self, atomlist, Nst=2, **kwds):
        """
        The DFTB module has many parameters which can be set from the command
        line using options, e.g. --scf_conv=1.0e-10. 
        During the initialization 
           - the command line arguments are parsed for options
           - Slater-Koster tables and repulsive potentials
             are loaded for the atom pair present in the molecule

        Parameters:
        ===========
        atomlist: list of tuples (Zi,[xi,yi,zi]) for each atom in the molecule
        Nst: number of electronic states (including ground state)
        """
        usage = "Type --help to show all options for DFTB"

        parser = optparse.OptionParserFuncWrapper(
            [DFTB2.__init__, 
             DFTB2.runSCC,
             SolventCavity.__init__,
             LR_TDDFTB.getEnergies
            ], usage, section_headers=["DFTBaby"],
            unknown_options="ignore")
        options, args = parser.parse_args(DFTB2.__init__)
        self.atomlist = atomlist
        self.tddftb = LR_TDDFTB(atomlist, **options)
        
        solvent_options, args = parser.parse_args(SolventCavity.__init__)
        solvent_cavity = SolventCavity(**solvent_options)
        self.tddftb.dftb2.setSolventCavity(solvent_cavity)
                                           
        self.grads = Gradients(self.tddftb)
        self.tddftb.setGeometry(atomlist, charge=kwds.get("charge", 0.0))
        self.options, args = parser.parse_args(self.tddftb.getEnergies)
        self.scf_options, args = parser.parse_args(self.tddftb.dftb2.runSCC)
        self.options.update(self.scf_options)
        self.Nst = Nst
#        # always use iterative diagonalizer for lowest Nst-1 excited states
        self.options["nstates"] = Nst-1
        # save geometry, orbitals and TD-DFT coefficients from 
        # last calculation
        self.last_calculation = None
        # save transition dipoles from last calculation
        self.tdip_old = None
    def getEnergies(self, x):
        """
        This functions computes the adiabatic energies ONLY. It should
        not be used during a dynamics simulation
        """
        self.atomlist = XYZ.vector2atomlist(x, self.atomlist)
        self.tddftb.setGeometry(self.atomlist, keep_charge=True, update_symmetry=False)
        self.tddftb.getEnergies(**self.options)
        # total ground state energy
        E0 = self.tddftb.dftb2.E_tot
        # excitation energies
        Exc = self.tddftb.getExcEnergies()
        # total energies of ground and excited states
        energies = np.hstack(([E0], E0+Exc))
        return energies[:self.Nst]
    def getEnergy_S0(self, x):
        """
        This function compute the ground state energy ONLY.
        """
        self.atomlist = XYZ.vector2atomlist(x, self.atomlist)
        self.tddftb.setGeometry(self.atomlist, keep_charge=True, update_symmetry=False)
        E0 = self.tddftb.dftb2.getEnergy(**self.scf_options)
        energies = np.array([E0])  # only ground state energy
        return energies
    @T.timer
    def getEnergiesAndGradient(self, x, gradient_state):
        """
        This function computes the adiabatic electronic energies
        for the nuclear geometry contained in 'x' and the gradient
        of the total energy of state 'gradient_state'.

        Parameters:
        ===========
        x: a numpy array with the nuclear coordinates, 
           x[3*i:3*(i+1)] should contain the cartesian coordinates of atom i
        gradient_state: index of the electronic state, whose gradient should 
           be calculated. 0 stands for the ground state, 1 for the 1st excited
           state, etc.

        Returns:
        ========
        energies: numpy array with total energies (in Hartree) for all electronic
           states
        gradTot: total gradient of the selected state,
           gradTot[3*i:3*(i+1)] contains the energy derivative w/r/t the 3 cartesian
           coordinates of atom i.
        """
        self.atomlist = XYZ.vector2atomlist(x, self.atomlist)
        self.tddftb.setGeometry(self.atomlist, keep_charge=True, update_symmetry=False)
        self.tddftb.getEnergies(**self.options)
        # total ground state energy
        E0 = self.tddftb.dftb2.E_tot
        # excitation energies
        Exc = self.tddftb.getExcEnergies()
        # total energies of ground and excited states
        energies = np.hstack(([E0], E0+Exc))

        gradVrep, gradE0, gradExc = self.grads.gradient(I=gradient_state, save_intermediates_CPKS=1)
        gradTot = gradVrep + gradE0 + gradExc
        # QM/MM
        if self.tddftb.dftb2.qmmm != None:
            # expand gradient to the full system
            gradTot = self.tddftb.dftb2.qmmm.getGradientFull(gradTot)
        # confining cavity
        if self.tddftb.dftb2.cavity != None:
            if self.tddftb.dftb2.qmmm != None:
                atomlist = self.tddftb.dftb2.qmmm.getGeometryFull()
            else:
                atomlist = self.atomlist
            gradTot += self.tddftb.dftb2.cavity.getGradient(atomlist)
        return energies, gradTot
    def getEnergyAndGradient_S0(self, x, has_scf=False):
        """
        obtain the ground state energy and gradient while skipping an excited state calculation.

        This function should be called with has_scf=True directly after a failed TD-DFTB calculation. Then the scf calculation is skipped as well, assuming that the converged ground state orbitals are available.

        Returns:
        ========
        E0: ground state energy
        gradTot: gradient on ground state
        """
        # ground state energy
        if has_scf == True:
            # SCF calculation has been done already
            E0 = self.tddftb.dftb2.E_tot
        else:
            self.atomlist = XYZ.vector2atomlist(x, self.atomlist)
            self.tddftb.setGeometry(self.atomlist, keep_charge=True, update_symmetry=False)
            # SCF calculation
            E0 = self.tddftb.dftb2.getEnergy(**self.scf_options)
        # try to compute gradient on ground state
        gradVrep, gradE0, gradExc = self.grads.gradient(I=0, save_intermediates_CPKS=1)
        gradTot = gradVrep + gradE0 + gradExc
        # QM/MM
        if self.tddftb.dftb2.qmmm != None:
            # expand gradient to the full system
            gradTot = self.tddftb.dftb2.qmmm.getGradientFull(gradTot)
        # confining cavity
        if self.tddftb.dftb2.cavity != None:
            if self.tddftb.dftb2.qmmm != None:
                atomlist = self.tddftb.dftb2.qmmm.getGeometryFull()
            else:
                atomlist = self.atomlist
            gradTot += self.tddftb.dftb2.cavity.getGradient(atomlist)

        return E0, gradTot
    def resetXY(self):
        # reset excitation vectors, so that in the next iteration
        # we start with a clean guess
        self.tddftb.XmY = None
        self.tddftb.XpY = None
    def resetSCF(self):
        # set Mulliken partial charges to zero and density matrix P to
        # reference density matrix P0, so that in the next iteration
        # we start with a clean guess
        self.tddftb.dftb2.dq *= 0.0
        self.tddftb.dftb2.P = self.tddftb.dftb2.density_matrix_ref()
    def getScalarCouplings(self, threshold=0.01):
        """
        Parameters:
        ===========
        threshold: smaller excitation coefficients are neglected in the calculation
           of overlaps

        Returns:
        ========
        coupl: antisymmetric matrix with scalar non-adiabatic coupling
          coupl[A,B] ~ <Psi_A|d/dt Psi_B>*dt

        The coupling matrix has to be divided by the nuclear time step dt

        It is important that this function is called
        exactly once after calling 'getEnergiesAndGradient'
        """
        atomlist2 = self.atomlist
        orbs2 = self.tddftb.dftb2.orbs
        C2 = self.tddftb.Cij[:self.Nst-1,:,:] # only save the states of interest
        # off-diagonal elements of electronic hamiltonian
        SC = ScalarCoupling(self.tddftb)
        # retrieve last calculation
        if self.last_calculation == None:
            # This is the first calculation, so use the same data
            atomlist1, orbs1, C1 = atomlist2, orbs2, C2
        else:
            atomlist1, orbs1, C1 = self.last_calculation
        Sci = SC.ci_overlap(atomlist1, orbs1, C1, atomlist2, orbs2, C2, self.Nst, threshold=threshold)
        # align phases
        # The eigenvalue solver produces vectors with arbitrary global phases
        # (+1 or -1). The orbitals of the ground state can also change their signs.
        # Eigen states from neighbouring geometries should change continuously.
        signs = np.sign(np.diag(Sci))
        self.P = np.diag(signs)
        # align C2 with C1
        C2 = np.tensordot(self.P[1:,:][:,1:], C2, axes=(0,0))
        # align overlap matrix
        Sci = np.dot(Sci, self.P)
        # The relative signs for the overlap between the ground and excited states at different geometries
        # cannot be deduced from the diagonal elements of Sci. The phases are chosen such that the coupling
        # between S0 and S1-SN changes smoothly for most of the states.
        if hasattr(self, "last_coupl"):
            s = np.sign(self.last_coupl[0,1:]/Sci[0,1:])
            w = np.abs(Sci[0,1:] - self.last_coupl[0,1:])
            mean_sign = np.sign( np.sum(w*s)/np.sum(w) )
            sign0I = mean_sign
            for I in range(1,self.Nst):
                Sci[0,I] *= sign0I
                Sci[I,0] *= sign0I
        #
        state_labels = [("S%d" % I) for I in range(0, self.Nst)]
        if self.tddftb.dftb2.verbose > 0:
            print "Overlap <Psi(t)|Psi(t+dt)>"
            print "=========================="
            print utils.annotated_matrix(Sci, state_labels, state_labels)
        # overlap between wavefunctions at different time steps
        olap = np.copy(Sci)
        coupl = Sci
        # coupl[A,B] = <Psi_A(t)|Psi_B(t+dt)> - delta_AB
        #            ~ <Psi_A(t)|d/dR Psi_B(t)>*dR/dt dt
        # TEST
        # The scalar coupling matrix should be more or less anti-symmetric
        # provided the time-step is small enough
        # set diagonal elements of coupl to zero
        coupl[np.diag_indices_from(coupl)] = 0.0
        err = np.sum(abs(coupl+coupl.transpose()))
        if err > 1.0e-1: 
            print "WARNING: Scalar coupling matrix is not antisymmetric, error = %s" % err
        # 
        # Because of the finite time-step it will not be completely antisymmetric,
        # so antisymmetrize it
        coupl = 0.5 * (coupl - coupl.transpose())
        # save coupling for establishing relative phases
        self.last_coupl = coupl
        state_labels = [("S%d" % I) for I in range(0, self.Nst)]
        if self.tddftb.dftb2.verbose > 0:
            print "Scalar Coupling"
            print "==============="
            print utils.annotated_matrix(coupl, state_labels, state_labels)

        # save last calculation
        self.last_calculation = (atomlist2, orbs2, C2)
        return coupl, olap
    def saveLastCalculation(self):
        """
        save geometry, MO coefficients and TD-DFTB coefficients of last calculation.
        This data is needed to align the phases between neighbouring steps.
        """
        atomlist = self.atomlist
        orbs = self.tddftb.dftb2.orbs
        # only save the states of interest
        C = self.tddftb.Cij[:self.Nst-1,:,:] 
        # save last calculation
        self.last_calculation = (atomlist, orbs, C)
        
    def alignCoefficients(self, c):
        """
        align phases of coefficients of adiabatic states
        """
        if hasattr(self, "P"):
            c = np.dot(self.P[1:,:][:,1:], c)
        return  c
    
    def getTransitionDipoles(self):
        """
        computes and aligns transition dipoles between all Nst states
        This function should be called only once directly after getScalarCouplings!
        """
        tdip = self.tddftb.getTransitionDipolesExc(self.Nst)
        # eigenvectors can have arbitrary phases
        # align transition dipoles with old ones
        tdip_old = self.tdip_old
        if type(tdip_old) == type(None):
            tdip_old = tdip
        for A in range(0, self.Nst):
            for B in range(0, self.Nst):
                sign = np.dot(tdip_old[A,B,:], tdip[A,B,:])
                if sign < 0.0:
                    tdip[A,B,:] *= -1.0
        self.tdip_old = tdip
        if self.tddftb.dftb2.verbose > 0:
            state_labels = [("S%d" % I) for I in range(0, self.Nst)]
            print "Transition Dipoles"
            print "=================="
            print "X-Component"
            print utils.annotated_matrix(tdip[:,:,0], state_labels, state_labels)
            print "Y-Component"
            print utils.annotated_matrix(tdip[:,:,1], state_labels, state_labels)
            print "Z-Component"
            print utils.annotated_matrix(tdip[:,:,2], state_labels, state_labels)
        #
        return tdip

