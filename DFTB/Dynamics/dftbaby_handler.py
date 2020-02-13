"""
This is the interface between the molecular dynamics code and the DFTBaby module, that
provides electronic excitation energies, gradients and non-adiabatic couplings.
"""
import sys
from DFTB.PES import PotentialEnergySurfaces
from DFTB import XYZ, AtomicData
from DFTB.DFTB2 import SCFNotConverged, SCFChargeInconsistency # exception
from DFTB.Solver import ExcitedStatesNotConverged, ExcitedStatesError # exceptions
from DFTB.OrbitalLocalization import OrbitalLocalization

import numpy as np
import os.path

class dftbaby_handler:
    def __init__(self,symbols,coordinates,tstep,nstates,charge, sc_threshold=0.001):
        # build list of atoms
        atomlist = []
        for s,xyz in zip(symbols, coordinates):
            Z = AtomicData.atomic_number(s)
            atomlist.append( (Z,xyz) )
        self.dt_nuc = tstep   # nuclear time step in a.u.
        self.nstates = nstates # number of electronic states including the ground state
        self.sc_threshold = sc_threshold # threshold for coefficients that are included in the
                                      # computation of the scalar coupling
        self.Nat = len(atomlist)
        self.masses = AtomicData.atomlist2masses(atomlist)
        self.pes = PotentialEnergySurfaces(atomlist, nstates, charge=charge)
        # save results from last step
        self.last_step = None
    def getBrightestState(self, coordinates):
        """
        following a call to __init__(...), this function determines which excited state has the largest
        oscillator strength. This state is taken as the initially excited state, on which
        the non-adiabatic dynamic starts.

        Returns
        -------
        state: index of the initial photoexcited state
        """
        x = np.ravel(coordinates).real
        # only excitation energies
        energies = self.pes.getEnergies(x)
        en_exc = energies[1:] - energies[0]
        # compute oscillator strengths
        tdip = self.pes.getTransitionDipoles()
        # oscillator strengths  f[I] = 2/3 omega[I] |<0|r|I>|^2
        f = 2.0/3.0 * en_exc * (tdip[0,1:,0]**2 + tdip[0,1:,1]**2 + tdip[0,1:,2]**2)
        print "The dynamics starts on the excited state with the largest oscillator strength." 
        state = 1 + np.argmax(f)
        print " State     Excitation energy / eV      Oscillator strength / e*bohr"
        for i in range(0, self.nstates-1):
            print " %d         %.7f                   %.7f" % (i+1, en_exc[i]*AtomicData.hartree_to_eV, f[i]),
            if i+1 == state:
                print "      (initial)"
            else:
                print ""
        print "selected initial state: %d" % state
        return state

    def getFragmentExcitation(self, coordinates, ifrag,iorb, afrag,aorb):
        
        x = np.ravel(coordinates).real
        # only excitation energies
        energies = self.pes.getEnergies(x)
        # We save the MO and TD-DFTB coefficients, so that the phase can be
        # aligned.
        self.pes.saveLastCalculation()
        
        # localize orbitals and project iorb->aorb excitation onto adiabatic
        # eigenstates
        LocOrb = OrbitalLocalization(self.pes.tddftb)
        proj = LocOrb.projectLocalizedExcitation(ifrag,iorb, afrag,aorb)
        
        return proj
    
    def getAll_S0(self, coordinates, has_scf=False):
        """
        ground state calculation only, in case the excited state calculation has failed.
        """
        x = np.ravel(coordinates).real
        energies = np.zeros(self.nstates)
        coupl = np.zeros((self.nstates, self.nstates))
        tdip = np.zeros((self.nstates, self.nstates, 3))
        olap = np.eye(self.nstates)
        # ground state energy
        try:
            E0, gradTot = self.pes.getEnergyAndGradient_S0(x, has_scf=has_scf)
            # The excitation energies are set to  0, this forces a hop to the ground state
            energies[:] = E0
            # convert everything into the format expected by Jens' program
            accel = - gradTot / self.masses
            accel = np.reshape(accel,(self.Nat,3))
        except SCFNotConverged as e:
            if self.last_step != None:
                print "WARNING: %s" % e
                print "Trying to continue with gradient and energy from last step!"
                energies, accel, coupl, tdip, olap = self.last_step
                self.pes.resetSCF()
            else:
                # If it fails in the first step, there is something wrong
                raise e

        return energies, accel, 0*coupl, 0*tdip, olap
        
    def getAll(self, coordinates, state):
        x = np.ravel(coordinates).real
        try:
            energies, gradTot = self.pes.getEnergiesAndGradient(x, state)
        except (SCFNotConverged, ExcitedStatesNotConverged, ExcitedStatesError) as e:
            # The SCF calculation has not converged or TDDFT has broken down
            # probably because we have hit a conical 
            # intersection or because the molecule has dissociated.
            # We try to continue on the ground state until we have passed the CI.
            print "WARNING: %s" % e
            print "Trying to continue with ground state gradient, energies=E0 and zero coupling!"
            if self.last_step != None or (state == 0):
                energies, accel, coupl, tdip, olap = self.getAll_S0(coordinates, has_scf=True)
                self.pes.resetXY()
                return energies, accel, 0*coupl, 0*tdip, olap
            else:
                # If the TD-DFT calculation fails already in the first step, there must be a serious problem
                raise e
            
        coupl, olap = self.pes.getScalarCouplings(threshold=self.sc_threshold)
        # olap = <Psi_A(t)|Psi_B(t+dt)>
        # coupl = <Psi_A|d/dR Psi_B>*dR/dt * dt
        # divide by nuclear time step
        coupl /= self.dt_nuc
        # 
        tdip  = self.pes.getTransitionDipoles()
        # convert everything into the format expected by Jens' program
        accel = - gradTot / self.masses
        accel = np.reshape(accel,(self.Nat,3))
        # save results...
        self.last_step = energies, accel, coupl, tdip, olap
        # ...and return them
        return energies, accel, coupl, tdip, olap

    #############
    # TIME SERIES:
    # To analyze trajectories it is helpful to write out additional quantities along the trajectory.
    # 
    # calculate quantities along the trajectory
    def writeTimeSeries(self, fish, time_series=[]):
        atomlist = self.pes.tddftb.dftb2.getGeometry()
        Nat = len(atomlist)  # number of QM atoms
        x = XYZ.atomlist2vector(atomlist)
        coordinates = np.reshape(x, (Nat, 3))
        symbols = [AtomicData.atom_names[Z-1] for (Z,pos) in atomlist]
        # fish is an instance of the class fish.moleculardynamics
        if "particle-hole charges" in time_series:
            ## ph charges are weighted by quantum populations |C(I)|^2
            #particle_charges, hole_charges = self.getAvgParticleHoleCharges(fish.c)
            # ph charges for coherent superposition of excited states
            particle_charges, hole_charges = self.getParticleHoleChargesSuperposition(fish.c)
            # append charges to file
            outfile = open("particle_hole_charges.xyz", "a")
            outfile.write("%d\n" % Nat)
            outfile.write(" time: %s fs\n" % (fish.time/fish.fs_to_au))
            tmp = fish.au_to_ang * coordinates
            for i in range(0, Nat):
                outfile.write("%s  %20.12f  %20.12f  %20.12f      %5.7f %5.7f\n" \
                              %(symbols[i],tmp[i][0],tmp[i][1],tmp[i][2],particle_charges[i], hole_charges[i]))
            outfile.close()
        if "particle-hole charges current" in time_series:
            # ph charges of the current electronic state
            if fish.state == 0:
                # ground state
                particle_charges = np.zeros(Nat)
                hole_charges = np.zeros(Nat)
            else:
                # excited state
                particle_charges, hole_charges = self.pes.tddftb.ParticleHoleCharges(fish.state-1)
            # append charges to file
            outfile = open("particle_hole_charges_current.xyz", "a")
            outfile.write("%d\n" % Nat)
            outfile.write(" time: %s fs\n" % (fish.time/fish.fs_to_au))
            tmp = fish.au_to_ang * coordinates
            for i in range(0, Nat):
                outfile.write("%s  %20.12f  %20.12f  %20.12f      %5.7f %5.7f\n" \
                              %(symbols[i],tmp[i][0],tmp[i][1],tmp[i][2],particle_charges[i], hole_charges[i]))
            outfile.close()
        if "Lambda2 current" in time_series:
            if fish.state == 0:
                # ground state, Lambda2 descriptor is only defined for excited states
                Lambda2 = -1.0
            else:
                # excited state
                Lambda2 = self.pes.tddftb.Lambda2[fish.state-1]

            if not os.path.isfile("lambda2_current.dat"):
                outfile = open("lambda2_current.dat", "w")
                # write header
                print>>outfile, "# TIME / fs       LAMBDA2 "
                print>>outfile, "#              0  -  charge transfer state"
                print>>outfile, "#              1  -  local excitation     "
                print>>outfile, "#             -1  -  ground state (Lambda2 undefined)"
            else:
                outfile = open("lambda2_current.dat", "a")
                
            print>>outfile, "%14.8f    %+7.4f" % (fish.time/fish.fs_to_au, Lambda2)
            outfile.close()
        if "Lambda2" in time_series:
            # Lambda2 descriptors for all excited states
            if not os.path.isfile("lambda2.dat"):
                outfile = open("lambda2.dat", "w")
                # write header
                print>>outfile, "# TIME / fs       LAMBDA2's OF EXCITED STATES"
                print>>outfile, "#              0  -  charge transfer state"
                print>>outfile, "#              1  -  local excitation     "
                print>>outfile, "#             -1  -  ground state (Lambda2 undefined)"
            else:
                outfile = open("lambda2.dat", "a")

            print>>outfile, "%14.8f    " % (fish.time/fish.fs_to_au),
            for I in range(1, self.nstates):
                print>>outfile, "  %+7.4f" % (self.pes.tddftb.Lambda2[I-1]),
            print>>outfile, ""
            
            outfile.close()
        if "transition charges current" in time_series:
            # transition charges for S0 -> current electronic state S_n
            if fish.state == 0:
                # ground state
                transition_charges = np.zeros(Nat)
            else:
                # excited state
                transition_charges = self.pes.tddftb.TransitionChargesState(fish.state-1)
            # append charges to file
            outfile = open("transition_charges_current.xyz", "a")
            outfile.write("%d\n" % Nat)
            outfile.write(" time: %s fs\n" % (fish.time/fish.fs_to_au))
            tmp = fish.au_to_ang * coordinates
            for i in range(0, Nat):
                outfile.write("%s  %20.12f  %20.12f  %20.12f      %5.7f\n" \
                              %(symbols[i],tmp[i][0],tmp[i][1],tmp[i][2], transition_charges[i]))
            outfile.close()

            
    def getAvgParticleHoleCharges(self, C):
        """
        This function computes the state-averaged particle and hole charges. 
        The density difference between and excited state and the ground state can be divided into 
        particle charges rho_p and hole charges rho_h so that
           d rho^I = rho^I - rho^0 = rho^I_p + rho^I_h
        The state averaged densities are
             ___
           d rho = sum_I |C_I|^2 * d rho^I

        and 
            ___
            rho_p/h = sum_I |C_I|^2 * d rho^I_p/h

        In tight-binding DFT The particle and hole charges are represented as spherical charge fluctuations
        around the atoms, so that it is enough to give the p-h charges on each atom
   
        Parameters:
        ===========
        C: coefficients of electronic wavefunction

        Returns:
        ========
        particle_charges, hole_charges: lists with the charges for each atom
        """
        # particle hole charges can only be calculated for QM atoms, not for
        # MM atoms. In a QM/MM calculation self.Nat is the total number of atoms,
        # but we need the number of QM atoms.
        atomlist = self.pes.tddftb.dftb2.getGeometry()
        Nat = len(atomlist)  # number of QM atoms
        particle_charges = np.zeros(Nat)
        hole_charges = np.zeros(Nat)
        for I in range(1, self.nstates):
            # I = 0 is ground state
            dqI_p, dqI_h = self.pes.tddftb.ParticleHoleCharges(I-1)
            # weight of charges is probability to be on that state
            wI = abs(C[I])**2
            particle_charges += wI * dqI_p
            hole_charges += wI * dqI_h
        # charges should sum to 0
        assert abs(np.sum(particle_charges + hole_charges)) < 1.0e-10

        return particle_charges, hole_charges

    def getParticleHoleChargesSuperposition(self, C):
        """
        This function computes the particle and hole charges for a superposition of adiabatic excited states. 

        In tight-binding DFT The particle and hole charges are represented as spherical charge fluctuations
        around the atoms, so that it is enough to give the p-h charges on each atom
   
        Parameters:
        ===========
        C: coefficients of electronic wavefunction

        Returns:
        ========
        particle_charges, hole_charges: lists with the charges for each atom
        """
        # particle hole charges can only be calculated for QM atoms, not for
        # MM atoms. In a QM/MM calculation self.Nat is the total number of atoms,
        # but we need the number of QM atoms.
        atomlist = self.pes.tddftb.dftb2.getGeometry()
        Nat = len(atomlist)  # number of QM atoms
        particle_charges = np.zeros(Nat)
        hole_charges = np.zeros(Nat)

        # exclude ground state
        coeffs = C[1:]
        # Eigenvalue solvers produce random phases, so we have to align the coefficients
        # between steps
        coeffs = self.pes.alignCoefficients(coeffs)
        
        particle_charges, hole_charges = self.pes.tddftb.ParticleHoleChargesSuperposition(coeffs)
        # charges should sum to 0
        assert abs(np.sum(particle_charges + hole_charges)) < 1.0e-10

        return particle_charges, hole_charges

