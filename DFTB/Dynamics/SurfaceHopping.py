#!/usr/bin/env python
"""
non-adiabatic molecular dynamics with surface hopping trajectories

run 
   SurfaceHopping.py --help
to see all options.
"""
import os
import time
import random

import numpy as np
import numpy.linalg as la
import scipy.linalg as sla

from DFTB import XYZ, AtomicData
from DFTB import optparse
from DFTB.Dynamics import  dftbaby_handler
from DFTB.MetaDynamics import meta
# We need to know the arguments of the __init__ functions for
# the following classes and expose them to the user as command line options.
from DFTB.DFTB2 import DFTB2
from DFTB.ImplicitSolvent import SolventCavity
from DFTB.LR_TDDFTB import LR_TDDFTB

class MolecularDynamics:
	def __init__(self,
                     charge=0, initial_state='0', nstates=0, 
                     nstep=1000, nuclear_step=0.1,
                     dyn_mode="E", temp=300.0, timecoupling=1.0, 
                     scalar_coupling_threshold=0.01,
                     switch_to_groundstate=1,
                     artificial_energy_conservation=0,
                     time_series=None,
                     output_step=1,
                     fragment_excitation=None,
                     decoherence_correction=0):
                """
                Parameters
                ==========
                Molecular Dynamics.charge: total charge of the molecule, the cation or anion has to be a closed shell species for TD-DFTB to work properly.
                Molecular Dynamics.nstates: number of excited states. Only the lowest states are calculated with TD-DFTB. For dynamics in the ground state `nstates` should be set to 0 to avoid the expensive calculation of excitation energies and non-adiabatic couplings.
                Molecular Dynamics.initial_state: initial electronic state of the trajectory, 0 for ground state. 'brightest' selects the state with the largest oscillator strength as the initial state. 'fragment' creates a superposition of adiabatic states corresponding to a localized excitation on a single fragment or between fragments (see option 'fragment_excitation').
                Molecular Dynamics.nstep: number of nuclear time steps.
                Molecular Dynamics.nuclear_step: length of nuclear time step for integration of Newton's equations (in fs).
                Molecular Dynamics.dyn_mode: 'T' for constant temperature, 'E' for constant energy. To equilibrate a trajectory on the ground state use 'T', non-adiabatic dynamics on excited states should be run at constant energy.
                Molecular Dynamics.temp: temperature in Kelvin, only needed if the dynamics is run at constant temperature. The temperature is controlled using a Berendsen thermostat.
                Molecular Dynamics.timecoupling: Time constant for Berendsen thermostat in fs. The strength of the coupling to the external heat bath is proportional to 1/timecoupling. 
                Molecular Dynamics.scalar_coupling_threshold: Excitation coefficients that are smaller than this threshold are neglected when calculating scalar couplings from the overlap between electronic states at successive time steps. For large molecules this value should be reduced.
                Molecular Dynamics.switch_to_groundstate: If set to 1, a hop to the ground state is forced if the S0-S1 energy gap drops below 0.1 eV. In TD-DFT(B) conical intersections to the ground state are not described correctly. If a point of point of degeneracy between S0 and S1 is reached, the TD-DFT(B) calculation usually breaks down. 
If something goes wrong in the excited state calculation, the trajectory continues with the ground state gradients until the excited state calculation starts working again.
                Molecular Dynamics.artificial_energy_conservation: Energy conservation can be enforced artificially by rescaling the velocities after each time step. This avoids the drift of the total energy observed for long dynamics simulations. This option only takes effect if dyn_mode=="E".
                Molecular Dynamics.time_series: To better analyze trajectories it is helpful to write out time series of additional quantities along the trajectory. You can specify a list of the desired quantities (so far only --time_series="['particle-hole charges']" is available).
                Molecular Dynamics.output_step: Output is written to disk only for every N'th time step. 
                Molecular Dynamics.fragment_excitation: The initial state does not have to be an adiabatic state. If initial_state='fragment', this option allows to specify a local excitation on a fragment or a charge-transfer excitation from one fragment to the other as the initial state. To this end the disconnected fragments are identified and the molecular orbitals are localized onto the fragments by the Pipek-Mezey method. The excitation is specified by a tuple of four integers '(ifrag,iorb, afrag,aorb)' with the following meaning, ifrag - index of fragment for occupied orbital (1-based), iorb - occupied orbital counted from the HOMO downward (so iorb=0 corresponds to HOMO, iorb=1 to HOMO-1 etc.), afrag - index of fragment for virtual orbital (1-based), aorb - virtual orbital counted from the LUMO upward (so aorb=0 corresponds to LUMO, aorb=1 to LUMO+1 etc.). For instance, if there are two fragments, the tuple '(1,0,1,0)' constitutes a local HOMO(1)->LUMO(1) excitation, while '(1,0,2,0)' constitutes a charge-transfer HOMO(1)->LUMO(2) excitation.
                Molecular Dynamics.decoherence_correction: If set to 1, the decoherence correction according to eqn. (17) in JCP 126, 134114 (2007) is turned on.
                """
		# units checked 06/08/2014
		self.fs_to_au=1.0/AtomicData.autime2fs # 41.34137333
		self.au_to_fs=1.0/self.fs_to_au
		self.ang_to_au=1.0/AtomicData.bohr_to_angs #1.0/0.52917721092
		self.au_to_ang=AtomicData.bohr_to_angs #0.52917721092
		self.kB=AtomicData.kBoltzmann  #3.166811382e-6
		random.seed()

		self.charge=charge
                self.nstep=nstep # number of time steps for dynamics
                self.tstep=nuclear_step*self.fs_to_au # time step for dynamics
                self.mode=dyn_mode
                self.temperature = temp

                self.decoherence_correction = decoherence_correction
                # use the recommended value for C in eqn. (17) of JCP 126, 134114 (2007)
                self.decoherence_constant = 0.1 # in Hartree
                
                # write information on coefficients and hopping:
                # 0: only write state.dat,
                # 1: also |c_i|^2 in coeff_$i.dat,
                # 2: also hopping probabilities in prob.dat and rejected hops in rej_hop.dat,
                # 3: also real and imaginary parts of coeffs, 4: coherences instead of real and imaginary parts
                self.printcoeff=2

                # write nonadiabatic couplings and transition dipole moments if available
                self.printcoup=1

		self.nstates=nstates+1   # internally nstates includes the ground state
                
		self.writeflag="xyz"

                self.scalar_coupling_threshold = scalar_coupling_threshold
                self.switch_to_groundstate = switch_to_groundstate
                self.artificial_energy_conservation = artificial_energy_conservation
                if time_series == None:
                        time_series = []
                self.time_series = time_series
                assert output_step > 0
                self.output_step = output_step
                
                self.time=0.0
                try:
                        # the initial state is given as an integer index
                        self.state = int(initial_state)
                except ValueError:
                        assert initial_state in ["brightest", "fragment"], "Allowed values for initial_state are integers or 'brightest' or 'fragment'"
                        self.state = initial_state
                        try:
                                ifrag,iorb, afrag,aorb = fragment_excitation
                                self.fragment_excitation = ifrag,iorb, afrag,aorb
                        except (ValueError,TypeError) as e:
                                print "ERROR: expected tuple of integers (ifrag,iorb, afrag,aorb) for option 'fragment_excitation'"
                                print "       but got '%s'." % fragment_excitation
                                raise e
                        
                if (self.state==0) and (self.nstates==1):
                        # DFTBaby program needs to compute excited states at least once
                        # to initialize all variable correctly that are needed for the gradient
                        # calculations. In the first step a single excited state is calculate, in
                        # later steps the excited state calculation is skipped if grounddyn==1.
                        self.nstates += 1
                        self.grounddyn=1
                else:
                        self.grounddyn=0


		if self.mode=="T" or self.mode=="M":
			self.currtemperature=200.0
			self.timecoupling=timecoupling*self.fs_to_au
                        
                if os.path.exists("dynamics.in0"):
                        print ""
                        print "'dynamics.in0' found => dynamics is RESTARTED from last step"
                        print ""
                        print "  Initial coordinates and velocities are read from   'dynamics.in0'"
                        print "  Electronic coefficients are read from              'coefficients.dat0'"
                        print "  Current state is read from                         'state.dat'"
                        print "  Non-adiabatic couplings cannot be recovered fully."
                        print ""
                        print "To avoid a restart you have to delete 'dynamics.in0'."
                        print ""
                        self.restart = 1
                else:
                        self.restart = 0

	###########################################################################################################################
	# MOLECULAR DYNAMICS ROUTINES
	###########################################################################################################################	

	def getCenterOfMass(self,coordi):
		X=0.0
		Y=0.0
		Z=0.0
		for i in range(len(coordi)):
			X=X+self.masses[i]*coordi[i][0]/self.totalmass
			Y=Y+self.masses[i]*coordi[i][1]/self.totalmass
			Z=Z+self.masses[i]*coordi[i][2]/self.totalmass
		return np.array([float(X),float(Y),float(Z)])

	def shiftToCenterOfMass(self,coordi):
		C=self.getCenterOfMass(coordi)
		self.centerOfMass=1.0*C
		return coordi-C

	def momentum(self,vel,mass):
		a=np.transpose(vel)*mass
		return np.array([sum(a[0]),sum(a[1]),sum(a[2])])

	def eliminateTranslation(self,vel,mom):
		return np.transpose(np.array([np.transpose(vel)[0]-mom[0]/self.totalmass,np.transpose(vel)[1]-mom[1]/self.totalmass,np.transpose(vel)[2]-mom[2]/self.totalmass]))

	def myCrossProduct(self,a,b):
		return np.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])

	def angularMomentum(self,coordi,vel,mass):
		L=np.array([0.0,0.0,0.0])
		for i in range(len(mass)):
			L=L+mass[i]*self.myCrossProduct(coordi[i],vel[i])
		return L	

	def momentOfInertia(self,coordi,mass):
		I=np.zeros((3,3))
		Emat=np.identity(3)
		for i in range(len(mass)):
			I=I+mass[i]*(np.inner(coordi[i],coordi[i])*Emat-np.outer(coordi[i],coordi[i]))
		return I
				
	def getAngularVelocity(self,Angmom,Inertia):
		return np.dot(la.inv(Inertia),(Angmom))

	def eliminateRotation(self,vel,omega,coordi):
		for i in range(len(vel)):
			vel[i]=vel[i]-self.myCrossProduct(omega,coordi[i])
		return vel

	def getKineticEnergy(self,vel,mass):
		kin=0.000000000000000
		for i in range(len(vel)):
			kin=kin+mass[i]*0.500000000000000*np.sum(vel[i]**2)
		return kin
	
	###### Velocity Verlet routines ######

	def getCoordVerlet(self):
		self.coordinates=self.coordinates+self.tstep*self.velocities+0.500000000000000*((self.tstep)**2)*self.accelerations

	def getVelVerlet(self):
		self.velocities=self.velocities+self.tstep*0.500000000000000*(self.old_accelerations+self.accelerations)

	#########################################################################################
	# INTERFACE TO DFTB PROGRAM
	#########################################################################################

	def initializeInterface(self):                
                assert self.nstates > 1, "You have to include at least 1 excited state."
                self.handler=dftbaby_handler.dftbaby_handler(self.symbols,self.coordinates, self.tstep, self.nstates, self.charge, sc_threshold=self.scalar_coupling_threshold)
                if self.restart == 0:
                        # If a trajectory is restarted the electronic coefficients have already been read from file.
                        self.c=np.array([0.0+0.0j for i in range(self.nstates)])
                        if self.state == "brightest":
                                self.state = self.handler.getBrightestState(self.coordinates)
                                self.c[self.state]=1.0+0.0j
                        elif self.state == "fragment":
                                ifrag,iorb, afrag,aorb = self.fragment_excitation
                                proj = self.handler.getFragmentExcitation(self.coordinates, ifrag,iorb, afrag,aorb)
                                self.c = proj[:self.nstates]
                                nrm = la.norm(self.c)
                                if nrm < 0.8:
                                        msg =  "\n\n ERROR: Norm of projection of local excitation onto adiabatic state is too small (%4.3f < 0.8)!\n" % nrm
                                        msg += "        Increase the active space and/or the number of excited states.\n"
                                        raise RuntimeError(msg)
                                # normalize coefficients to 1
                                self.c /= la.norm(self.c)
                                # select initial state randomly with probability |C_i|^2
                                prob = abs(self.c)**2
                                r = np.random.rand(1)[0]
                                probsum = np.cumsum(prob)
                                for i in range(0, self.nstates):
                                        if r <= probsum[i]:
                                                self.state = i
                                                break
                                #print "coefficients of electronic wavefunction"
                                #print " C = %s" % self.c
                                print ""
                                print "initial state (randomly selected) = %d" % self.state
                        
                        else:
                                self.c[self.state]=1.0+0.0j
                self.energy,self.accelerations,self.nonad_scalar,self.dipole,self.olap = self.handler.getAll(self.coordinates, self.state)

	def getData(self):
                self.dipole_old=1.0*self.dipole
                self.nonad_scalar_old=1.0*self.nonad_scalar
                if (self.grounddyn == 1) and (self.time > 0.0):
                        self.energy,self.accelerations,self.nonad_scalar,self.dipole,self.olap = self.handler.getAll_S0(self.coordinates)
                else:
                        self.energy,self.accelerations,self.nonad_scalar,self.dipole,self.olap = self.handler.getAll(self.coordinates, self.state)

        ############################################################################################################
        # HOPPING ROUTINES
        ############################################################################################################

        def getLocalDiabatization(self):
                """
                The coefficients of the electronic wavefunction are propagated in the local diabatic basis as
                explained in 
                   [1]  JCP 114, 10608 (2001) and 
                   [2]  JCP 137, 22A514 (2012)
                """
                S = self.olap  # overlap matrix <Psi_A(t)|Psi_B(t+dt)>
                # Loewding orthogonalization of the S matrix
                # see eqns. (B5) and (B6) in [2]
                StS = np.dot(S.transpose(), S)
                # eigenvalues of S^T.S
                L,O = la.eigh(StS)
                Lm12 = np.sqrt(1.0/L)
                # unitary transformation matrix, see eqn. (B5) in [1]
                T = np.dot(S, np.dot(O, np.dot(np.diag(Lm12), O.transpose())))
                Tinv = T.transpose()
                #assert np.sum(abs( np.dot(T, Tinv) - np.eye(T.shape[0]) )) < 1.0e-10, "T.Tinv != Id"
                #
                c0 = np.copy(self.c)  # electronic coefficients c(t)
                # 
                E0 = self.energy_last   # adiabatic energies at the beginning of the time step, E(t)
                E1 = self.energy        # adiabatic energies at the end of the time step, E(t+dt)
                H = np.dot(T, np.dot(np.diag(E1), Tinv)) # diabatic hamiltonian H(t+dt)
                Hinterp = (np.diag(E0) + H)/2.0
                # subtract lowest energy from diagonal
                Hinterp[np.diag_indices_from(Hinterp)] -= Hinterp[0,0]
                # propagator in diabatic basis, see eqn. (11) in [1]
                U = sla.expm(-1.0j * Hinterp * self.tstep)
                # at the beginning of the time step the adiabatic and diabatic basis is assumed to coincide
                c1 = np.dot(Tinv, np.dot(U, c0))  # new electronic coefficients c(t+dt) in the adiabatic basis
                self.c = c1
                # norm of electronic wavefunction
                nrmC = np.sum(abs(self.c)**2)
                #print "norm = %s" % nrmC
                assert abs(nrmC - 1.0) < 1.0e-3, "norm of electronic coefficients not conserved! norm = %s" % nrmC

                # save the diabatic hamiltonian along the trajectory
                if not hasattr(self, "Ttot_last"):
                        self.Ttot_last = np.eye(T.shape[0])
                # the transformations are concatenated to obtain the diabatic Hamiltonian relative to the
                # first time step
                Ttot = np.dot(T, self.Ttot_last)
                Ttot_inv = Ttot.transpose()
                self.Hdiab = np.dot(Ttot, np.dot(np.diag(E1), Ttot_inv))
                
                # save last transformation matrix
                self.Ttot_last = Ttot

	def getNewState(self):
		occupations=np.zeros(self.nstates)
		derivatives=np.zeros(self.nstates)
		for i in range(self.nstates):
			occupations[i]=abs(self.c[i])**2
			derivatives[i]=(occupations[i]-abs(self.c_old[i])**2)/self.tstep
                if derivatives[self.state]<0.0:
                        hopping_probabilities={}
                        prob_k=0.0
                        states_to_hop=[]
                        for k in range(self.nstates):
                                if derivatives[k]>0.0:
                                        states_to_hop.append(k)
                                        prob_k=prob_k+derivatives[k]
                        for state in states_to_hop:
                                hopping_probabilities[state]=-(derivatives[self.state]/(occupations[self.state]))*derivatives[state]*self.tstep/prob_k
                        random_number=random.random()
                        sum=0.0
                        for item in hopping_probabilities:
                                sum=sum+hopping_probabilities[item]
                                if random_number<sum:
                                        self.state=item
					if self.printcoeff>1:
						f=open("prob.dat","a")
						f.write(str(self.time*self.au_to_fs)+" "+str(hopping_probabilities[item])+" "+str(random_number)+"\n")
						f.close()
                                        break
                # added by A.Humeniuk
                # If the energy gap between the first excited state and the ground state
                # approaches zero, because the trajectory has hit a conical intersection to
                # the ground state, TD-DFT will break down. In this case, a transition
                # to the ground state is forced.
                threshold = 0.1/27.211  # 0.1 eV
                if self.state > 0 and self.switch_to_groundstate == 1:
                        gap = self.energy[self.state]-self.energy[0]
                        #print "gap = %s" % gap
                        if gap < threshold:
                                print "Conical intersection to ground state reached."
                                print "The trajectory will continue on the ground state."
                                self.state = 0
                                # trajectory will remain on ground state forever, excited states and couplings won't be calculated anymore
                                self.grounddyn = 1
                #
	def normalizeCoeff(self):
		norm=0.0
		for i in range(len(self.c)):
			norm+=abs(self.c[i])**2
		self.c=self.c/np.sqrt(norm)


	##########################################################################################################

        def getUniformlyRescaledVelocities(self):
                #hop is rejected when kinetic energy is too low
                if (self.state>self.oldstate) and ((self.energy[self.state]-self.energy[self.oldstate])>self.kinetic_energy):
			if self.printcoeff>1:
				out=open("rej_hop.dat","a")
                        	out.write(str(self.time*self.au_to_fs)+" "+str(self.oldstate)+" "+str(self.state)+"\n")
                        	out.close()
                        self.state=self.oldstate
		else:
	                velscale=np.sqrt((self.kinetic_energy+(self.energy[self.oldstate]-self.energy[self.state]))/self.kinetic_energy)
        	        self.velocities=velscale*self.velocities

	def getDecoherenceE(self):
                """
                decoherence correction according to eqn. (17) in

                G. Granucci, M. Persico,
                "Critical appraisal of the fewest switches algorithm for surface hopping",
                J. Chem. Phys. 126, 134114 (2007)

                If the trajectory is in the current state K, the coefficients of the other 
                states J != K are made to decay exponentially, C'_J = exp(-dt/tau_JK) C_J.
                The decay time is proportional to the inverse of the energy gap |E_J-E_K|,
                so that the coherences C_J*C_K decay very quickly if the energy gap between
                the two states is large. The electronic transitions become irreversible.
                """
                sm=0.0+0.0j
                for i in range(self.nstates):
                    if i<>self.state:
                        tauij=1.0/abs(self.energy[i]-self.energy[self.state])*(1.0+self.decoherence_constant/self.kinetic_energy)
                        self.c[i]*=np.exp(-self.tstep/tauij)
                        sm+=abs(self.c[i])**2
                self.c[self.state]*=np.sqrt((1.0-sm)/abs(self.c[self.state])**2)
                        
	########################################################################################################
	# MAIN TRAJECTORY ROUTINES
	########################################################################################################

	def eliminateTransRot(self):
		self.velocities=self.eliminateTranslation(self.velocities,self.momentum(self.velocities,self.masses))
                angmom=self.angularMomentum(self.coordinates,self.velocities,self.masses)
                inertia=self.momentOfInertia(self.coordinates,self.masses)
                angvel=self.getAngularVelocity(angmom,inertia)
                self.velocities=self.eliminateRotation(self.velocities,angvel,self.coordinates)

	def scaleVelocitiesT(self):
                """
                velocities are rescaled using the Berendsen thermostat (JCP (1984) 81, 3684)
                """
		self.currtemperature=(23209./float(3*self.nat-6))*self.kinetic_energy*27.2114
                # equation (11) in Berendsen's article
                scalingfactor=np.sqrt(1.+(self.tstep/self.timecoupling)*(self.temperature/self.currtemperature-1.))
                self.velocities=scalingfactor*self.velocities
        # added by A.Humeniuk
        def scaleVelocitiesEconst(self):
                """
                The velocities are rescaled so that energy conservation between two time-steps
                is fulfilled exactly. 
                """
                if self.state != self.oldstate:
                        # surface hop, velocities are rescaled somewhere else
                        return
                Ekin0 = self.kinetic_energy_old
                Epot0 = self.potential_energy_old
                Ekin1 = self.kinetic_energy
                Epot1 = self.energy[self.state]
                # If the integration step is small enough an the analytical gradients
                # are correct, the energy change should be zero
                scalingfactor = np.sqrt( (Ekin0 + (Epot0 - Epot1))/Ekin1 )
                #print "scalingfactor = %s" % scalingfactor
                assert abs(scalingfactor - 1.0) < 1.0e-1, "Total energy is not conserved, Etot(t) = %s, Etot(t+dt) = %s. The velocities would have to be scaled with a factor s=%s to artificially enforce energy conservation!" % (Ekin0+Epot0, Ekin1+Epot1, scalingfactor)
                #self.kinetic_energy *= scalingfactor**2
                self.velocities = scalingfactor * self.velocities

	def initiateTrajectory(self):

                if self.restart:
                        self.restartTrajectory()
                else:
                        # 
                        self.readCoordVel("dynamics.in")
                        
                # origin of energy axis
                self.zero=0.0 #1.0*self.energy[0]
                        
		if self.nat>1:        
			self.coordinates=self.shiftToCenterOfMass(self.coordinates)
			self.velocities=self.eliminateTranslation(self.velocities,self.momentum(self.velocities,self.masses))

			if self.writeflag=="Full":
				self.writeCoordVel_COM_noTrans()
			
                	angmom=self.angularMomentum(self.coordinates,self.velocities,self.masses)
                	inertia=self.momentOfInertia(self.coordinates,self.masses)
                	angvel=self.getAngularVelocity(angmom,inertia)

                	self.velocities=self.eliminateRotation(self.velocities,angvel,self.coordinates)
			
			if self.writeflag=="Full":
                                if self.restart == 0:
                                        self.writeVel_noRot()
		self.initializeInterface()
                
                self.kinetic_energy=self.getKineticEnergy(self.velocities,self.masses)
	
	def verlet(self):
		# velocity Verlet integation with constant energy; if self.mode=="T", with constant temperature
                self.initiateTrajectory()
                self.getData()
                ### metadynamics
                if self.mode=="M":
                        # if vg_center.dat is present, the metadynamics calculation is restarted
                        if os.path.exists("vg_centers.dat"):
                                restart = True
                        else:
                                restart = False
                        metaob = meta.metadynamics(symbols=self.symbols, restart=restart)
                        metaob.read_input()
                        metaob.set_pes_object(self.handler.pes)
                        metaob.initialize_vg()
                        
                ###
                        
		self.coord_old=self.coordinates

		self.getCoordVerlet()

		if self.nat>1:
                	self.coordinates=self.shiftToCenterOfMass(self.coordinates)

                if self.restart == 0:
                        # When restarting the trajectory, do not duplicate the last step
                        self.writeStep()

                self.time+=self.tstep

                if self.mode == "M":
                        print "== METADYNAMICS =="
                elif self.mode == "T":
                        print "== DYNAMICS AT CONSTANT TEMPERATURE =="
                else:
                        print "== DYNAMICS AT CONSTANT ENERGY =="
                        
		for i in range(1,self.nstep+1):
                        os.system("/bin/date >> timing")
                        self.old_accelerations=1.0*self.accelerations
                        self.oldenergy=1.0*self.energy[self.state]+self.kinetic_energy
                        # added by A.Humeniuk
                        self.kinetic_energy_old = self.kinetic_energy
                        self.potential_energy_old = self.energy[self.state]
                        self.energy_last = 1.0*self.energy
                        #
                        if i > 1:
			        self.getData()
                        ### metadynamics
                        if self.mode=="M":
                                tmp = metaob.get_forces(self.coordinates,i)
                                tmp = (tmp.transpose()/self.masses).transpose() # now we have accelerations
                                self.accelerations += tmp
                        ###
                        
			self.oldstate=self.state

                        if self.grounddyn == 1:
                                # The calculation of hopping probabilities is skipped, because
                                # the trajectory is forced on the ground state
                                pass
                        else:
                                self.c_old=1.0*self.c
                                self.getLocalDiabatization()
                                self.getNewState()
                                if self.decoherence_correction == 1:
                                        self.getDecoherenceE()
                                if (self.state<>self.oldstate):
                                        self.getUniformlyRescaledVelocities()

			self.getVelVerlet()

                        if self.nat>1:
				self.eliminateTransRot()
			self.kinetic_energy=self.getKineticEnergy(self.velocities,self.masses)

			if self.mode=="T" or self.mode=="M":
				self.scaleVelocitiesT()
                        # added by A.Humeniuk
                        if self.mode=="E" and self.artificial_energy_conservation == True:
                                self.scaleVelocitiesEconst()
			
                        self.coord_old=self.coordinates

			self.getCoordVerlet()

                        if self.nat>1:
				self.coordinates=self.shiftToCenterOfMass(self.coordinates)
                                    
                        self.time+=self.tstep

                        if i % self.output_step == 0:
                                # Don't write output for every time step, otherwise
                                # the disk gets fill up quickly.
                                self.writeStep()

                            
                print "FINISHED"

        ############################################################################################################
        # INPUT ROUTINES
        ############################################################################################################
        def readCoordVel(self, dynamics_in):
                # read initial conditions (coordinates and velocities) and assign atom symbols and masses
                atomlist_coords, atomlist_vels = XYZ.read_initial_conditions(dynamics_in, units="bohr")
                self.nat = len(atomlist_coords)
                self.coordinates = np.zeros((self.nat, 3))
                self.velocities = np.zeros((self.nat, 3))
                self.masses = np.zeros(self.nat)
                self.symbols = []
                for i,(Zi,posi) in enumerate(atomlist_coords):
                        atname = AtomicData.atom_names[Zi-1]
                        self.symbols.append(atname)
                        self.masses[i] = AtomicData.atom_masses[atname]
                        self.coordinates[i,:] = posi
                        self.velocities[i,:] = atomlist_vels[i][1]

		self.totalmass=float(sum(self.masses))

        def restartTrajectory(self):
                # function for reading last line of file
                def last_line(filename):
                        with open(filename) as fh:
                                lines = fh.readlines()
                                last = lines[-1]
                                vals = map(float, last.strip().split())
                                return vals
                # read coordinates and velocities of last step
                self.readCoordVel("dynamics.in0")
                # read coefficients of electronic wavefunction
                self.c = np.loadtxt("coefficients.dat0").view(complex)
                # normalize coefficients
                self.c /= la.norm(self.c)
                # read current time step and state
                self.time, self.state = last_line("state.dat")
                # convert time from fs to au
                self.time = self.time * self.fs_to_au 
                self.state = int(self.state)

        ############################################################################################################
        # OUTPUT ROUTINES
        ############################################################################################################

	def writeCoordVel_COM_noTrans(self):
		output=open("dynamics.out","a")
                output.write("           Coordinates in the center of mass frame:\n")
                for i in range(self.nat):
                         output.write("%s   %20.12f  %20.12f  %20.12f\n" %(self.symbols[i],self.coordinates[i][0],self.coordinates[i][1],self.coordinates[i][2]))

                output.write("           Velocities after eliminating translation:            \n")
                for i in range(self.nat):
                        output.write("   %20.12f  %20.12f  %20.12f\n" %(self.velocities[i][0],self.velocities[i][1],self.velocities[i][2]))
		output.close()

	def writeVel_noRot(self):
		output=open("dynamics.out","a")
                output.write("           Velocities after eliminating rotation:\n")
                for i in range(self.nat):
                         output.write("   %20.12f  %20.12f  %20.12f\n" %(self.velocities[i][0],self.velocities[i][1],self.velocities[i][2]))
		output.close()

        def writeStep(self):

		if self.writeflag=="Full":
			self.writeXYZ()
		elif self.writeflag=="xyz":
			self.writeXYZ()

                self.writeActualEnergies()

                if self.printcoup==True:
                        self.writeNACoupling()
                self.writeEnergies()
                if hasattr(self, "Hdiab"):
                        self.writeLocalDiabaticEnergies()
                        self.writeLocalDiabaticCouplings()
                self.writeState()
                self.writeTimeSeries()

                if self.printcoeff>0:
                        self.writeCoefficients()
                if self.printcoeff==3:
                        self.writeCoefficientsReIm()
                if self.printcoeff==4:
                        self.writeCoherences()
		
	def writeXYZ(self):
		outfile=open("dynamics.xyz","a")
		outfile.write(str(self.nat)+"\n")
		outfile.write("\n")
		tmp=self.au_to_ang*self.coordinates
		for i in range(len(tmp)):
                   outfile.write("%s  %20.12f  %20.12f  %20.12f\n" %(self.symbols[i],tmp[i][0],tmp[i][1],tmp[i][2]))
		outfile.close()
		self.writeRestart()
                
	def writeRestart(self):
		"""
		save coordinates and velocities (in au) and electronic coefficients for the current step
		so that the calculation can be restarted 
		"""
		outfile=open("dynamics.in0","w")
		outfile.write(str(self.nat)+"\n")
		for i in range(len(self.coordinates)):
                   outfile.write("%s  %20.12f  %20.12f  %20.12f\n" %(self.symbols[i],self.coordinates[i][0],self.coordinates[i][1],self.coordinates[i][2]))
		for i in range(len(self.velocities)):
                   outfile.write("%20.12f  %20.12f  %20.12f\n" %(self.velocities[i][0],self.velocities[i][1],self.velocities[i][2]))
		outfile.close()

                # electronic coefficients
                np.savetxt("coefficients.dat0", self.c.view(float))
                
	def writeNACoupling(self):
               #for i in range(self.firststate,self.firststate+4):
               #         for j in range(i+1,self.firststate+4):
               for i in range(self.nstates):
                        for j in range(i+1,self.nstates):
                                file=open("nonadiabatic"+str(i)+str(j)+".dat","a")
                                file.write(str(self.time/self.fs_to_au)+" "+str(self.nonad_scalar[i][j])+"\n")
                                file.close()

	def writeEnergies(self):
		for i in range(self.nstates):
			file=open("energy_"+str(i)+".dat","a")
			file.write(str(self.time/self.fs_to_au)+" "+str(self.energy[i]-self.zero)+"\n")
			file.close()

        def writeActualEnergies(self):
                file=open("curr_energy.dat","a")
                etot = self.kinetic_energy + self.energy[self.state]
                file.write(str(self.time/self.fs_to_au)+" "+str(self.kinetic_energy)+" "+str(self.energy[self.state])+" "+str(etot)+"\n")
                file.close()

	def writeCoefficients(self):
		norm=0.0
		for i in range(self.nstates):
			c_sq=abs(self.c[i])**2
			norm+=c_sq
			file=open("coeff_"+str(i)+".dat","a")
			file.write(str(self.time/self.fs_to_au)+" "+str(c_sq)+"\n")
			file.close()
		file=open("coeff_norm.dat","a")
		file.write(str(self.time/self.fs_to_au)+" "+str(norm)+"\n")
		file.close()

        def writeCoherences(self):
                for i in range(self.nstates):
                        for j in range(i+1,self.nstates):
                                c_sq=abs(conjugate(self.c[i])*self.c[j])
                                file=open("coherence_"+str(i)+"_"+str(j)+".dat","a")
                                file.write(str(self.time/self.fs_to_au)+" "+str(c_sq)+"\n")
                                file.close()

        def writeCoefficientsReIm(self):
                for i in range(self.nstates):
                        file1=open("coeff_re_"+str(i)+".dat","a")
                        file1.write(str(self.time/self.fs_to_au)+" "+str(self.c[i].real)+"\n")
                        file1.close()
                        file2=open("coeff_im_"+str(i)+".dat","a")
                        file2.write(str(self.time/self.fs_to_au)+" "+str(self.c[i].imag)+"\n")
                        file2.close()
                        
	def writeState(self):
		file=open("state.dat","a")
                file.write(str(self.time/self.fs_to_au)+" "+str(self.state)+"\n")
                file.close()

        def writeLocalDiabaticEnergies(self):
                fh = open("local_diabatic_energies.dat", "a")
                line = str(self.time/self.fs_to_au)+" "
                for i in range(self.nstates):
                        line += " "+str(self.Hdiab[i,i] - self.zero)
                line += "\n"
                fh.write(line)
                fh.close()
        def writeLocalDiabaticCouplings(self):
                fh = open("local_diabatic_couplings.dat", "a")
                line = str(self.time/self.fs_to_au)+" "
                for i in range(self.nstates):
                        for j in range(i+1,self.nstates):
                                line += " "+str(self.Hdiab[i,j])
                line += "\n"
                fh.write(line)                
                fh.close()
        def writeTimeSeries(self):
                # This function writes additional quantities
                # such as particle-hole charges along the trajectory
                # The function call is delegated to the handler
                if hasattr(self.handler, "writeTimeSeries"):
                        self.handler.writeTimeSeries(self, time_series=self.time_series)
#####################################################################################################################

if __name__ == "__main__":
        import sys
        
        usage  = "Usage: %s\n" % os.path.basename(sys.argv[0])
        usage += "  run surface hopping dynamics.\n\n"
        usage += "  Required Input Files:\n"
        usage += "    - `dynamics.in` with the initial coordinates and velocities (in bohr)\n"
        usage += "    - `dftbaby.cfg` with the options controlling the TD-DFTB calculation and dynamics\n"
        usage += "  Output Files:\n"
        usage += "    - `dynamics.xyz`: geometries for each time step (in Angstrom)\n"
        usage += "    - `state.dat`: current electronic state\n"
        usage += "    - `energy_I.dat`: total energy of electronic state I in Hartree vs. time in fs.\n"
        usage += "        The ground state energy at time t=0 is subtracted from all later energies.\n"
        usage += "    - `nonadiabaticIJ.dat`: scalar non-adiabatic coupling between states I and J.\n"

        # This wrapper makes the optional parameters of the python function __init__ visible
        # as optional command line argument.
        parser = optparse.OptionParserFuncWrapper(
                [
                        # options for MD
                        MolecularDynamics.__init__,
                        # electronic structure options
                        DFTB2.__init__, 
                        DFTB2.runSCC,
                        SolventCavity.__init__,
                        LR_TDDFTB.getEnergies
                ],
                usage, section_headers=["SurfaceHopping", "DFTBaby"])
        # extract optional parameters from command line
        (options,args) = parser.parse_args(MolecularDynamics.__init__)
        
        md=MolecularDynamics(**options)
        # run molecular dynamics - this takes a while
        md.verlet()

        print "FINISHED"
