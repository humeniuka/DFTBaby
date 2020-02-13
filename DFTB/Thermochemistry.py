"""
computes thermodynamics quantities such as
  internal thermal energy E
  entropy S
  and heat capacity

The formulae are taken from a Gaussian whitepaper (http://www.gaussian.com/g_whitepap/thermo.htm)
and Atkins' 'Physical Chemistry'. Sometimes the same quantities are computed from different
functional forms and are compared as a double check.

TODO: rewrite formulae in terms of beta=1/(kB*T) and compare them with those in Atkins' book
"""
import numpy as np
import numpy.linalg as la

from DFTB import XYZ, AtomicData
from DFTB.Modeling import MolecularCoords as MolCo

class Thermochemistry:
    def __init__(self, atomlist, E0, vib_freq, symmetry_group,
                 pressure=AtomicData.atm_pressure,
                 temperature=AtomicData.satp_temperature):
        """
        temperature in Kelvin
        """
        self.E0 = E0
        self.P = pressure
        self.T = temperature
        self.vib_freq = vib_freq
        self.symmetry_group = symmetry_group
        
        self.atomlist = atomlist
        Nat = len(atomlist)
        self.masses = AtomicData.atomlist2masses(self.atomlist)
        # compute rotational constants from tensor of inertia
        x0 = XYZ.atomlist2vector(atomlist)
        # shift the origin to the center of mass
        x0_shift = MolCo.shift_to_com(x0, self.masses)
        # diagonalize the tensor of inertia to obtain the principal moments and 
        # the normalized eigenvectors of I
        Inert = MolCo.inertial_tensor(self.masses, x0_shift)
        principle_moments, X = la.eigh(Inert)
        Iaa,Ibb,Icc = np.sort(abs(principle_moments))
        print "principle moments of inertia (in a.u.): %s %s %s" % (Iaa,Ibb,Icc)
        # In a linear triatomic molecule we have Icc = Ibb > Iaa = 0
        self.is_linear = False
        if abs(Icc-Ibb)/abs(Icc) < 1.0e-6 and abs(Iaa)/abs(Icc) < 1.0e-6:
            self.is_linear = True
            print "Molecule is linear"
        # Determine the type of rotor
        if self.is_linear == True:
            self.rotor_type = "linear rotor"
        else:
            if abs(Icc-Ibb)/abs(Icc) < 1.0e-4 and abs(Icc-Iaa)/abs(Icc) < 1.0e-4:
                # three equal moments of inertia
                self.rotor_type = "spherical rotor"
            elif abs(Icc-Ibb)/abs(Icc) < 1.0e-4 or abs(Ibb-Iaa)/abs(Icc) < 1.0e-4:
                # two equal moments of inertia
                self.rotor_type = "symmetric rotor"
            else:
                self.rotor_type = "asymmetric rotor"
                
        # rotational constants
        self.rotational_constants = 1.0/(2.0 * principle_moments + 1.0e-20) # avoid division by zero error for a linear molecule, the invalid values are not actually used
        # symmetry number
        if symmetry_group == None:
            print "Symmetry group unknown, setting rotational symmetry number to 1"
            self.sigma_rot = 1
        else:
            self.sigma_rot = symmetry_group.rotational_symmetry_number()
        # beta = 1/(kB*T)
        kB = AtomicData.kBoltzmann
        self.beta = 1.0/(kB*self.T)
    def translational(self):
        kB = AtomicData.kBoltzmann
        kBT = kB * self.T
        # total mass of the molecule
        M = np.sum(self.masses[::3])
        print "Molecular mass = %s amu" % (M * AtomicData.aumass2amu)
        # translational partition function
        self.qt = pow(2.0*np.pi, -3.0/2.0) * pow(M * kBT,3.0/2.0) * kBT / self.P
        self.thermal_wavelength = 2.0*np.pi/np.sqrt(2.0*np.pi * M * kBT)  # in atomic units, hbar is set to 1, so h has to be set to 2pi
        # compute volume per particle from pressure for an ideal gas
        self.volume = kBT / self.P
        # translational partition function, again, computed differently
        qt = self.volume / pow(self.thermal_wavelength, 3)
        assert np.abs(self.qt-qt)/abs(self.qt) < 1.0e-10
        # contribution from translations
        self.Et = 3.0/2.0 * kBT                     # internal energy
        self.St = kB * (np.log(self.qt) + 5.0/2.0)  # entropy
        # Sackur-Tetrode equation for indistinguishable particles, where Q = q^N/N!
        St = self.Et / self.T + kB * (np.log(qt) + 1.0)
        assert np.abs(self.St-St)/abs(self.St) < 1.0e-10
        self.Ct = 3.0/2.0 * kB                      # heat capacity at constant volume
    def rotational(self):
        # The rotational partition function is defined as qR = sum_J (2J+1) exp[-beta * hc * B*J*(J+1)]
        # but the approximation is made that when kBT is much larger than the separation between neighbouring
        # rotational energy levels, the sum over J can be converted into an integral.
        kB = AtomicData.kBoltzmann
        kBT = kB * self.T
        if self.is_linear == True:
            C = self.rotational_constants[2]
            print "rotational constant:  %10.6f cm^-1" % (C*AtomicData.hartree_to_wavenumbers)
            # rotational partition function (integral approximation)
            self.qr = kBT / (self.sigma_rot * C)  #
            # compute partition function exactly by summing over J
            self.qr_sum = 0.0
            J = 0
            while True:
                enJ = J*(J+1)*C  # C = 1/2I
                dqr = (2*J+1) * np.exp(-self.beta * enJ)
                self.qr_sum += dqr
                if J > 10:
                    if dqr/self.qr_sum < 1.0e-10:
                        break
                J += 1
            self.qr_sum /= self.sigma_rot

            self.Er = kBT
            self.Sr = kB * (np.log(self.qr) + 1.0)  # ?
            self.Cr = kB
        else:
            # asymmetric top
            A,B,C = self.rotational_constants/kB
            f = AtomicData.hartree_to_wavenumbers * kB
            print "rotational constants:  %10.6f  %10.6f  %10.6f  cm^-1" % (A*f,B*f,C*f)
            self.qr = np.sqrt(np.pi)/self.sigma_rot * pow(self.T,3.0/2.0) / np.sqrt(A*B*C)
            self.qr_sum = -1.0 # How dow I calculate the partition function for general A,B and C?
            #
            self.Er = 3.0/2.0 * kBT
            self.Sr = kB * (np.log(self.qr) + 3.0/2.0)  # ?
            self.Cr = 3.0/2.0 * kB
    def vibrational(self):
        kB = AtomicData.kBoltzmann
        kBT = kB * self.T

        self.Ev = 0.0
        self.Sv = 0.0
        self.Cv = 0.0
        #
        self.Evib0 = 0.0   # vibrational zero-point energy measured from the bottom of the potential energy surface
        nvib = len(self.vib_freq)
        self.qv = 1.0
        # iterate over all vibrational modes
        for i in range(0, nvib):
            # energy levels of QM harmonic oscillator   En = h*nu * (1/2 + n)
            hnu = self.vib_freq[i].real
            # vibrational partition function for i-th mode
            qvi = 1.0/(1.0 - np.exp(-self.beta * hnu))
            # qvi = sum_n exp(-b * (En-E0)) = sum_n exp(-b*hnu * n)) = 1/(1 - exp(-b*hnu))
            # The partition function is the product of the partition functions for each vibration
            # qv = qv(1)*qv(2)*...*qv(Nvib)
            self.qv *= qvi
            self.Evib0 += hnu / 2.0
            # Boltzmann factor
            fi = np.exp(-self.beta * hnu)
            # U-U0 = - 1/qv * [dq/dbeta]_(constant volume)
            self.Ev += hnu * fi / (1.0 - fi)
            # Cv = (dU/dT)_{V=const} = -k * b^2 dU/db
            self.Cv += 1.0/4.0 * (self.beta*hnu)**2 / np.sinh(self.beta*hnu/2.0)**2
            # same formula
            #self.Cv += pow(hnu/kBT,2) * pow(np.exp(-hnu/(2*kBT))/(1.0 - np.exp(-hnu/kBT)),2)
            
        # S = (U-U0)/T + kB * log(qv)
        self.Sv = self.Ev / self.T + kB * np.log(self.qv)

        self.Sv *= kB
        self.Ev *= kB
        self.Cv *= kB
    def electronic(self):
        kB = AtomicData.kBoltzmann
        g0 = 1.0                       # The degeneracy of the ground state is 1, since my DFTB can only calculated Singlet ground states.
        self.qe = g0  # only the ground state is thermally accessible,
                       # assuming 1st excited state is much higher in energy (Born-Oppenheimer approximation)
                       # so that it is not accessible at room temperature.
        self.Ee = 0.0  
        self.Se = kB * np.log(self.qe)
        self.Ce = 0.0
    def calculate(self):
        print "**********************************"
        print "* Thermochemistry (not tested)   *"
        print "**********************************"
        print ""
        self.translational()
        self.rotational()
        self.vibrational()
        self.electronic()
        
        self.E = self.Et + self.Er + self.Ev + self.Ee
        self.S = self.St + self.Sr + self.Sv + self.Se
        self.C = self.Ct + self.Cr + self.Cv + self.Ce

        f = AtomicData.hartree_to_kcalmol
        if self.symmetry_group != None:
            print "point group: %s" % self.symmetry_group.name()
        else:
            print "point group: unknown"
        print "rotational symmetry number: %s (check this, it's probably wrong!)" % self.sigma_rot
        print "Molecule is a %s"  % self.rotor_type
        print "Temperature: %s K       Pressure: %s atm" % (self.T, self.P / AtomicData.atm_pressure)
        print "thermal wavelength:  %10.6f bohr" % self.thermal_wavelength
        print "volume per particle:   %.2e bohr^3" % self.volume
        print "                  Partition Function Q"
        print "                   integral approx.            summation"
        print "translational qt:  %e                    - " % (self.qt)
        print "rotational qr:     %e                %e    " % (self.qr, self.qr_sum)
        print "vibrational qv:    %e" % (self.qv)
        print "total qt*qr*qv:    %e" % (self.qt*self.qr*self.qv)
        print "                  Thermal Energy E              Heat Capacity Cv              Entropy S"
        print "                    kcal/mol                      kcal/(mol K)                 kcal/(mol K)"
        print "translational:  %10.6f                    %10.6f                    %10.6f" % (self.Et*f, self.Ct*f, self.St*f)
        print "rotational:     %10.6f                    %10.6f                    %10.6f" % (self.Er*f, self.Cr*f, self.Sr*f)
        print "vibrational:    %10.6f (incl. Evib0)      %10.6f                    %10.6f" % ((self.Ev+self.Evib0)*f, self.Cv*f, self.Sv*f)
        print "electronic:     %10.6f                    %10.6f                    %10.6f" % (self.Ee*f, self.Ce*f, self.Se*f)
        print "total:          %10.6f                    %10.6f                    %10.6f" % (self.E*f, self.C*f, self.S*f)
        print ""
        print "electronic energy E0 (full DFTB energy):     %10.6f kcal/mol" % (self.E0*f)
        print "vibrational zero-point energy Evib0:         %10.6f kcal/mol" % (self.Evib0*f)
        print "electronic and zero-point energy E0+Evib0:   %10.6f kcal/mol" % ((self.E0+self.Evib0)*f)
        H = self.E0+self.Evib0+self.E
        print "enthalpy/heat H=E0+Evib0+E:                  %10.6f kcal/mol" % (H*f)
        G = H - self.T*self.S
        print "Gibbs free energy G=H-TS:                    %10.6f kcal/mol" % (G*f)
        print ""
        
