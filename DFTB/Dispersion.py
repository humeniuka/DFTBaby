"""
Grimme's semiempirical dispersion correction

"Semiempirical GGA-Type Density Functional Constructred with a Long-Range Dispersion Correction"
Grimme,S. J. Comput. Chem. 2006 27 1787-1799

The C6 parameters and vdW radii can be found in DFTB.AtomicData
"""
from DFTB import XYZ, AtomicData

import numpy as np
import numpy.linalg as la

class DispersionCorrection:
    def __init__(self, atomlist):
        """
        This function is called at the start of the program. It may initialize necessary parameters based on the atoms present in the molecule.

        Parameters:
        ===========
        atomlist: a list of tuples (Zi, [xi,yi,zi]) with the atomic number
          and cartesian coordinates in bohr for each atom
        """
        # find unique atom types
        atomtypes = list(set([Zi for (Zi,posi) in atomlist]))
        atomtypes.sort()
        # C6 coefficients in [J * nm^6 * mol^(-1)]
        C6 = np.array([AtomicData.Grimme_C6[AtomicData.atom_names[Zi-1]] for Zi in atomtypes])
        # convert to [Hartree * bohr^6]
        C6 *= 17.34525539265301
        # atomic van der Waals radii in Angstrom
        R0 = np.array([AtomicData.Grimme_R0[AtomicData.atom_names[Zi-1]] for Zi in atomtypes])
        # convert to bohr
        R0 /= AtomicData.bohr_to_angs
        
        self.C6_pairs = {}  # C6 coefficients between atom pair (Zi,Zj)
        self.Rr = {}
        for i,Zi in enumerate(atomtypes):
            for j,Zj in enumerate(atomtypes):
                # Grimme uses the geometric mean in eqn. (13)
                self.C6_pairs[(Zi,Zj)] = np.sqrt(C6[i]*C6[j])
                # Rr is the sum of atomic vdW radii
                self.Rr[(Zi,Zj)] = R0[i] + R0[j]
        # exponent in the damping function, f_dmp(Rij) = 1/(1+exp(-d*(Rij/Rr-1)))
        self.d = 20.0
        # The scaling factor really depends on the functional, because DFTB is parametrized based on PBE
        # the value 0.75 should be used. However, this leads to very strong attraction between pentacene and C60.
        self.s6 = 0.75
    def getEnergy(self, atomlist):
        """
        Parameters:
        ===========
        atomlist: a list of tuples (Zi, [xi,yi,zi]) with the atomic number
          and cartesian coordinates in bohr for each atom

        Returns:
        ========
        energy: dispersion correction in Hartree
        """
        #print "energy of dispersion correction"

        en = 0.0 # energy from dispersion correction in Hartree
        for i,(Zi,posi) in enumerate(atomlist):
            for j,(Zj,posj) in enumerate(atomlist):
                if j <= i:
                    continue
                C6 = self.C6_pairs[(Zi,Zj)]
                Rr = self.Rr[(Zi,Zj)]
                Rij = la.norm( np.array(posj) - np.array(posi) )
                a = self.d * (Rij/Rr - 1.0)
                f = 1.0/(1.0 + np.exp(-a))
                en -= C6/pow(Rij,6) * f
        # scaling factor
        en *= self.s6
        return en
    def getGradient(self, atomlist, distances, directions):
        """
        Parameters:
        ===========
        atomlist: a list of tuples (Zi, [xi,yi,zi]) with the atomic number
          and cartesian coordinates in bohr for each atom
        distances: a Nat x Nat dimensional matrix with the distances between
          the atoms in bohr, distances[i,j] = |R_i - R_j|
        directions: directions[i,j,:] is the unit vector pointing from
          atom j to atom i

        Returns:
        ========
        grad: a 1D numpy array with the 3*Nat components of the gradient of the dispersion energy in Hartree/bohr
        """
        #print "gradient of dispersion correction"
        Nat = len(atomlist)
        grad = np.zeros(3*Nat)
        for i,(Zi,posi) in enumerate(atomlist):
            for j,(Zj,posj) in enumerate(atomlist):
                if j == i:
                    continue
                Rij = distances[i,j]
                eij = directions[i,j,:]
                C6 = self.C6_pairs[(Zi,Zj)]
                Rr = self.Rr[(Zi,Zj)]
                # damping function
                a = self.d * (Rij/Rr - 1.0)
                f = 1.0/(1.0 + np.exp(-a))
                Vij_deriv = - C6/pow(Rij,6) * f * (self.d/Rr * (1.0-f) - 6.0/Rij)
                grad[i*3:i*3+3] += Vij_deriv * eij
        grad *= self.s6
        return grad

