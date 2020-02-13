"""
A cavity potential confines the molecules to a sphere of radius
r0 and pushes molecules back that would otherwise fly away.

           /  0                      for r <= r0
Vcav(r) = {
           \  1/2 * k * m * (r-r0)^2     for r > r0

The potential is proportionate to the mass m of the particle so that all particles experience
the same acceleration at the boundary of the cavity.
"""

import numpy as np
import numpy.linalg as la

from DFTB.AtomicData import atom_names, atom_masses

class Cavity:
    def __init__(self, r0=0.0, k=0.0):
        self.r0 = r0   # cavity radius
        self.k = k     # force constant
    def getEnergy(self, atomlist):
        enCav = 0.0
        for (Zi,posi) in atomlist:
            ri = la.norm(posi)
            if ri > self.r0:
                mi = atom_masses[atom_names[Zi-1]]
                enCav += 0.5*self.k * mi * (ri-self.r0)**2
        return enCav
    def getGradient(self, atomlist):
        Nat = len(atomlist)
        grad = np.zeros(3*Nat)
        for i,(Zi,posi) in enumerate(atomlist):
            ri = la.norm(posi)
            if ri > self.r0:
                mi = atom_masses[atom_names[Zi-1]]
                # unit vector in r direction
                ei = np.asarray(posi) / ri
                grad[3*i:3*(i+1)] = self.k * mi * (ri-self.r0) * ei
        return grad
            
