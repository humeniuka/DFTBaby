"""
ring puckering coordinates as defined in

  Cremer,D.; Pople,J.
  'A General Definition of Ring Puckering Coordinates'
  JACS, 1975, 97, 1354-1358.
"""
from DFTB import XYZ, AtomicData
from cmath import phase

import numpy as np
import copy

class Ring:
    def __init__(self, atomlist, ring_atoms=[], units="bohr"):
        """
        creates a molecule with a ring

        Parameters
        ----------
        atomlist: list of tupes (Zi,[xi,yi,zi]) with atomic positions

        Optional
        --------
        ring_atoms: indeces of atoms forming the ring
        bohr: units of coordinates [xi,yi,zi]
        """
        if units == "Angstrom":
            x = XYZ.atomlist2vector(atomlist) / AtomicData.bohr_to_angs
            self.atomlist = XYZ.vector2atomlist(x, atomlist)
        else:
            self.atomlist = atomlist
        self.ring_atoms = ring_atoms
    def getGeometry(self):
        return copy.deepcopy(self.atomlist)
    def set_ring_atoms(ring_atoms):
        self.ring_atoms = ring_atoms
    def center_of_coords(self):
        """find the center of coordinates of the ring atoms,
        equals the center of mass in a homocyclic ring"""
        M = 0.0
        R = np.zeros(3)
        for ra in self.ring_atoms:
            Zr,posr = self.atomlist[ra]
            R += np.array(posr)
            M += 1.0
        R /= M
        return R
    def relative_coords(self):
        """find relative positions of ring atoms relative to the center of coordinates"""
        R = self.center_of_coords()
        Rrel = [np.array(self.atomlist[ra][1]) - R for ra in self.ring_atoms]
        return Rrel
    def ring_plane(self):
        """find the normal that defines the ring plane"""
        N = len(self.ring_atoms)
        Rrel = self.relative_coords()
        R1 = reduce(lambda r1,r2: r1+r2, [Rrel[i]*np.sin(2.0*np.pi*i/float(N)) for i in range(0, N)])
        R2 = reduce(lambda r1,r2: r1+r2, [Rrel[i]*np.cos(2.0*np.pi*i/float(N)) for i in range(0, N)])
        cr = np.cross(R1,R2)
        n = cr/np.sqrt(np.dot(cr,cr))
        return n
    def displacements(self):
        n = self.ring_plane()
        Rrel = self.relative_coords()
        displZ = [np.dot(rrel,n) for rrel in Rrel]
        N = len(self.ring_atoms)
        assert(abs(sum(displZ)) < 10e-10)
        assert(abs(sum([displZ[i]*np.cos(2.0*np.pi*i/float(N)) for i in range(0, N)])) < 10e-10)
        assert(abs(sum([displZ[i]*np.sin(2.0*np.pi*i/float(N)) for i in range(0, N)])) < 10e-10)
        """displacements along the normal of the ring atoms"""
        return displZ
    def get_puckering_coords(self):
        """calculate the puckering amplitudes and phases from the geometry of the ring atoms"""
        N = len(self.ring_atoms)
        z = self.displacements()
        mmax = maximum_m(N)
        puck_ampl = []
        puck_phase = []
        for m in range(2,mmax+1):
            qm = np.sqrt(2.0/N)*reduce(lambda x,y:x+y, \
                                     [z[i]*np.exp(-1.0j * 2.0*np.pi*m*i/float(N)) for i in range(0, N)])
            puck_ampl.append(abs(qm))
            puck_phase.append(phase(qm))
        if N % 2 == 0:
            """if N is even one has an additional puckering coordinate"""
            qN2 = 1.0/np.sqrt(N)*reduce(lambda x,y:x+y, [pow(-1,i)*z[i] for i in range(0, N)])
            puck_ampl.append(qN2)
            puck_phase.append(None)
        return puck_ampl, puck_phase
    def set_puckering_coords(self, puck_coords):
        puck_ampl, puck_phase = puck_coords
        """displace atoms according to puckering coordinates"""
        N = len(self.ring_atoms)
        #assert N-3 == 2*len(puck_ampl), "wrong number of puckering amplitudes"
        z = np.zeros(N)
        mmax = maximum_m(N)
        """z-displacements"""
        pre_fac = np.sqrt(2.0/N)
        for i in range(0, N):
            for m in range(2, mmax+1):
                z[i] += pre_fac*puck_ampl[m-2]*np.cos(puck_phase[m-2] + 2.0*np.pi*m*i/float(N))
        if N % 2 == 0:
            for i in range(0, N):
                z[i] += np.sqrt(1.0/N)*puck_ampl[-1]*pow(-1,i)
        zold = self.displacements()
        dz = map(lambda x,y: x-y, z, zold)
        """displace the ring atoms by an amount dz along the normal
        to achieve the desired puckering geometry"""
        n = self.ring_plane()
        atomlist = copy.deepcopy(self.atomlist)
        for i in range(0, N):
            Zr,posr  = self.atomlist[self.ring_atoms[i]]
            posr = np.array(posr) + dz[i]*n
            atomlist[self.ring_atoms[i]] = (Zr,posr)
        self.atomlist = atomlist
            
def maximum_m(N):
    if N % 2 == 0:
        """N is even"""
        mmax = N/2-1
    else:
        """N is odd"""
        mmax = (N-1)/2
    return mmax

def test_puckering_coords():
    O = 8
    C = 6
    """reproduce puckering coordinates of Furanoid ring in the paper by Cremer/Pople"""
    Furanoid = Ring(    [(O, (0.,     1.2111,  -0.0189)),
                         (C, (1.1622, 0.4349,  0.1461 )),
                         (C, (0.7425, -1.0012, -0.2174)),
                         (C, (-0.7221,-1.0309, 0.2057)),
                         (C, (-1.1826, 0.3861, -0.1154))],
                        ring_atoms = [0,1,2,3,4], units="Angstrom")
    puck_ampl, puck_phase = Furanoid.get_puckering_coords()
    assert(abs(puck_ampl[0]*AtomicData.bohr_to_angs - 0.353) < 0.001)
    assert(abs(puck_phase[0]+np.pi - 4.627 < 0.001))
    print "q_2 = %s"  % (puck_ampl[0]*AtomicData.bohr_to_angs)     # should be 0.353
    print "phi_2 = %s" % (puck_phase[0]+np.pi)       # should be 265.1 deg = 4.627 rad

    Furanoid.set_puckering_coords( (puck_ampl, puck_phase) )
    """reproduce puckering for pyranoid ring of sucrose"""
    Pyranoid = Ring([(O, (0.0,     1.3839,    0.1976)),
                     (C, (1.1997,  0.7624,    -0.2106)),
                     (C, (1.2356,  -0.7040,   0.2393)),
                     (C, (0.0110,  -1.4564,  -0.2550)),
                     (C, (-1.2300, -0.7208,  0.2420)),
                     (C, (-1.2164, 0.7350,   -0.2133))],
                    ring_atoms = [0,1,2,3,4,5], units="Angstrom")
    puck_ampl, puck_phase = Pyranoid.get_puckering_coords()
    print "q_2 = %s" % (puck_ampl[0]*AtomicData.bohr_to_angs)       # should be 0.050
    print "q_3 = %s" % (puck_ampl[1]*AtomicData.bohr_to_angs)       # should be 0.554
    print "phi_2 = %s" % (puck_phase[0]+2.0*np.pi)    # should be 183.7 deg = 3.2062 rad

if __name__ == "__main__":
    import sys
    test_puckering_coords()

    # Test puckering on Furan
    O = 8  # atomic numbers of oxygen and carbon
    C = 6
    ring = Ring(    [(O, (0.,     1.2111,  -0.0189)),
                     (C, (1.1622, 0.4349,  0.1461 )),
                     (C, (0.7425, -1.0012, -0.2174)),
                     (C, (-0.7221,-1.0309, 0.2057)),
                     (C, (-1.1826, 0.3861, -0.1154))],
                    ring_atoms = [0,1,2,3,4], units="Angstrom")
    
    puckered_geometries = [ring.getGeometry()]

    for phi2 in np.linspace(0.0, 2.0*np.pi, 10):
        for q2 in np.linspace(0.0, 2.0, 10):
            puck_ampl = [q2]
            puck_phase = [phi2]
            ring.set_puckering_coords( (puck_ampl, puck_phase) )
            puckered_geometries.append( ring.getGeometry() )

    xyz_file = "/tmp/furan_puckering.xyz"
    XYZ.write_xyz(xyz_file, puckered_geometries)
    print "different puckered geometries of furan written to '%s'" % xyz_file

