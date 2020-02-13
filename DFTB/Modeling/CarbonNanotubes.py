"""
construct capped carbon nanotubes
"""
from DFTB import XYZ, AtomicData

import numpy as np
import numpy.linalg as la
from scipy import optimize

class CarbonNanotube:
    def __init__(self, n1,n2, a=4.6487, L=4*18.897):
        """
        Parameters:
        ===========
        a: magnitude of the primitive lattice vectors in bohr
        L: length of tube in bohr
        n1, n2: chiral indeces
        """
        # 
        self.n1 = n1
        self.n2 = n2
        # primitive lattice vectors 
        # This lattice is triangular 
        self.a1 = a*np.array([np.sqrt(3.0)/2.0, -1.0/2.0, 0.0])
        self.a2 = a*np.array([np.sqrt(3.0)/2.0,  1.0/2.0, 0.0])
        # position of carbon atoms
        self.c1 = (self.a1 + self.a2)/3.0
        self.c2 = (self.a1 + self.a2)*2.0/3.0
        # reciprocal lattice vectors satisfy ai*bj = 2 pi * delta_ij
        self.b1 = 2*np.pi/a * np.array([1.0/np.sqrt(3.0),-1.0, 0.0])
        self.b2 = 2*np.pi/a * np.array([1.0/np.sqrt(3.0), 1.0, 0.0])
        # chiral vector => circumference
        self.chiral_vector = self.n1*self.a1 + self.n2*self.a2
        self.chiral_norm = la.norm(self.chiral_vector)
        self.chiral_axis = self.chiral_vector / self.chiral_norm
        # chiral angle eqn. 2.1 from Reich's book
        self.theta = np.arccos( (n1+0.5*n2)/np.sqrt(n1**2+n1*n2+n2**2) )
        # tube radius
        self.R = la.norm(self.chiral_vector)/(2*np.pi)
        # the vector orthogonal to the chiral vector
        if self.n1 != 0:
            m2 = 1.0
            m1 = - self.n2/float(self.n1) * m2
        else:
            m2 = 0.0
            # m1 can be anything
            m1 = 1.0
        self.tube_vector = m1*self.b1 + m2*self.b2
        # normalze tube axis to L
        self.tube_axis = self.tube_vector / la.norm(self.tube_vector)
        self.L = L
        self.tube_vector = self.L*self.tube_axis
        assert abs(np.dot(self.chiral_axis, self.tube_axis)) < 1.0e-10
        print "a1            = %s" % self.a1
        print "a2            = %s" % self.a2
        print "chiral_vector = %s" % self.chiral_vector
        print "tube vector   = %s" % self.tube_vector 
    def createTube(self):
        # 
        def plane2cylinder(coords):
            u = np.dot(coords,self.chiral_axis)/self.chiral_norm
            v = np.dot(coords, self.tube_axis)/self.L
            if (0 <= u < 1) \
                    and (0 <= v < 1):
                # map (u,v) to cylinder
                z = v*self.L
                x = self.R*np.cos(2*np.pi*u)
                y = self.R*np.sin(2*np.pi*u)
                return np.array([x,y,z])
            # outside
            return None
        atomlist = []
        thomson_tube = [] # Thomson points on the tube
        for i in range(-30, 30):
            for j in range(-30, 30):
                # lattice vector
                vij = i*self.a1 + j*self.a2
                # carbon atoms in sublattices A and B
                cA = vij + self.c1
                cB = vij + self.c2
                # carbon lattice
                for coords_plane in [cA, cB]:
                    # check whether atoms are inside the rectangle
                    # spanned by the chiral vector and the tube vector
                    # t are the cartesian coordinates on the tube
                    coords_tube = plane2cylinder(coords_plane)
                    if coords_tube != None:
                        atomlist.append( (6, coords_tube) )
                # Thomson points
                pt = plane2cylinder(vij)
                if pt != None:
                    thomson_tube.append( (12, pt) )
        return atomlist, thomson_tube
    def createCap(self, thomson_tube, nrT):
        # create initial Thomson points on the cap
        thomson_cap = []
        # number of points in theta direction
        nr_th = int(np.floor(np.sqrt(nrT)))
        # number of points in phi direction
        nr_ph = nr_th
#        assert nr_th*nr_ph == nrT
        angles = []
        """
        for ph in np.linspace(0.001, 2.0*np.pi, nr_ph):
            for th in np.linspace(0.001, np.pi/2.0, nr_th):
                angles.append( (ph, th) )
        # distribute the other points randomly
        nr_rest = nrT - nr_th*nr_ph
        """
        nr_rest = nrT
        for i in range(0, nr_rest):
            th = np.random.rand(1)[0]*np.pi/2.0
            ph = np.random.rand(1)[0]*2*np.pi
            angles.append( (ph, th) )
        for (ph, th) in angles:
            x = self.R*np.sin(th)*np.cos(ph)
            y = self.R*np.sin(th)*np.sin(ph)
            z = self.L + self.R*np.cos(th)
            thomson_cap.append( (12, np.array([x,y,z])) )
        assert len(thomson_cap) == nrT
        thomson_cap_min = solve_constrained_Thomson(thomson_tube, thomson_cap, self.L, self.R)
        return thomson_cap
#        return thomson_tube + thomson_cap
        
def solve_constrained_Thomson(thomson_tube, thomson_cap, L, R):
    def potential_energy(thomson_pts):
        en = 0.0
        for i,(Zi,posi) in enumerate(thomson_pts):
            for j,(Zj,posj) in enumerate(thomson_pts):
                if i == j:
                    continue
                rij = la.norm(posi-posj)
                if rij == 0.0:
                    raise ValueError("rij = 0")
                en += 1.0/rij
        en *= 0.5
        return en
    def gradient_energy(thomson_pts):
        grad = []
        for k,(Zk,posk) in enumerate(thomson_pts):
            grad_k = np.zeros(3)
            for j,(Zj,posj) in enumerate(thomson_pts):
                if j == k:
                    continue
                rjk = posk-posj
                grad_k += rjk/pow(la.norm(rjk),3)
            grad.append( (Zk, grad_k) )
        return grad
    def enforce_constaints(thomson_cap):
        thomson_cap_rescaled = []
        for i,(Zi, posi) in enumerate(thomson_cap):
            x,y,z = posi
            z -= L
            if x >= 0:
                # keep point glued to a sphere of radius R
                r = np.sqrt(x*x+y*y+z*z)
                sc = R/r
                x,y,z = sc*x,sc*y,sc*z+L
            else:
                # point has moved into the tube
                # keep it stuck to a cylinder of radius R
                r = np.sqrt(x*x+y*y)
                sc = R/r
                x,y = sc*x,sc*y
                # z is left as it is
            thomson_cap_rescaled.append( (Zi,np.array([x,y,z])) )
        return thomson_cap_rescaled
    # function to minimize
    thomson_cap_ref = thomson_cap
    XYZ.write_xyz("/tmp/thomson_minimization.xyz", [thomson_cap_ref], mode="w")
    def f(x):
        thomson_cap = XYZ.vector2atomlist(x, thomson_cap_ref)
#        thomson_cap = enforce_constaints(thomson_cap)
        print "thomson_cap"
        print thomson_cap
        thomson_pts = thomson_tube + thomson_cap
        XYZ.write_xyz("/tmp/thomson_minimization.xyz", [thomson_cap], mode="a")
        en = potential_energy(thomson_pts)
        # gradient of potential energy
        grad = XYZ.atomlist2vector(gradient_energy(thomson_pts))
        print "en = %s" % en
        return en #, grad
    def grad(x):
        thomson_cap = XYZ.vector2atomlist(x, thomson_cap_ref)
        thomson_pts = thomson_tube + thomson_cap
        gr = XYZ.atomlist2vector(gradient_energy(thomson_pts))
        return gr
    x0 = XYZ.atomlist2vector(thomson_cap)
    err = optimize.check_grad(f, grad, x0)
    assert err < 1.0e-10, "err = %s" % err
    xmin,fmin,d = optimize.fmin_l_bfgs_b(f, x0)
    thomson_cap_min = XYZ.vector2atomlist(xmin, thomson_cap_ref)
    return thomson_cap_min

def dual_lattice(atomlist):
    Con = XYZ.connectivity_matrix(atomlist)
    # place and atom between all closest 
    pass

if __name__ == "__main__":
    cnt = CarbonNanotube(10,5)
    atomlist_tube, thomson_tube = cnt.createTube()
    XYZ.write_xyz("/tmp/cnt.xyz", [atomlist_tube])
    XYZ.write_xyz("/tmp/thomson.xyz", [thomson_tube])
    thomson_cap = cnt.createCap(thomson_tube, 3)
    XYZ.write_xyz("/tmp/thomson_cap.xyz", [thomson_cap])
