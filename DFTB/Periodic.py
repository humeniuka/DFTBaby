import numpy as np
import numpy.linalg as la
import itertools

class CrystalSystem:
    dim = 3
    def setLatticeVectors(self, vecs):
        # primitive lattice vectors of space lattice a[i,:]
        self.a = np.vstack(vecs)
        # primitive reciprocal lattice vectors b[i,:]
        self.b = 2.0 * np.pi * la.inv(self.a.transpose())
        # volume of real (Wigner-Seitz) and reciprocal (1st Brillouin zone) unit cell
        self.Va = abs(la.det(self.a))
        self.Vb = abs(la.det(self.b))
    def getLatticeVectors(self, space="real"):
        if space == "real":
            return list(self.a)
        else:
            return list(self.b)
    def setBasis(self, atomlist):
        self.basis = atomlist
    def T(self, n):
        # lattice vector
        t = np.zeros(self.dim)
        for i,ni in enumerate(n):
            t += ni*self.a[i,:]
        return t
    def translation_vectors(self, nmax, space="real", include_origin=True):
        #klm = list(itertools.product(range(-nmax+1, nmax), repeat=self.dim))
        klm = list(itertools.product(range(-nmax[0]+1, nmax[0]),range(-nmax[1]+1,nmax[1]), range(-nmax[2]+1,nmax[2])))
        Ts = []
        for int_disp in klm:
            if include_origin == False:
                # exlude the vector t=(0,0,0) 
                # which does not correspond to any translation
                if all(np.array(int_disp) == 0):
                    continue
            if space == "real":
                t = self.T(int_disp)
            else:
                t = self.G(int_disp)
            print "%s   ->  %s" % (str(int_disp), t)
            Ts.append(t)
        print "NR translation vectors = %s" % len(Ts)
        return Ts
    def G(self, m):
        # reciprocal lattice vector
        G = np.zeros(self.dim)
        for i,mi in enumerate(m):
            G += mi*self.b[i,:]
        return G
    def nearest_neighbours_r(self, Rmax):
        Nmax = [int(la.norm(self.b[i,:])/(2.0*np.pi)) for i in range(0, self.dim)]
        ns = list(itertools.product(*[range(-Nmax[i], Nmax[i]+1) for i in range(0, self.dim)]))
        print ns
        Rs = np.array([self.T(n) for n in ns])
        ds = np.array([la.norm(R) for R in Rs])
        sort_index = np.argsort(ds)
        ds = ds[sort_index]
        Rs = Rs[sort_index]
        shells = []
        print "Nearest neighbours"
        print "=================="
        for ishell, (k, g) in enumerate(itertools.groupby(Rs, key=lambda R: np.around(la.norm(R), decimals=8))):
            print "Shell %s" % ishell
            print "distance = %s" % k
            gl = list(g)
            for R in gl:
                print "    %s" % R
            shells.append(gl)
        return shells
    def nearest_neighbours_k(self, order):
        pass
    def getWignerSeitz(self, rs):
        # find the vertices of the Wigner-Seitz cell, this are the points
        # that are equidistant from 3 lattice points.
        def corner(R1,R2,R3):
            # 
            M = np.vstack( (R2-R1, R3-R1, R3-R2) )
            y = 0.5 * np.array([ \
                    np.dot(R1+R2, R2-R1), np.dot(R1+R3, R3-R1), np.dot(R2+R3,R3-R2)])
            x = la.lstsq(M, y, rcond=1.0e-10)[0] 
            return x
        # nearest neighbour shell
        nn = self.nearest_neighbours_r(pow(self.Va, 1.0/3.0)*2)[1]
        print nn
        R1 = np.array([0,0,0])
        Cs = []
        print "Corners of Wigner-Seitz polyhedron"
        print "=================================="
        for i2, R2 in enumerate(nn):
            for i3, R3 in enumerate(nn[i2+1:]):
                try:
                    c = corner(R1, R2, R3)
                    # make sure corner is closer to the origin than to any other lattice point
                    if 1.0e-10 < la.norm(c) <= min([la.norm(c - Rn) for Rn in nn]):
                        Cs.append(c)
                        print c
                except la.LinAlgError:
                    pass
        self.plotPolyhedron(Cs, nn)
    def points_inside_BZ(self, ks):
        # select those points ki in ks that lie inside the Brillouin zone
        # ki has to lie closer to the origin than to any other point in the reciprocal lattice
#        nnk = self.nearest_neighbours_k(pow(self.Va
        pass
    def plotPolyhedron(self, corners, nn):
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        pts = np.vstack(corners)
        ax.scatter(pts[:,0], pts[:,1], pts[:,2])
        # connect triangles by lines if they lie on a bisecting plane
        for i1,c1 in enumerate(corners):
            for i2,c2 in enumerate(corners):
                edge1 = c2-c1
                for neigh in nn:
                    if abs(np.dot(edge1, neigh) < 1.0e-10) and la.norm(neigh) > la.norm(c1-neigh) and la.norm(neigh) > la.norm(c2-neigh):
                        # edge is perpendicular to 0-nn
                        n = neigh
                        l = np.vstack([c1,c2])
#                        ax.plot(l[:,0],l[:,1],l[:,2], color="blue", lw=1)
                        break
                else:
                    pass
        plt.show()
    def getBrillouinZone(self):
        pass
    def sample_BZ(self, nmax):
        # Monkhorst-Pack
        q = nmax
        def u(i,r):
            return (2.0*r - q[i] - 1.0)/(2.0*q[i])
        kpoints = []
        for p in range(1, q[0]+1):
            for r in range(1, q[1]+1):
                for s in range(1, q[2]+1):
                    k_prs = u(0,p)*self.b[0,:] + u(1,r)*self.b[1,:] + u(2,s)*self.b[2,:]
                    kpoints.append(k_prs)
        print "NR kpoints = %s" % len(kpoints)
        print "lattice vectors:"
        print self.b
        return kpoints
    def integrate_over_BZ(self, func, nmax):
        """
        integrate a function f(k) over the Brillouin zone
        """
        Ts = self.translation_vectors(nmax, space="real")
        ks = self.sample_BZ(nmax)
        fk = func(Ts,ks)
        print "self.Vb = %s" % self.Vb
        avg = sum(fk) * self.Vb / float(len(ks))
        return avg
    def __str__(self):
        txt = ""
        txt += "Lattice Vectors:\n"
        txt += "================\n"
        for i in range(0, self.dim):
            x,y,z = list(self.a[i,:]/(2.0*np.pi))
            txt += "a%d: %.6f %.6f %.6f * 2.0*pi\n" % (i+1,x,y,z)
        return txt
    def symmetry_kpoints(self):
        return {}
    def uvw2kxkykz(self, symmetry_kpoints):
        b1,b2,b3 = self.getLatticeVectors(space="reciprocal")
        kpoints = {}
        for name,(u,v,w) in symmetry_kpoints.iteritems():
            k = u*b1 + v*b2 + w*b3
            kpoints[name] = k
        return kpoints


# 3D lattice system

class SimpleCubic(CrystalSystem):
    def __init__(self, a):
        a1 = a * np.array([1,0,0]) 
        a2 = a * np.array([0,1,0]) 
        a3 = a * np.array([0,0,1]) 
        self.setLatticeVectors([a1,a2,a3])

class SimpleHexagonal(CrystalSystem):
    def __init__(self, a, c):
        a1 = a * np.array([1,0,0]) 
        a2 = a * np.array([1.0/2.0, np.sqrt(3.0)/2.0, 0]) 
        a3 = a * np.array([0,0,c/a]) 
        self.setLatticeVectors([a1,a2,a3])
    
class FaceCenteredCubic(CrystalSystem):
    def __init__(self, a):
        a1 = a * np.array([0,       1.0/2.0, 1.0/2.0]) 
        a2 = a * np.array([1.0/2.0,    0.0,  1.0/2.0]) 
        a3 = a * np.array([1.0/2.0, 1.0/2.0,  0.0   ]) 
        self.setLatticeVectors([a1,a2,a3])
    def symmetry_kpoints(self):
        # coordinates of symmetry points in the basis [b1,b2,b3]
        sympoints = {
                   "W": (0.25,0.75,0.50),
                   "X": (0.00,0.50,0.50,),
                   "Gamma": (0.00,0.00,0.00),
                   "L": (0.50,0.50,0.50),
                   "U": (0.25,5.0/8.0, 5.0/8.0),
                   "K": (3.0/8.0, 3.0/4.0, 3.0/8.0)}
        return self.uvw2kxkykz(sympoints)

    
class BodyCenteredCubic(CrystalSystem):
    def __init__(self, a):
        a1 = a * np.array([-1.0/2.0,  1.0/2.0, 1.0/2.0]) 
        a2 = a * np.array([1.0/2.0,  -1.0/2.0, 1.0/2.0]) 
        a3 = a * np.array([1.0/2.0, 1.0/2.0,  -1.0/2.0]) 
        self.setLatticeVectors([a1,a2,a3])
        
# 2D Lattice systems are modelled as 3D lattices with a very long period
# along the 3rd lattice vector
class PorpheneTetragonal(CrystalSystem):
    def __init__(self, a):
        a1 = a * 1.0/np.sqrt(2.0) * np.array([1.0, 1.0, 0.0])    
        a2 = a * 1.0/np.sqrt(2.0) * np.array([1.0, -1.0, 0.0])    
        a3 = 1000.0 * a * np.array([0.0, 0.0, 1.0])
        self.setLatticeVectors([a1,a2,a3])
    def symmetry_kpoints(self):
        # coordinates of symmetry points in the basis [b1,b2,b3]
        sympoints = {
                   "Gamma": (0.00,0.00,0.00),
                   "Sigma": (0.25,0.00,0.00),
                   "M":     (0.50,0.00,0.00),
                   "K":     (0.50,0.50,0.00),
                   "Lambda":(0.25,0.25,0.00)}
        return self.uvw2kxkykz(sympoints)

# 1D lattice system
class PorpheneWire(CrystalSystem):
    def __init__(self, a):
        a1 = a * 1.0/np.sqrt(2.0) * np.array([1.0, 1.0, 0.0])    
        a2 = 1000.0 * a * np.array([1.0, -1.0, 0.0])    
        a3 = 1000.0 * a * np.array([0.0, 0.0, 1.0])
        self.setLatticeVectors([a1,a2,a3])
    def symmetry_kpoints(self):
        # coordinates of symmetry points in the basis [b1,b2,b3]
        sympoints = {
                   "Gamma": (0.00,0.00,0.00),
                   "Sigma": (0.25,0.00,0.00),
                   "M":     (0.50,0.00,0.00),
                   "E":     (1.00,0.00,0.00)
                   }
        return self.uvw2kxkykz(sympoints)


if __name__ == "__main__":
    lat = BodyCenteredCubic(1.0)
#    lat = FaceCenteredCubic(1.0)
    print lat
    lat.nearest_neighbours_r(1.0)
#    lat.getWignerSeitz(None)
    print lat.integrate_over_BZ(lambda t,k: np.ones(len(t)), 5)
