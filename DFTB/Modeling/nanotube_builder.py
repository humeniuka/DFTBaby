#!/usr/bin/env python
"""
This program helps in constructing capped carbon nanotubes

This implementation is based on the ideas in the following paper:

'Generalized method for constructing the atomic coordinates of nanotube caps'
by Robinson,M.; Suarez-Martinez,I.; Marks,N.   Phys. Rev. B 87, 155430 (2013)
"""
from os.path import join, expandvars, expanduser

import numpy as np
import numpy.linalg as la
from scipy import optimize

from DFTB import XYZ, AtomicData
from DFTB.extensions import thomson
from DFTB.Formats.Gaussian2py import Gaussian, Checkpoint

class CarbonNanotube:
    def __init__(self, n1, n2, L, a=4.6487, verbose=True):
        """
        Parameters:
        ===========
        n1, n2: chiral indeces of the nanotube
        L: length of the nanotube in bohr
        a: lattice constant in bohr
        """
        self.n1, self.n2 = n1,n2
        self.a = a
        self.L = L
        self.chiral_vector, self.tube_vector, self.a1, self.a2, self.c1, self.c2 = chiral_tube_vectors(self.n1,self.n2, self.L, self.a)
        # area of unit cell
        self.Aunit = la.norm(np.cross(self.a1,self.a2))
        circ = la.norm(self.chiral_vector) # circumferance
        self.R = circ/(2.0*np.pi) # radius of nanotube
        if verbose == True:
            print ""
            print "Carbon-Nanotube with Cap"
            print "========================"
            print "Chirality: (%s,%s)" % (self.n1,self.n2)
            print "Length:     %7.4f bohr (%7.4f Ang)" % (self.L, self.L*AtomicData.bohr_to_angs)
            print "Diameter:   %7.4f bohr (%7.4f Ang)" % (2*self.R, 2*self.R*AtomicData.bohr_to_angs)
            print ""
    def NrCapAtoms(self):
        """
        This method determines the number of carbon atoms that should
        be placed on a spherical cap assuming that the atoms density
        is the same as in graphene.
        """
        # area of hemisphere
        Ahemi = 0.5 * 4.0 * np.pi * pow(self.R,2) 
        # Determine how many Thomson points should be placed on the hemisphere
        # number of carbon atoms on the hemisphere, there are 2 atoms per unit cell
        nC = Ahemi / self.Aunit * 2
        # number of Thomson points
        nT = int((nC + 4)/2.0)  # formula (3) of Robinson (2013)
        return nC, nT

class ThomsonProblem:
    def __init__(self, cnt, Ncap, north_south="N"): #18.897): #4*18.897):
        """
        Parameters:
        ===========
        cnt: instance of CarbonNanotube
        Ncap: number of Thomson points that should be placed on the cap
        north_south: 'N' or 'S' depending on which side of the tube should be capped
        """
        self.cnt = cnt
        self.Ncap = Ncap
        self.north_south = north_south
    def addCap(self):
        # 
        x_cap = initial_distribution_cap(self.Ncap) # positions of Thomson points on the cap
        # for points in the cap x[i,:] = (theta_i, phi_i)
        # and for points in the tube x[i,:] = (z_i, phi_i)
        s_cap = ["C" for i in range(0, self.Ncap)] # segment of carbo-nanotube
        # to which the Thomson point belongs, 'C' for cap and 'T' for tube
        o_cap = np.ones(self.Ncap, dtype=int)
        # o[i] == 1, position of i-th Thomson point is optimized
        # o[i] == 0, position is held constant
        #
        x_tube, o_tube, s_tube = thomson_tube(self.cnt.chiral_vector, 
                                  self.cnt.tube_vector, self.cnt.a1, self.cnt.a2)

        # These unit vectors define the orientation of the
        # coordinate system
        ex = np.array([1,0,0])
        ey = np.array([0,1,0]) 
        ez = np.array([0,0,1]) 
        self.axes = [ex,ey,ez]    
        self.Z0 = 0.0
        
        # combine Thomson point of cap with those of the tube
        self.x = np.vstack((x_cap, x_tube))
        self.s = np.array(s_cap + s_tube)
        self.o = np.hstack((o_cap, o_tube))
    
        self.N = len(self.x)

        if self.north_south == "N":
            pass
        elif self.north_south == "S":
            self.x = reflect_tube_xy(self.x, self.s, self.cnt.L)

        atomlist = self.thomson2atomlist()
        XYZ.write_xyz("/tmp/thomson_opt.xyz", [atomlist], mode="w")

        active, inactive = active_layers(self.x,self.s,self.cnt.a1,self.cnt.a2)

        def savexyz():
            atomlist = self.thomson2atomlist()
            XYZ.write_xyz("/tmp/thomson_opt.xyz", [atomlist], mode="a")

        self.x[active,:], self.s[active], self.o[active], self.en, self.grad = optimize_thomson(self.x[active,:], self.s[active], self.o[active], R=self.cnt.R, callback=savexyz)
        return (self.en, self.grad)
    def dual_lattice(self):
        """
        construct the dual hexagonal lattice by placing a carbon
        atom in the center of each triangle and scaling it to lie
        on the hemisphere or tube
        """
        R = self.cnt.R
        triangles = self.triangulate()
        atomlist = []
        for (i,j,k) in triangles:
            posi = self.getCoords(i)
            posj = self.getCoords(j)
            posk = self.getCoords(k)
            # center of triangle
            center = (posi+posj+posk)/3.0
            if center[2] > 0.0:
                # scale to hemisphere
                center *= R / la.norm(center)
            else:
                # scale to tube
                x,y,z = center
                r = np.sqrt(x*x+y*y)
                x *= R/r
                y *= R/r
                center = np.array([x,y,z])
            #
            atomlist.append( (6, center) )
        if self.north_south == "S":
            atomlist = reflect_atomlist_xy(atomlist, self.cnt.L)
        return atomlist
    def triangulate(self):
        # for each vertex determine the atoms that are closer than 
        # (1+eps)*(lattice constant a)
        r2,dr2d1,dr2d2 = thomson.thomson.rij2(self.x,self.s,self.o, self.cnt.R, self.Z0, 0)
#        r2,dr2d1,dr2d2 = compute_rij2(self.x,self.s,self.o, R=self.R, opt=False)
        a2 = self.cnt.a**2
        r2min = a2 * 0.5
        r2max = a2 * 1.5
        # build connectivity matrix
        Con = np.zeros((self.N, self.N), dtype=int)
        for i in range(0, self.N):
            for j in range(i+1, self.N):
                if r2min < r2[i,j] < r2max:
                    # There is a bond between i and j
                    Con[i,j] = 1
        # 
#        print "Connectivity"
#        print Con
        faces = []
        for i in range(0, self.N):
            # select 
            bonded_i = np.where(Con[i,:] == 1)[0]
            for j in bonded_i:
                for k in bonded_i:
                    if Con[j,k] == 1:
                        # bonds i-j  i-k and j-k
                        face = set((i,j,k))
                        if not (face in faces):
                            faces.append(face)
#        print "Faces"
#        print faces
        return faces
    def getCoords(self, i):
        """
        This finds the cartesian coordinates of the i-th Thomson point.
        If the point lies on the cap, one has to transform from spherical to
        cartesian coordinates, if it lies on the tube from cylindrical to
        cartesian. 
        """
        R = self.cnt.R
        X = self.x
        if self.s[i] == "C":
            x = R*np.sin(X[i,0])*np.cos(X[i,1])
            y = R*np.sin(X[i,0])*np.sin(X[i,1])
            z = R*np.cos(X[i,0]) + self.Z0
        elif self.s[i] == "T":
            x = R*np.cos(X[i,1])
            y = R*np.sin(X[i,1])
            z = X[i,0]
        else:
            raise Exception("????")
        ex,ey,ez = self.axes
        pos = x*ex + y*ey + z*ez
        return pos
    def thomson2atomlist(self):
        """
        Find the coordinates of the Thoms points and save them as they belonged to a molecule.
        To be able to visualize the points with molden, the points on the cap are encoded as 
        Ca atoms, and the points on the tube as Sr atoms. The points that are excluded from
        the optimization are represented as Mg atoms.
        """
        atomlist = []
        for i in range(0, self.N):
            pos = self.getCoords(i)
            if self.o[i] == 1:
                if self.s[i] == "C":
                    atnum = 20 # Z(Ca) = 20
                else:
                    atnum = 38 # Z(Sr) = 38
            else:
                atnum = 12 # Z(Mg) = 12
            atomlist.append( (atnum, pos) )
        return atomlist

def initial_distribution_cap(n):
    """
    Initially the n Thomson points are randomly distributed on the cap. This means that one might 
    have to run several optimizations with different initial distributions to find the best cap
    (one without adjacent pentagons) and that some initial distributions won't work.
    """
    nr_th = int(np.floor(np.sqrt(n)))
    # number of points in phi direction
    nr_ph = nr_th
    
    x = np.zeros((n,2))
    i=0
    """
    for ph in np.linspace(0.0, 2.0*np.pi*(1-1.0/float(nr_ph)), nr_ph):
        for th in np.linspace(np.pi/2.0*(0.5/float(nr_th)), np.pi/2.0*(1.0-0.5/float(nr_th)), nr_th):
            x[i,0] = th
            x[i,1] = ph
            i+=1
    # distribute the other points randomly
    nr_rest = n - nr_th*nr_ph
    """
    #
    nr_rest = n
    #
    for j in range(0, nr_rest):
        th = np.random.rand(1)[0]*np.pi/2.0
        ph = np.random.rand(1)[0]*2*np.pi
        x[i,0] = th
        x[i,1] = ph
        i+=1
    assert i == n
    return x

def chiral_tube_vectors(n1,n2, L, a):
    """
    For a set of chiral indeces (n1,n2) and
    the length of the tube L find the chiral vector and
    the vector perpendicular to it, which points along the tube axis

    Parameters:
    ===========
    n1, n2: chiral indeces
    L: length of the tube
    a: lattice constant
    """
    # primitive lattice vectors 
    # This lattice is triangular 
    a1 = a*np.array([np.sqrt(3.0)/2.0, -1.0/2.0, 0.0])
    a2 = a*np.array([np.sqrt(3.0)/2.0,  1.0/2.0, 0.0])
#    # position of carbon atoms
    c1 = (a1 + a2)/3.0
    c2 = (a1 + a2)*2.0/3.0
    # reciprocal lattice vectors satisfy ai*bj = 2 pi * delta_ij
    b1 = 2*np.pi/a * np.array([1.0/np.sqrt(3.0),-1.0, 0.0])
    b2 = 2*np.pi/a * np.array([1.0/np.sqrt(3.0), 1.0, 0.0])
    # chiral vector => circumference
    chiral_vector = n1*a1 + n2*a2
#    # chiral angle eqn. 2.1 from Reich's book
#    theta = np.arccos( (n1+0.5*n2)/np.sqrt(n1**2+n1*n2+n2**2) )
    # tube radius
    R = la.norm(chiral_vector)/(2*np.pi)
    # the vector orthogonal to the chiral vector
    if n1 != 0:
        m2 = 1.0
        m1 = - n2/float(n1) * m2
    else:
        m2 = 0.0
        # m1 can be anything
        m1 = 1.0
    tube_vector = m1*b1 + m2*b2
    # normalze tube axis to L
    tube_axis = tube_vector / la.norm(tube_vector)
    tube_vector = L*tube_axis

    return chiral_vector, tube_vector, a1, a2, c1, c2

def thomson_tube(chiral_vector, tube_vector, a1, a2):
    chiral_norm = la.norm(chiral_vector)
    chiral_axis = chiral_vector / chiral_norm
    L = la.norm(tube_vector)
    tube_axis = tube_vector / L
    x = []  # coordinates of Thomson points on the tube (z,phi)
    o = []  # flag (0 if thomson points are to be fixed in place or
            # 1 if they are allowed to move)
    s = []  # Type of the Thomson point:
            #   'T': point lies on the tube, the coordinate pair is interpreted a cylindrical coordinates (z, phi)
            #   'C': point lies on the cap, the coordiante pair is interpreted as (theta, phi)
    #
    da = min(la.norm(a1), la.norm(a2))
    ncirc = int(chiral_norm/da)
    nL = int(L/da)
    nmax = 2*max(ncirc,nL)
    for i in range(-nmax,nmax):
        for j in range(-nmax, nmax):
            # lattice vector
            vij = i*a1 + j*a2
            # The normalized coordinates (u,v) map the square [0,1]x[0,1]
            # to the surface of the cylinder
            u = np.dot(vij, chiral_axis)/chiral_norm
            v = np.dot(vij, tube_axis)/L
#            if (0 <= u < 1) and (0 <= v <= 1):
            if (-da/L <= u < 1+da/L) and (0 <= v <= 1):
                # Only the points that will lie on the cylinder of length
                # L are added
                z = -v*L
                phi = 2.0*np.pi*u
                x.append([z,phi])
                o.append(0)
                s.append("T")

    x = np.array(x)
    o = np.array(o)
    return x, o, s

def reflect_tube_xy(x, s, L):
    """
    reflect thomson points on the tube at the xyz plane
    """
    x_refl = np.copy(x) # reflected points
    n, dummy = x.shape
    for i in range(n):
        if s[i] == "T":
            z,phi = x[i,:]
            x_refl[i,:] = np.array([-z-L,phi])
    return x_refl

def active_layers(x,s,a1,a2):
    """
    Split Thomson points into a set containing the cap + the first layer of the tube
    and the remaining tube. For optimizing the cap shape only the repulsion from
    the closest layer is important. Neglecting the remaining layers speeds up the
    calculation
    """
    active = []    # indeces of active Thomson points (cap + first layers)
    inactive = []  # indeces of inactive Thomson points (rest of the tube)
    dz_layer = 10.0 #2*la.norm(a1+a2) # approximate width of one layer
    n,dummy = x.shape
    for i in range(n):
        if s[i] == "T":
            z,phi = x[i,:]
            if z > -dz_layer:
                active.append(i)
            else:
                inactive.append(i)
        else:
            # all caps states should be active
            active.append(i)
    return active, inactive
    

def reflect_atomlist_xy(atomlist, L):
    atomlist_refl = []
    for (Zi, posi) in atomlist:
        x,y,z = posi
        atomlist_refl.append( (Zi, np.array([x,y,-z-L])) )
    return atomlist_refl

def delta(i,j):
    if i == j:
        return 1.0
    else:
        return 0.0

def compute_rij2(x,s,o, R=1.0, Z0=0.0, opt=True):
    n,dummy = x.shape
    # distance^2
    r2 = np.zeros((n,n))
    for i in range(0, n):
        for j in range(0, n):
            if i==j:
                continue
            if o[i] == 0 and o[j] == 0 and opt==True:
                # energy between fixed points is not interesting
                r2[i,j] = 100000000000000.0
                continue
            if s[i] == "T" and s[j] == "T":
                # Tube-Tube
                r2[i,j] = 2*R**2*(1.0-np.cos(x[i,1]-x[j,1])) + (x[i,0]-x[j,0])**2
            elif s[i] == "T" and s[j] == "C":
                # Tube-Cap
                r2[i,j] = R**2*(1+np.sin(x[j,0])**2 - 2*np.sin(x[j,0])*np.cos(x[i,1]-x[j,1])) + (x[i,0]-Z0-R*np.cos(x[j,0]))**2
            elif s[i] == "C" and s[j] == "T":
                # Cap-Tube
                r2[i,j] = R**2*(1+np.sin(x[i,0])**2 - 2*np.sin(x[i,0])*np.cos(x[j,1]-x[i,1])) + (x[j,0]-Z0-R*np.cos(x[i,0]))**2
            elif s[i] == "C" and s[j] == "C":
                # Cap-Cap
                r2[i,j] = 2*R**2*(1-(np.sin(x[i,0])*np.sin(x[j,0])*np.cos(x[i,1]-x[j,1]) + np.cos(x[i,0])*np.cos(x[j,0])))
            else:
                raise Exception("BUG: s[%s] = %s   s[%s] = %s ?" % (i, j, s[i], s[j]))
            assert r2[i,j] > 0.0

    # gradients
    dr2d1 = np.zeros((n,n,n)) # with respect to theta or z
    dr2d2 = np.zeros((n,n,n)) # with respect to phi

    for i in range(0, n):
        for j in range(0, n):
            for a in [i, j]:
                if i == j:
                    continue
                if o[i] == 0 and o[j] == 0:
                    # energy between fixed points remains constant
                    continue               
                if s[i] == "T" and s[j] == "T":
                    # Tube-Tube
                    # derivative with respect to z
                    dr2d1[i,j,a] = 2*R**2*np.sin(x[i,1]-x[j,1])*(delta(a,i)-delta(a,j))
                    # derivative with respect to phi
                    dr2d2[i,j,a] = 2*(x[i,0]-x[j,0])*(delta(a,i)-delta(a,j))
                elif s[i] == "C" and s[j] == "C":
                    # Cap-Cap
                    # derivative with respect to theta
                    dr2d1[i,j,a] = 2*R**2*(np.sin(x[i,0])*np.cos(x[j,0])-np.cos(x[i,0])*np.sin(x[j,0])*np.cos(x[i,1]-x[j,1]))*delta(a,i) \
                                  +2*R**2*(np.cos(x[i,0])*np.sin(x[j,0])-np.sin(x[i,0])*np.cos(x[j,0])*np.cos(x[i,1]-x[j,1]))*delta(a,j)
                    # derivative with respect to phi
                    dr2d2[i,j,a] = 2*R**2*np.sin(x[i,0])*np.sin(x[j,0])*np.sin(x[i,1]-x[j,1])*(delta(a,i)-delta(a,j))
                elif s[i] == "T" and s[j] == "C":
                    # Tube-Cap
                    # derivative with respect to z
                    dr2d1[i,j,a] += 2*(x[i,0]-Z0-R*np.cos(x[j,0]))*delta(a,i)
                    # derivative with respect to theta
                    dr2d1[i,j,a] += 2*R*(-R*np.cos(x[j,0])*np.cos(x[i,1]-x[j,1]) + (x[i,0]-Z0)*np.sin(x[j,0]) )*delta(a,j)
                    # derivative with respect to phi
                    dr2d2[i,j,a] += 2*R**2*np.sin(x[j,0])*np.sin(x[i,1]-x[j,1])*(delta(a,i)-delta(a,j))
                elif s[i] == "C" and s[j] == "T":
                    # Tube-Cap
                    # derivative with respect to z
                    dr2d1[i,j,a] += 2*(x[j,0]-Z0-R*np.cos(x[i,0]))*delta(a,j)
                    # derivative with respect to theta
                    dr2d1[i,j,a] += 2*R*(-R*np.cos(x[i,0])*np.cos(x[j,1]-x[i,1]) + (x[j,0]-Z0)*np.sin(x[i,0]) )*delta(a,i)
                    # derivative with respect to phi
                    dr2d2[i,j,a] = 2*R**2*np.sin(x[i,0])*np.sin(x[j,1]-x[i,1])*(delta(a,j)-delta(a,i))
                else:
                    raise Exception("s[i] has to be C or T")
    return r2,dr2d1,dr2d2
        

def thomson_energy(x,s,o, R=1.0, Z0=0):
    n,dummy = x.shape
    r2,dr2d1,dr2d2 = thomson.thomson.rij2(x,s,o, R, Z0, 1)
# check against python implementation
#    r2_py,dr2d1_py,dr2d2_py = compute_rij2(x,s,o, R=R, Z0=Z0)
#    err_r2 = np.sum(abs(r2[o==1]-r2_py[o==1]))
#    assert err_r2 < 1.0e-8, "err(r2) = %s" % err_r2
#    err_dr2d1 = np.sum(abs(dr2d1-dr2d1_py))
#    assert err_dr2d1 < 1.0e-8, "err(dr2d1) = %s" % err_dr2d1
#    err_dr2d2 = np.sum(abs(dr2d2-dr2d2_py))
#    assert err_dr2d2 < 1.0e-8, "err(dr2d2) = %s" % err_dr2d2

    coul_f, grad_f = thomson.thomson.coulomb_energy(o,r2,dr2d1,dr2d2)
# compare with python implementation    
    """
    # build Coulomb energy
    coul = 0.0
    grad = np.zeros((n,2))
    for i in range(0, n):
        for j in range(0, n):
            if i == j:
                continue
            coul += 0.5*pow(r2[i,j], -0.5)
            if r2[i,j] == 0.0:
                print "i = %s  j = %s" % (i,j)
            # gradient
            for a in range(0, n):
                if o[a] == 0:
                    continue
                grad[a,0] -= 1.0/4.0 * pow(r2[i,j],-1.5) * dr2d1[i,j,a]
                grad[a,1] -= 1.0/4.0 * pow(r2[i,j],-1.5) * dr2d2[i,j,a]
    print "coul = %s" % coul
#    print "grad = %s" % grad
    print "|grad| = %s" % la.norm(grad.ravel())

    # compare with Fortran
    err_coul = abs(coul - coul_f)
    err_grad = np.sum(abs(grad-grad_f))
    assert err_coul < 1.0e-8, "err(coul) = %s" % err_coul
    assert err_grad < 1.0e-8, "err(grad) = %s" % err_grad
    """
#    print "|grad| = %s" % la.norm(grad_f.ravel())
    return coul_f, grad_f


def optimize_thomson(x, s, o, R=1.0, callback=None):
    """
    Find the arrangement of the Thomson points that maximizes the distance between them.
    This is equivalent to minimizing the electrostatic energy of charged particles sitting
    at the Thomson points:
    
    E = sum_ij 1/(ri-rj)
    
    During the optimization Thomson points can move between the cap
    and the tube. When a point switches from the cap to the tube
    its type is changed and from there on its motion is constrained to 
    the cylinder. On the other hand, a point that moves onto the cap
    is constrained to a sphere. 

    This function adjusts the types (in the array self.s) of the points and transforms their
    coordinates (in the array self.x) from cylindrical to spherical coordinates at the boundary
    between the cap and the tube. Thomson points that are held constant (self.o[i] == 0)
    are excluded.
    """
#    print "OPTIMIZE"
    N = len(x)
    # target function to minimize
    def f(xo):
        x[o == 1,:] = np.reshape(xo, (len(xo)/2,2))
            # 
        en,gr = thomson_energy(x, s, o, R)
        grad = gr[o == 1,:].flatten()
        return en, grad
    # move thomson points between tube and cap if necessary
    def adjust_assignment(xo):
        x[o == 1,:] = np.reshape(xo, (len(xo)/2,2))
        for i in range(0, N):
            if o[i] == 1:
                if s[i] == "C" and x[i,0] > np.pi/2:
                    # point has moved into tube
#                    print "C => T"
                    s[i] = "T"
                    x[i,0] = 0.0 #self.R * np.cos(self.x[i,0]) #0.0 # set z_i = -Z0
                elif s[i] == "T" and x[i,0] > 0.0:
                    # point has moved back onto sphere
#                    print "T => C"
                    s[i] = "C"
                    x[i,0] = np.pi/2.0 #np.arccos(self.x[i,0]/self.R) #np.pi/2.0
        if callback != None:
            callback()
    # 
    x0 = x[o == 1,:]
    #  # Test gradient
    # err = optimize.check_grad(lambda x: f(x)[0], lambda x: f(x)[1], x0.flatten())
    # assert err < 1.0e-4, "err = %s" % err
    xmin, fmin, d = optimize.fmin_l_bfgs_b(f, x0.flatten(), \
          callback=lambda x: adjust_assignment(x), \
          approx_grad=0, \
          factr=0.01)
    x[o == 1,:] = np.reshape(xmin, (len(xmin)/2,2))
    return x,s,o, fmin, d["grad"]

def remove_overlapping_atoms(atomlist, tol=1.0e-4):
    print "remove overlapping atoms"
    nat = len(atomlist)
    # The i-th atom should be removed if duplicate[i] == 1
    # Fortran
    rs = XYZ.atomlist2vector(atomlist)
    rs = np.reshape(rs, (nat, 3))
    duplicate = thomson.thomson.remove_overlapping(rs, tol)
    # python
    """
    duplicate_py = np.zeros(nat, dtype=int)
    for i in range(0, nat):
        for j in range(i+1,nat):
            posi = np.array(atomlist[i][1])
            posj = np.array(atomlist[j][1])
            Rij = la.norm(posi-posj)
            if Rij < tol:
                duplicate_py[j] = 1
    assert np.sum(abs(duplicate - duplicate_py)) < 1.0e-10
    """
    # only copy unique atoms
    atomlist_unique = []
    for i in range(0, nat):
        if duplicate[i] == 0:
            atomlist_unique.append( atomlist[i] )
    return atomlist_unique

def check_connectivity(atomlist):
    """
    In a capped carbon nanotube, every atoms should be connected
    to 6 or 5 neighbours. If this is not the case, there must a
    hole somewhere.
    """
    Con = XYZ.connectivity_matrix(atomlist)
    nr_neighbours = np.sum(Con, axis=1)
    con_min = nr_neighbours.min()
    con_max = nr_neighbours.max()
    print "nr_neighbours = %s" % nr_neighbours
    print "%s <= connectivity <= %s" % (con_min, con_max)
    if (3 <= con_min) and (con_max <= 3):
        return True
    else:
        print "Not all C-atoms are connected to 3. There must be a hole somewhere!"
        return False
        
###############################################################


################################################################################

def find_best_cap(cnt, Ncap, Ntrial=10, north_south="N"):
    """
    This function creates Ntrial caps starting with different randomly distributed
    points on a sphere and select the cap with the lowest Coulomb
    energy.

    Parameters:
    ===========
    cnt: instance of class CarbonNanotube
    """
    caps = []
    ens = []
    print "CAP FOR %s-HEMISPHERE" % north_south
    for n in range(0, Ntrial):
        tp = ThomsonProblem(cnt, Ncap, north_south=north_south)
        en,grad = tp.addCap()
        grad_nrm = la.norm(grad)
        if grad_nrm < 1.0e-4:
            atomlist = tp.dual_lattice()
            caps.append(atomlist)
            ens.append(en)
            print "Trial %s (%s points)   => cap with Coulomb en = %s" % (n,Ncap,en)
            break
        else:
            print "Trial %s (%s points)   optimization did not converge, |grad| = %s" % (n,Ncap,grad_nrm)
    nmin = np.argmin(ens)
    print "I select cap %s" % nmin
    return caps[nmin]

def capped_nanotube(cnt,NcapN, NcapS):
    import time
    np.random.seed(int(1000*time.time()))
    # north pole
    atomlistN = find_best_cap(cnt,NcapN,Ntrial=100, north_south="N")
    # south pole
    atomlistS = find_best_cap(cnt,NcapS,Ntrial=100, north_south="S")
    # combine the two lattices and remove overlapping atoms, actually
    # only the atoms on the tube should overlap
    atomlist_carbons = remove_overlapping_atoms(atomlistN + atomlistS)
    return atomlist_carbons

def capped_uff_nanotube(cnt,NcapN=0, NcapS=0, optimize_uff=1, out_xyz="/tmp/cnt.xyz"):
    """
    build capped nanotube and optimize it with UFF
    """
    nC,nT = cnt.NrCapAtoms()
    print "Expected number of Thomson points nT = %s" % nT
    print "If the caps have holes or very close atoms, you should try increasing"
    print "or reducing the number of Thomson points (by changing the options NcapN or NcapS)"
    if NcapN == 0:
        NcapN = nT
    if NcapS == 0:
        NcapS = nT

    """
    Different number of Thomson points are tested to find
    a pair (NcapN,NcapS) 
    """
    dcap = []
    dnmin=0
    dnmax=1
    ntrial=5
    for dnN in range(dnmin, dnmax):
        for dnS in range(dnmin, dnmax):
            # try for 5 times
            for n in range(0, ntrial):
                dcap.append( (dnN,dnS) )
    dcap = sorted(dcap, key=lambda (a,b): abs(a)+abs(b))

    for (dnN,dnS) in dcap:
        print "north cap: %s points, south cap: %s points" % (NcapN+dnN,NcapS+dnS)
        atomlist_carbons = capped_nanotube(cnt,NcapN+dnN,NcapS+dnS)
        if optimize_uff == 1:
            print "CNT will be optimized further with the Universal Force Field of G09"
            tmp_dir="/tmp"
            com_file = join(tmp_dir, "cnt_uff_opt.com")
            chk_file = join(tmp_dir, "cnt_uff_opt.chk")
            fchk_file = chk_file.replace(".chk", ".fchk")
            Gaussian.write_input(com_file, atomlist_carbons, \
                route="# UFF Opt", \
                title="Optimize (%s,%s)-Nanotube with UFF" % (cnt.n1,cnt.n2), chk=chk_file,mem=1000)
            try:
                Gaussian.run(com_file)
                # format checkpoint file
                Gaussian.formchk(chk_file, fchk_file)
                Data = Checkpoint.parseCheckpointFile(fchk_file)
                pos = Data["_Current_cartesian_coordinates"]
                atomlist_uff = XYZ.vector2atomlist(pos, atomlist_carbons)
                atomlist_carbons = atomlist_uff
            except Gaussian.GaussianError:
                print "UFF-optimization failed!"
                continue
        if check_connectivity(atomlist_carbons) == True:
            # no holes, correct connectivity
            break
        XYZ.write_xyz("/tmp/trials.xyz", [atomlist_carbons], mode="a")
    else:
        print ""

    XYZ.write_xyz(expandvars(expanduser(out_xyz)), [atomlist_carbons])
    print "Geometry of capped CNT written to %s" % out_xyz


if __name__ == "__main__":
    import sys
    import optparse

    usage = "Usage: python %s\n" % sys.argv[0]
    usage += "  builds the geometry for a capped single-walled carbon-nanotube (CNT)\n\n"

    parser = optparse.OptionParser(usage)
    parser.add_option("--n1", dest="n1", help="Chiral index n1 [default: %default]", default=6, type=int)
    parser.add_option("--n2", dest="n2", help="Chiral index n2 [default: %default]", default=5, type=int)
    parser.add_option("--L", dest="L", help="Length of CNT in bohr [default: %default]", default=100.0, type=float)
    parser.add_option("--NcapN", dest="NcapN", help="Number of Thomson points on the northern cap, if set to 0 the number is estimated from the density of carbon atoms in graphene [default: %default]", default=0, type=int)
    parser.add_option("--NcapS", dest="NcapS", help="Number of Thomson points on the southern cap, if set to 0 the number is estimated from the density of carbon atoms in graphene [default: %default]", default=0, type=int)
    parser.add_option("--out_xyz", dest="out_xyz", help="Save the CNT geometry to this xyz-file [default: %default]", default="/tmp/cnt.xyz") 
    parser.add_option("--optimize_uff", dest="optimize_uff", help="Optimize the built CNT with the UFF force field. This requires that the program g09 is installed. [default: %default]", default=1, type=int)

    (opts,args) = parser.parse_args()

    cnt = CarbonNanotube(opts.n1, opts.n2, opts.L)
    capped_uff_nanotube(cnt,
                        NcapN=opts.NcapN, NcapS=opts.NcapS, 
                        optimize_uff=opts.optimize_uff, out_xyz=opts.out_xyz)
