import numpy as np
import numpy.linalg as la
from scipy import optimize

from DFTB import XYZ

#################### CARTESIAN COORDINATES WITH CONSTRAINTS #######################
# objective function
def f_cart(x):
    x = rescale(x)
    n = x.shape[0]
    npts = n/3
    coul = 0.0
    for i in range(0, npts):
        for j in range(0, npts):
            if i==j:
                continue
            coul += 1.0/la.norm(x[3*i:3*(i+1)]-x[3*j:3*(j+1)])
    #
    coul *= 0.5
    print "coul = %s" % coul
    return coul

# gradient of objective function
def grad_cart(x):
    n = x.shape[0]
    npts = n/3
    gr = np.zeros(n)
    for k in range(0, npts):
        for j in range(0, npts):
            if j == k:
                continue
            rjk = x[3*k:3*(k+1)]-x[3*j:3*(j+1)]
            gr[3*k:3*(k+1)] -= rjk/pow(la.norm(rjk),3)
    return gr

def rescale(x):
    # scale position vectors to unit sphere
    n = x.shape[0]
    npts = n/3
    for i in range(0, npts):
        x[3*i:3*(i+1)] /= la.norm(x[3*i:3*(i+1)])
    return x

def initial_thomson_points_sph(nrT):
    # number of points in theta direction
    nr_th = int(np.floor(np.sqrt(nrT)))
    # number of points in phi direction
    nr_ph = nr_th

    angles = []
    for ph in np.linspace(0.12, 2.0*np.pi-0.11, nr_ph):
        for th in np.linspace(0.13, np.pi-0.14, nr_th):
            angles += [ph, th]
    # distribute the other points randomly
    nr_rest = nrT - nr_th*nr_ph
    for i in range(0, nr_rest):
        th = np.random.rand(1)[0]*np.pi
        ph = np.random.rand(1)[0]*2*np.pi
        print "th = %s   ph = %s" % (th, ph)
        angles += [ph, th]
    assert len(angles) == 2*nrT
    return np.array(angles)

def initial_thomson_points_cart(nrT):
    angles = initial_thomson_points_sph(nrT)
    pts = np.zeros(3*nrT)
    nang = len(angles)/2
    for i in range(0, nang):
        th,ph = angles[2*i:2*(i+1)]
        x = R*np.sin(th)*np.cos(ph)
        y = R*np.sin(th)*np.sin(ph)
        z = R*np.cos(th)
        pts[3*i:3*(i+1)] = np.array([x,y,z])
    return pts

def cart2atomlist(x):
    n = x.shape[0]
    npts = n/3
    atomlist = []
    for i in range(0, npts):
        posi = x[3*i:3*(i+1)]
        atomlist.append( (6, posi) )
    return atomlist

########################## SPHERICAL COORDINATES WITHOUT CONSTRAINTS
def spherical2cartesian(R, angles):
    th, ph = angles
    sth = np.sin(th)
    cth = np.cos(th)
    sph = np.sin(ph)
    cph = np.cos(ph)
    x = R*sth*cph
    y = R*sth*sph
    z = R*cth
    return np.array([x,y,z])

def f_sph_old(x):
    n = x.shape[0]
    npts = n/2
    coul = 0.0
    for i in range(0, npts):
        for j in range(0, npts):
            if i==j:
                continue
            # R is a global variable
            ri = spherical2cartesian(R, x[2*i:2*(i+1)])
            rj = spherical2cartesian(R, x[2*j:2*(j+1)])
            coul += 1.0/la.norm(ri-rj)
    #
    coul *= 0.5
    print "coul = %s" % coul
    return coul

def f_sph(x):
    n = x.shape[0]
    npts = n/2
    coul = 0.0
    for i in range(0, npts):
        for j in range(0, npts):
            if i==j:
                continue
            th_i, ph_i = x[2*i:2*(i+1)]
            th_j, ph_j = x[2*j:2*(j+1)]

            rij2 = 2*(1.0-( np.sin(th_i)*np.sin(th_j)*np.cos(ph_i-ph_j)\
                            +np.cos(th_i)*np.cos(th_j)))
            rij = R*np.sqrt(rij2)
            coul += 1.0/rij
    coul *= 0.5
    err = abs(coul - f_sph_old(x))/abs(coul)
    assert err < 1.0e-8, "err = %s" % err
    return coul

def sph2atomlist(x, R):
    n = x.shape[0]
    npts = n/2
    assert 2*npts == n
    atomlist = []
    for i in range(0, npts):
        angles = x[2*i:2*(i+1)]
        posi = spherical2cartesian(R, angles)
        atomlist.append( (6, posi) )
    return atomlist
    
def delta(i,j):
    if i == j:
        return 1.0
    else:
        return 0.0

def compute_rij2(x):
    n = x.shape[0]
    npts = n/2
    r2 = np.zeros((npts, npts))
    r2dth = np.zeros((npts,npts,npts))
    r2dph = np.zeros((npts,npts,npts))
    for i in range(0, npts):
        for j in range(0, npts):
            if i==j:
                continue
            th_i, ph_i = x[2*i:2*(i+1)]
            th_j, ph_j = x[2*j:2*(j+1)]

            r2[i,j] = 2*R**2*(1.0-( np.sin(th_i)*np.sin(th_j)*np.cos(ph_i-ph_j)\
                            +np.cos(th_i)*np.cos(th_j)))
            if abs(r2[i,j]) < 1.0e-10:
                raise ValueError("rij^2 = %s      i=%s j=%s" % (r2[i,j], i, j))
            # derivatives
            for a in range(0, npts):
                if not (a == i or a == j):
                    continue
                # derivative of r_ij^2 with respect to theta_a
                r2dth[i,j,a] = 2*R**2*(delta(a,i)*(-np.cos(ph_i-ph_j)*np.cos(th_i)*np.sin(th_j) \
                                                   +np.sin(th_i)*np.cos(th_j)) \
                                      +delta(a,j)*(-np.cos(ph_i-ph_j)*np.sin(th_i)*np.cos(th_j) \
                                                   +np.cos(th_i)*np.sin(th_j)))
                
                # derivative of r_ij^2 with respect to phi_a
                r2dph[i,j,a] = 2*R**2*np.sin(th_i)*np.sin(th_j)*np.sin(ph_i-ph_j) \
                    * (delta(a,i)-delta(a,j))
    return r2,r2dth,r2dph

def f_rij2(x):
    n = x.shape[0]
    npts = n/2
    r2,r2dth,r2dph = compute_rij2(x)
    # build Coulomb energy
    coul = 0.0
    grad = np.zeros(n)
    for i in range(0, npts):
        for j in range(0, npts):
            if i == j:
                continue
            coul += 0.5*pow(r2[i,j], -0.5)
            # gradient
            for a in range(0, npts):
                grad[2*a]   -= 1.0/4.0 * pow(r2[i,j],-1.5) * r2dth[i,j,a]
                grad[2*a+1] -= 1.0/4.0 * pow(r2[i,j],-1.5) * r2dph[i,j,a]
    print "coul = %s" % coul
#    print "grad = %s" % grad
    print "|grad| = %s" % la.norm(grad)
    return coul, grad

if __name__ == "__main__":
    """
    # CARTESIAN
    x0 = initial_thomson_points_cart(12)
    print "err = %s" % optimize.check_grad(f_cart, grad_cart, x0)
#    xmin, fmin, d = optimize.fmin_l_bfgs_b(f_cart, x0, grad)
    xmin = optimize.fmin(f_cart, x0,maxfun=100000000000000000000, maxiter=100000000000000000)
    atomlist = cart2atomlist(4*xmin)
    XYZ.write_xyz("/tmp/thomson_spherical.xyz", [atomlist])
    """
    # SPHERICAL
    R = 5.0
    x0 = initial_thomson_points_sph(16)
#    print "err = %s" % optimize.check_grad(f_cart, grad_cart, x0)
#    xmin = optimize.fmin(f_sph, x0,maxfun=100000000000000000000, maxiter=100000000000000000)

    xmin, fmin, d = optimize.fmin_l_bfgs_b(f_rij2, x0,maxfun=100000000000000000000, maxiter=100000000000000000, approx_grad=0, pgtol=1.0e-8)
    atomlist = sph2atomlist(xmin, R)
    XYZ.write_xyz("/tmp/thomson_spherical.xyz", [atomlist])
 
    print "fmin = %s" % (2*fmin)

#    print "err = %s" % optimize.check_grad(lambda x: f_rij2(x)[0], lambda x: f_rij2(x)[1], x0)
