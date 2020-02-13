"""
Tranformations between 
   cartesian and internal (Z-matrix)
coordinates

copied from FCIQMC-0.0.6
"""
from DFTB import AtomicData
from DFTB import XYZ
from DFTB import utils

from numpy import array, dot, sum, zeros, arccos, cross, arctan2, argsort, pi, identity, cos, sin, copy, hstack, sqrt, arccos, loadtxt
from numpy.linalg import norm, eigh, inv, det
import numpy as np
import numpy.linalg as la

def center_of_mass(masses, pos):
    """
    find the center of mass
    R = (sum_i m_i * r_i)/sum(m_i)
    """
    xm = pos*masses
    com = np.array([np.sum(xm[0::3]), np.sum(xm[1::3]), np.sum(xm[2::3])])/np.sum(masses[0::3])
    return com

def shift_to_com(pos,masses):
    """
    shift center of mass to the origin
    """
    com = center_of_mass(masses, pos)
    pos_shifted = np.zeros(pos.shape)
    for i in xrange(0, len(pos)/3):
        pos_shifted[3*i:3*i+3] = pos[3*i:3*i+3] - com
    return pos_shifted
    
def inertial_tensor(masses, pos):
    """
    find tensor of inertia which relates angular velocity omega and angular momentum
    of a rigid body: L = I*w
    """
    I = zeros((3,3))
    Nat = len(pos)/3
    assert len(masses) == len(pos)
    for i in xrange(0, Nat):
        xi,yi,zi = pos[3*i:3*i+3]
        ri2 = xi*xi+yi*yi+zi*zi
        mi = masses[3*i]
        I[0,0] += mi*(ri2-xi*xi)
        I[0,1] -= mi*xi*yi
        I[0,2] -= mi*xi*zi
        I[1,1] += mi*(ri2-yi*yi)
        I[1,2] -= mi*(yi*zi)
        I[2,2] += mi*(ri2-zi*zi)
    I[1,0] = I[0,1]
    I[2,0] = I[0,2]
    I[2,1] = I[1,2]
#    print "Inertial Tensor"
#    print "==============="
#    print I
    return I

def angular_momentum(masses, pos, vel):
    """
    compute total angular momentum
    """
    L = np.zeros(3)
    Nat = len(pos)/3
    for i in range(0, Nat):
        L += masses[3*i] * np.cross(pos[3*i:3*(i+1)], vel[3*i:3*(i+1)])
    return L

def linear_momentum(masses, vel):
    P = np.zeros(3)
    Nat = len(vel)/3
    for i in range(0, Nat):
        P += masses[3*i:3*(i+1)] * vel[3*i:3*(i+1)]
    return P
        
def angular_velocity(L,I):
    """
    angular velocity vector that relates angular momentum L and the
    tensor of inertia I
    
          L = I.w
    """
    omega = np.dot(la.inv(I), L)
    return omega

def eliminate_rotation(vel, omega, pos):
    Nat = len(vel)/3
    for i in range(0, Nat):
        vel[3*i:3*(i+1)] -= np.cross(omega,pos[3*i:3*(i+1)])
    return vel
        
def eliminate_translation(masses, vel, P):
    Nat = len(vel)/3
    total_mass = np.sum(masses[::3])
    for i in range(0, Nat):
        vel[3*i:3*(i+1)] -= P/total_mass
    return vel
        

def rotation_matrix(axis,theta):
    """
    compute the rotation matrix that rotates a 3D vector around
    an axis k=(kx,ky,kz) by an angle theta

    Parameters:
    ===========
    axis: numpy array [kx,ky,kz]
    theta: rotation angle in radians

    Returns:
    ========
    rotation matrix R

    stolen from http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    """
    axis = axis/sqrt(dot(axis,axis))
    a = cos(theta/2)
    b,c,d = -axis*sin(theta/2)
    return array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def reflection_matrix(axis):
    """
    compute the reflection matrix that reflects a 3D vector through
    the plane orthogonal to axis.

    Parameters:
    ===========
    axis: numpy array [kx,ky,kz]

    Returns:
    ========
    reflection matrix R
    """
    a = axis
    # Rij = delta_ij - 2 * (ai*aj)/|a|^2
    # see http://en.wikipedia.org/wiki/Reflection_(mathematics)
    R = np.identity(3) - 2*np.outer(a,a)/np.dot(a,a)
    return R

#The following functions are used to define internal coordinates

def distance(atom1, atom2):
    """
    distance between two atomic centers

    Parameters:
    ===========
    atom1:  tuple (Z1,(x1,y1,z1))
    atom2:  tuple (Z2,(x2,y2,z2))

    Returns:
    ========
    d = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
    """
    (Z1,pos1) = atom1
    (Z2,pos2) = atom2
    d = norm(array(pos1) - array(pos2))
    return d

def angle(atom1, atom2, atom3):
    """
    angle between vectors atom1->atom2 and atom1->atom3
    """
    v12 = array(atom2[1]) - array(atom1[1])
    v13 = array(atom3[1]) - array(atom1[1])
    a = arccos(dot(v12,v13)/(norm(v12)*norm(v13)))
    return a

def dihedral_angle(atom1, atom2, atom3, atom4):
    """
    dihedral angle between the planes spanned by atom1->atom2, atom2->atom3
    and atom2->atom3, atom3->atom4
    """
    b1 = array(atom2[1]) - array(atom1[1])
    b2 = array(atom3[1]) - array(atom2[1])
    b3 = array(atom4[1]) - array(atom3[1])
    n12 = cross(b1,b2)
    n23 = cross(b2,b3)
    dihedral = arctan2(dot(cross(n12,n23),b2/norm(b2)), dot(n12,n23))
    return dihedral

# copied from DFTB-0.0.2/Optimization/InternalCoordinates.py

def tripod2EulerAngles(xaxis,yaxis,zaxis):
    """
    Find the Euler angles a,b,g that would rotate the right-handed orthonormal
    set of axes x,y,z into the unrotated set of axes [1,0,0],[0,1,0],[0,0,1].

    see section "Geometric derivation" in https://en.wikipedia.org/wiki/Euler_angles

    Parameters:
    ===========
    xaxis,yaxis,zaxis: 3 vectors specifying the axes
    
    Returns:
    ========
    a,b,g: Euler angles such that R(a,b,g)*[1,0,0] = x and similarly for the y- and z-axis.
    """
#    assert abs(dot(xaxis,yaxis)) < 1.0e-10
#    assert abs(dot(yaxis,zaxis)) < 1.0e-10
#    assert abs(dot(xaxis,zaxis)) < 1.0e-10
#    assert abs(norm(xaxis) - 1.0) < 1.0e-10
#    assert abs(norm(yaxis) - 1.0) < 1.0e-10
#    assert abs(norm(zaxis) - 1.0) < 1.0e-10

    """
    print "Zaxis = %s" % zaxis[2]
    b = arccos(zaxis[2])
    # nodal line
    N = cross(array([0,0,1]), zaxis)
    assert(abs(norm(N) - 1.0) < 1.0e-10)
    N /= norm(N)
    a = arctan2(N[1],N[0])-pi/2.0 # projection of nodal line on space fixed x- and y-axis
    g = arctan2(dot(xaxis,N), dot(yaxis,N))
    
    # make sure angles stay in range
    a = a % (2*pi)
#    b = b % pi
    g = g % (2*pi)
    
    return (a,b,g)
    """
#    print "X-axis = %s" % xaxis
#    print "Y-axis = %s" % yaxis
#    print "Z-axis = %s" % zaxis

    b = arccos(min(zaxis[2],1.0))
    if (b == 0.0):
        a = arctan2(xaxis[1],xaxis[0])
        g = 0.0
    elif (b == pi):
        a = -arctan2(xaxis[1],xaxis[0])
        g = 0.0        
    else:
        a = arctan2(zaxis[0], -zaxis[1])
        g = arctan2(xaxis[2], yaxis[2])
#    print "Euler angles (a,b,g) = %.4f, deg  %.4f deg  %.4f deg" % (a*180.0/pi, b*180.0/pi, g*180.0/pi)
    
    return (a,b,g)

def rotation2EulerAngles(R, convention="z-y-z"):
    """
    extract the Euler angles from a rotation matrix (in z-y-z convention)
    """
    if convention == "z-y-z":
        b = np.arccos(R[2,2])
        if b == 0:
            a = np.arctan2(R[0,1], R[0,0])
            g = 0.0
        elif (b == np.pi):
            a = np.arctan2(R[0,1], -R[0,0])
            g = 0.0
        else:
            a = np.arctan2(R[1,2], -R[0,2])
            g = np.arctan2(R[2,1], R[2,0])
    else:
        raise ValueError("Convention '%s' not implemented!" % convention)
    return (a,b,g)

def EulerAngles2Rotation(a,b,g, convention="my own"):
    """
    construct the rotation matrix around the axes x,y,z that are fixed in space (extrinsic rotation):
#        Rz(-a)*Ry(-b)*Rz(-g)
        Rz(g)*Rx(b)*Rz(a)

    Parameters:
    ===========
    a,b,g: Euler angles
    convention: 'z-x-z', 'z-y-z' or 'my own'

    Returns:
    ========
    3x3 rotation matrix
    """
    def Rx(phi):
        """clock-wise rotation by angle phi around the x-axis"""
        R = identity(3)
        R[1,1] = cos(phi)
        R[1,2] = sin(phi)
        R[2,1] = -sin(phi)
        R[2,2] = cos(phi)
        return R
    def Ry(phi):
        """clock-wise rotation by angle phi around the y-axis"""
        R = identity(3)
        R[0,0] = cos(phi)
        R[0,2] = -sin(phi)
        R[2,0] = sin(phi)
        R[2,2] = cos(phi)
        return R
    def Rz(phi):
        """clock-wise rotation by angle phi around the z-axis"""
        R = identity(3)
        R[0,0] = cos(phi)
        R[0,1] = sin(phi)
        R[1,0] = -sin(phi)
        R[1,1] = cos(phi)
        return R
    #return dot(Rz(a),dot(Rx(b),Rz(g)))
    if convention == "z-x-z":
        return dot(Rz(a),dot(Rx(b),Rz(g)))
        #return dot(Rz(a),dot(Rx(b),Rz(g)))
    elif convention == "z-y-z":
        #return dot(Rz(g),dot(Ry(b),Rz(a)))
        return dot(Rz(a),dot(Ry(b),Rz(g)))
    elif convention == "my own":
        # Initially I did not use one of the usual conventions.
        # Now, unfortunately the code for determining molecular symmetry
        # only works when the rotations are concatenated in the
        # wrong order.
        return dot(Rz(g),dot(Rx(b),Rz(a)))
    else:
        raise ValueError("Convention '%s' not understood!" % convention)
    
def construct_tripod(rA,rB,rC):
    """
    construct a right-handed orthonormal coordinate system from three
    position vectors. 

    The first axis lies along the vector joining rA and rB. 
    The second axis is perpendicular the the first axis
    and spans the plane defined by the three points. The third axis forms
    a right hand tripod with the first two axes.

    Parameters:
    ===========
    rA, rB, rC:  tuple of 3 3D numpy arrays

    Returns:
    ========
    n, a, b: tuple of 3 3D numpy arrays, normalized axes
    """
    # construct tripod n,a,b
    N = rB-rA
    n = N/norm(N)
    A = rC-rB - dot(rC-rB,n)*n
    a = A/norm(A)
    b = cross(n,a)
    return n,a,b

def cartesian2Zmatrix(atomlist):
    """
    convert a list of cartesian positions of atoms into
    relative distances, bond angles and dihedral angles

    Parameters:
    ===========
    atomlist: list of tuples (Zi,(xi,yi,zi)) with cartesian positions of atoms

    Returns:
    ========
    zmat: for 1 atom    []
          for 2 atoms   [r12]
          for 3 atoms   [r12, 
                         r23, ang123]
          for more than 4 atoms
                        [r12, 
                         r23, ang123,
                         ...
                         rDC, angBCD, dihedralABCD,
                         ...]
    """
    zmat = []
    Nat = len(atomlist)
    print "Internal Coordinates"
    print "===================="
    if Nat >= 1:
        print "%s" % atomlist[0][0]
    if Nat >= 2:
        r = distance(atomlist[0], atomlist[1])
        zmat += [r]
        print "%s   %.7f" % (atomlist[1][0], r)
    if Nat == 3:
        # for 3 atoms it is more intuitive to define the bond angle
        # as ang(rBA,rBC) instead of ang(rAB,rAC)
        r = distance(atomlist[0], atomlist[2]) # ??
        # r = distance(atomlist[1], atomlist[2])
        ang = angle(atomlist[0], atomlist[1], atomlist[2])
        zmat += [r, ang]
        print "%s   %.7f   %.3f deg" % (atomlist[2][0], r, ang*180.0/pi)            
    elif Nat > 3:
        r = distance(atomlist[1], atomlist[2])
        ang = angle(atomlist[1], atomlist[0], atomlist[2])
        zmat += [r, ang]
        print "%s   %.7f   %.3f deg" % (atomlist[2][0], r, ang*180.0/pi)
    if Nat >= 4:
        for i in range(3,len(atomlist)):
            r = distance(atomlist[i-1], atomlist[i])
            ang = angle(atomlist[i-1], atomlist[i-2], atomlist[i])
            dihedral = dihedral_angle(atomlist[i-3], atomlist[i-2], atomlist[i-1], atomlist[i])
            zmat += [r,ang,dihedral]
            print "%s   %.7f   %.3f deg   %.3f deg" % (atomlist[i][0], r, ang*180.0/pi, dihedral*180.0/pi)
    return zmat

def Zmatrix_labels(atomlist):
    """
    A Z-matrix is a list of internal coordinates. This function assigns names to each internal coordinate, depending
    on whether they belong to bond lengths, angles or dihedral angles. Changing the order of the atoms in the cartesian
    xyz-file will change the Z-matrix, too.

    Parameters:
    ===========
    atomlist: list of tuples (Zi,(xi,yi,zi)) with cartesian positions of atoms

    Returns:
    ========
    labels: labels[i] gives a meaning to the i-th item in the Z-matrix, It can be
      "r(A-B)" or "angle(A-B-C)" or "dihedral(A-B-C-D)"
    """
    def r_name(i, j):
        ati = AtomicData.atom_names[atomlist[i][0]-1]
        atj = AtomicData.atom_names[atomlist[j][0]-1]
        name = "r(%s%d-%s%d)" % (ati, i, atj, j)
        return name
    def angle_name(i, j, k):
        ati = AtomicData.atom_names[atomlist[i][0]-1]
        atj = AtomicData.atom_names[atomlist[j][0]-1]
        atk = AtomicData.atom_names[atomlist[k][0]-1]
        name = "angle(%s%d-%s%d-%s%d)" % (ati, i, atj, j, atk, k)
        return name
    def dihedral_name(i, j, k, l):
        ati = AtomicData.atom_names[atomlist[i][0]-1]
        atj = AtomicData.atom_names[atomlist[j][0]-1]
        atk = AtomicData.atom_names[atomlist[k][0]-1]
        atl = AtomicData.atom_names[atomlist[l][0]-1]
        name = "dihedral(%s%d-%s%d-%s%d-%s%d)" % (ati, i, atj, j, atk, k, atl, l)
        return name
        
    labels = []
    Nat = len(atomlist)
    if Nat >= 1:
        pass
    if Nat >= 2:
        labels += [r_name(0,1)]
    if Nat == 3:
        # for 3 atoms it is more intuitive to define the bond angle
        # as ang(rBA,rBC) instead of ang(rAB,rAC)
        labels += [r_name(0,2), angle_name(0,1,2)]
    elif Nat > 3:
        labels += [r_name(1,2), angle_name(1,0,2)]
    if Nat >= 4:
        for i in range(3,len(atomlist)):
            labels += [r_name(i-1,i), angle_name(i-1,i-2,i), dihedral_name(i-3,i-2,i-1,i)]
    return labels

def Zmatrix2cartesian(zmat):
    """
    Convert a Z-matrix to cartesian coordinates. The first atom is placed at the origin,
    the second on the x-axis and the third in the xy-plane. The positions of the other
    atoms are fixed by the internal coordinates.

    Parameters:
    ===========
    zmat: Z-matrix as returned by cartesian2Zmatrix()

    Returns:
    ========
    pos: 3*Natoms-dimensional array with cartesian positions of the atoms
    """

    positions = []
    if len(zmat) >= 0:
        posA = array([0,0,0])
        positions.append(posA)
    if len(zmat) >= 1:
        r12 = zmat[0] 
        posB = array([r12,0,0])
        positions.append(posB)
    if len(zmat) >= 2:
        r23 = zmat[1]
        ang123 = zmat[2]
        if len(zmat) == 3:
            # for triatomics A-B-C the angle is ang(rBA,rBC) and not ang(AB,BC)
            posC = array([r23*cos(ang123), r23*sin(ang123), 0])
        else:
            posC = array([r12 - r23*cos(ang123), r23*sin(ang123), 0])
        positions.append(posC)
    if len(zmat) >= 3:
        for i in range(1, len(zmat)/3):
            r = zmat[3*i]
            ang = zmat[3*i+1]
            dihedral = zmat[3*i+2]
            # b1 and b2 span the first plane
            b1 = posB - posA
            b1 /= norm(b1)
            b2 = posB - posC
            b2 /= norm(b2)
            # create the second plane by rotating b1 around b2 by the dihedral angle
            b3 = dot(rotation_matrix(b2, dihedral),b1)
            # n is orthogonal to the second plane
            n = cross(b2,b3)
            # finally, b3 is rotated such that it has the correct angle with b2
            b3 = dot(rotation_matrix(n, ang),b2)
            # shift atom from position C along b3
            posD = posC + r*b3
            positions.append(posD)

            # check that bond lengths, angles and dihedral angles are correct
            r_check = norm(posD - posC)
            assert(abs(r_check - r) < 1.0e-10)
            ang_check = angle((0,posC), (0,posB), (0,posD))
            assert(abs(ang_check - ang) < 1.0e-10)
            dih_check = dihedral_angle((0,posA), (0,posB), (0,posC), (0,posD))
#            assert(abs(dih_check - dihedral) < 1.0e-10)

            posA = posB
            posB = posC
            posC = posD
            
    pos = hstack(positions)
    return pos
            
            

def euler_angles_inertia(masses, pos):
    """
    Find the Euler angles describing the orientation of
    the atoms relative to the principal axes of inertia.
    The angles are obtained from the rotation matrix that
    diagonalizes the tensor of inertia.

    Parameters:
    ===========
    masses: list with masses of each atom in atomic units
    pos: pos[3*i:3*i+3] are the cartesian positions of atom i

    Returns:
    ========
    Euler angles (a,b,g)   such that 
       R(a,b,g)*(1st largest principle moment of inertia) = [1,0,0]
       R(a,b,g)*(2nd largest principle moment of inertia) = [0,1,0]
       R(a,b,g)*(3rd largest principle moment of inertia) = [0,0,1]
    """
    I = inertial_tensor(masses, pos)
    if abs(det(I)) < 1.0e-10:
        # detect linear molecule
        print "LINEAR MOLECULE"
        zaxis = pos[3:6] - pos[0:3]
        zaxis /= norm(zaxis)
        # projection on xaxis
        projz = 0.0
        Nat = len(masses)/3
        for i in range(0, Nat):
            projz += (i+1) * dot(zaxis, pos[i*3:i*3+3])
#        print "projz = %s" % projz
        if projz < 0.0:
            zaxis *= -1
#            print "FLIP Z"
#        print "zaxis = %s" % zaxis
        beta = arccos(zaxis[2])
        if (beta % pi == 0.0):
            alpha = 0.0
        else:
            alpha = arctan2(zaxis[0], -zaxis[1])
        gamma = 0.0
    else:
        principal_moments, principal_axes = eigh(I)
        sort_indx = argsort(principal_moments)
        moments = principal_moments[sort_indx[::-1]]
        axes = principal_axes[:,sort_indx[::-1]]
#        print "Principal Axes of Inertia"
#        print "========================="
#        for i in [0,1,2]:
#            print "I[%s] = %e    axis[%s] = %s" % (i, moments[i], i, axes[:,i])
        # phases of eigen vectors are not unique, fix them by some convention
        # find projection of atom positions on axes
        proj = array([0.0,0.0,0.0])
        Nat = len(masses)/3
        for xyz in [0,1,2]:
            for i in range(0, Nat):
                proj[xyz] += (i+1) * dot(axes[:,xyz], pos[i*3:i*3+3])
#        print "proj = %s" % proj
        if proj[2] < 0.0:
            axes[:,2] *= -1
#            print "FLIP Z"
        if proj[1] < 0.0:
            axes[:,1] *= -1
#            print "FLIP Y"
        if dot(axes[:,2], cross(axes[:,0], axes[:,1])) < 0.0:
            # create right handed coordinate system
            axes[:,0] *= -1
#            print "FLIP X"

        assert(dot(axes[:,2], cross(axes[:,0], axes[:,1])) > 0.0)
        alpha,beta,gamma = tripod2EulerAngles(axes[:,0], axes[:,1],axes[:,2])
        #
        R = EulerAngles2Rotation(alpha, beta, gamma)
#        print "R"
#        print R
#        print "inv(axes)"
#        print inv(axes)
#        assert(sum(abs(R-inv(axes))) < 1.0e-6)
#        print "R * I"
#        print dot(R, axes[:,0])
#        print dot(R, axes[:,1])
#        print dot(R, axes[:,2])
        #
#    print "Euler angles     %.4f  %.4f  %.4f" % (alpha, beta, gamma)
    return (alpha,beta,gamma)

def molpro_standard_orientation(masses, pos):
    """
    reorient the molecule so that center of mass is at the origin 
    and that the principal axes of inertia correspond to the axes.

    Parameters:
    ===========
    masses: list with masses of each atom in atomic units
    pos: pos[3*i:3*i+3] are the cartesian positions of atom i

    Returns:
    ========
    pos: reoriented cartesian positions
    """
#    print "STANDARD ORIENTATION"
    Nat = len(pos)/3
    pos_std = zeros(pos.shape)
    # shift center of mass to origin
    cm = center_of_mass(masses, pos)
    for i in range(0, Nat):
        pos_std[3*i:3*i+3] = pos[3*i:3*i+3] - cm
    (a,b,g) = euler_angles_inertia(masses, pos_std)
    R = EulerAngles2Rotation(a,b,g)
    for i in xrange(0, Nat):
        pos_std[3*i:3*i+3] = dot(R, pos_std[3*i:3*i+3])
    #
    print "\nMolecule in standard orientation"
    print "================================\n"
    print "The origin is shifted to the center of mass"
    print "and the axes are rotated to coincide with the"
    print "principal axes of inertia.\n"
    print "                in bohr:"
    print     "Atom          X           Y           Z"
    for i in range(0, Nat):
        print "  %s  " % str(i).ljust(3),
        print "%+4.7f  %+4.7f  %+4.7f" % tuple(pos_std[3*i:3*(i+1)])
    print "                in Angstrom:"
    print     "Atom          X           Y           Z"
    for i in range(0, Nat):
        print "  %s  " % str(i).ljust(3),
        print "%+4.7f  %+4.7f  %+4.7f" % tuple(pos_std[3*i:3*(i+1)]*AtomicData.bohr_to_angs)
    print ""
    #
    I = inertial_tensor(masses, pos_std)
#    print "Inertial tensor in standard orientation"
#    print I
    (a0,b0,g0) = euler_angles_inertia(masses, pos_std)
#    print "Euler angles in standard orientation"
#    print "%.8f deg   %.8f deg   %.8f deg" % (a0 * 180.0/pi, b0 * 180.0/pi, g0 * 180.0/pi)
#    assert(abs(a0)+abs(b0)+abs(g0) < 1.0e-6)
    cm0 = center_of_mass(masses, pos_std)
#    print "Center of mass in standard orientation"
#    print cm0
    return pos_std

def molecular_frame_transformation(atomlist):
    """
    The molecule is shifted to the center of mass and its principle axes of inertia are aligned 
    with the coordinate axes. This standard orientation defines the molecular frame.
    The translation vector and Euler angles needed to transform the geometry from the
    molecular frame to the original frame are also returned. 

    Returns:
    ========
    atomlist_std: molecular geometry in standard orientation
    (a,b,g): Euler angles in z-y-z convention
    cm: 3D vector with center of mass
    """
    pos = XYZ.atomlist2vector(atomlist)
    masses = AtomicData.atomlist2masses(atomlist)
    # shift center of mass to origin
    Nat = len(atomlist)
    pos_std = np.zeros(pos.shape)
    cm = center_of_mass(masses, pos)
    for i in range(0, Nat):
        pos_std[3*i:3*i+3] = pos[3*i:3*i+3] - cm
    (a,b,g) = euler_angles_inertia(masses, pos_std)
    R = EulerAngles2Rotation(a,b,g)
    for i in xrange(0, Nat):
        pos_std[3*i:3*i+3] = np.dot(R, pos_std[3*i:3*i+3])
    atomlist_std = XYZ.vector2atomlist(pos_std, atomlist)
    # The MolecularCoord module uses a strange convention for Euler angles.
    # Therefore we extract the Euler angles in z-y-z convention directly
    # from the rotation matrix, that rotates the molecule from the standard
    # orientation into the original orientation
    Rinv = R.transpose()
    a,b,g = rotation2EulerAngles(Rinv, convention="z-y-z")
    return atomlist_std, (a,b,g), cm

def transform_molecule(atomlist, euler_angles, translation_vector):
    """
    rotates and shifts the molecular geometry.
    Each atomic position x is transformed according to:
        x ->  R*x + T

    Parameters:
    ===========
    atomlist: list of (Zi,[xi,yi,zi]) 
    euler_angles: tuple (a,b,g) of Euler angles in z-y-z convention
    translation_vector: 3D vector T

    Returns:
    ========
    atomlist_transformed: transformed molecular geometry
    """
    #
    a,b,g = euler_angles
    R = EulerAngles2Rotation(a,b,g, convention="z-y-z")
    atomlist_transformed = []
    # shift and rotate nuclear geometry
    for Zi,posi in atomlist:
        posi_transf = np.dot(R, np.array(posi)) + translation_vector
        atomlist_transformed.append( (Zi,posi_transf) )
    
    return atomlist_transformed

def internal2cartesian(masses, internal):
    Nat = len(internal)/3
    zmat = internal[:-6]
    cm = internal[-6:-3]
    a,b,g = internal[-3:]
    pos = Zmatrix2cartesian(zmat)
    pos_std = molpro_standard_orientation(masses, pos)
    # rotate geometry from the standard orientation back to the original orientation
    R = inv(EulerAngles2Rotation(a,b,g))    
    for i in range(0, Nat):
        pos[3*i:3*i+3] = dot(R, pos_std[3*i:3*i+3]) + cm
    return pos


def cartesian2internal(masses, pos, ref_atomlist, debug=0):
    """
    convert cartesian coordinates to internal coordinates + translation + rotation.
    Because global rotation and translation are included the number
    of internal coordinates equals the number of cartesian coordinates.

    Parameters:
    ===========
    masses: list with masses of each atom in atomic units
    pos: pos[3*i:3*i+3] are the cartesian positions of atom i
    ref_atomlist: list of tuples (Zi, (xi,yi,zi), positions do not have to be 
       the same as in pos, only the atomic numbers are needed

    Optional:
    =========
    debug: if > 0, 

    Returns:
    ========
    numpy array with zmat, center of mass and Euler angles
    """
    atomlist = XYZ.vector2atomlist(pos, ref_atomlist)
    Nat = len(atomlist)

    print "Original positions"
    print "=================="
    for i in range(0, Nat):
        print "%s   " % i,
        print "%.4f %.4f %.4f" % tuple(pos[3*i:3*i+3])
#    print "Euler angles for rotating original geometry into standard orientation"
    cm2 = center_of_mass(masses, pos)
#    print "Center of mass of original positions"
#    print cm2
    pos1 = zeros(pos.shape)
    for i in range(0, Nat):
        pos1[3*i:3*i+3] = pos[3*i:3*i+3] - cm2
    (a2,b2,g2) = euler_angles_inertia(masses, pos1)
    # test transformation
    pos_std1 = molpro_standard_orientation(masses, pos)
    pos = copy(pos)
    zmat = cartesian2Zmatrix(atomlist)
    pos2 = Zmatrix2cartesian(zmat)
    pos_std2 = molpro_standard_orientation(masses, pos2)
#    XYZ.write_xyz("std1.xyz", [XYZ.vector2atomlist(pos_std1, ref_atomlist)])
#    XYZ.write_xyz("std2.xyz", [XYZ.vector2atomlist(pos_std2, ref_atomlist)])
    pos2 = pos_std2
#    R2 = inv(EulerAngles2Rotation(a2,b2,g2))
    R2 = inv(EulerAngles2Rotation(a2,b2,g2))
#    print "R2"
#    print R2
#    print "pos2"
#    print pos2
    posrecovered = zeros(pos.shape)
#    print "cm2 = %s" % cm2
    for i in range(0, Nat):
        posrecovered[3*i:3*i+3] = dot(R2, pos2[3*i:3*i+3]) + cm2
#    print "Recovered positions"
#    print "==================="
#    for i in range(0, Nat):
#        print "%s   " % i,
#        print "%.4f %.4f %.4f" % tuple(posrecovered[3*i:3*i+3])
    #
#    XYZ.write_xyz("recovered.xyz", [XYZ.vector2atomlist(posrecovered, ref_atomlist)])
    cm3 = center_of_mass(masses, posrecovered)
#    print "center of mass of recovered positions"
#    print cm3
    (a3,b3,g3) = euler_angles_inertia(masses, posrecovered)
#    print "Euler angles of recovered positions"
#    print "%.4f deg   %.4f deg   %.4f deg" % (a3 * 180.0/pi, b3 * 180.0/pi, g3 * 180.0/pi)

    err = sum(abs(posrecovered - pos))
#    print "error = %s" % err
    assert(err < 1.0e-5)
    #
    cm = list(center_of_mass(masses, pos))
    # shift center of mass to origin
    for i in range(0, len(ref_atomlist)):
        pos[3*i:3*i+3] -= cm
    (a,b,g) = euler_angles_inertia(masses, pos)
    internal = zmat + cm + [a,b,g]
    return array(internal)

def jacobian_cartesian2internal(masses, pos, ref_atomlist, h=1.0e-8):
    """
    numerically find Jacobian matrix for transformation from cartesian to internal coordinates

    J_ij = d( chi_j )/d( x_i )

    where chi_1,chi_2,...,chi_3N are the internal coordinates
    are x_1,x_2,...,x_3N are the cartesian coordinates of N atoms
    
    Parameters:
    ===========
    masses: list of N masses in atomic units for each atom
    pos: pos[3*i:3*i+3] are the cartesian coordinates of atom i
    h: step for difference quotient
    
    Returns:
    ========
    3Nx3N numpy array with J_ij
    """
    Nat = len(pos)/3
    chi0 = cartesian2internal(masses, pos, ref_atomlist)
    Ncartesian = 3*Nat
    Ninternal = len(chi0)
    Jacobian = zeros((Ncartesian,Ninternal))
    print "Internal Coordinates"
    print chi0
    for i in range(0, Ncartesian):
        # forward difference quotient
        xplus = copy(pos)
        xplus[i] += h
        # 
        chiplus  = cartesian2internal(masses, xplus, ref_atomlist)
        print "Displaced Internal Coordinates"
        print chiplus
        print "i = %s" % i
        print "********"
        print "xplus"
        print xplus
        print "chiplus"
        print chiplus
        print "Center of mass"
        print chiplus[-6:-3]
        print "Euler angles"
        print chiplus[-3:]
        print "chi[x+h*%d] - chi[x]" % i
        print (chiplus - chi0) 
        Jacobian[i,:] = (chiplus - chi0) / h
    # inverse of jacobian
    invJacobian = zeros((Ninternal,Ncartesian))
    chi0 = cartesian2internal(masses, pos, ref_atomlist)
    x0 = pos
    for i in range(0, Ninternal):
        # forward difference quotient
        chiplus = copy(chi0)
        chiplus[i] += h
        # 
        xplus  = internal2cartesian(masses, chiplus)
        invJacobian[i,:] = (xplus - x0) / h
    return Jacobian, invJacobian

class MolecularCoords:
    def __init__(self, atomlist):
        self.atomlist = atomlist
        self.masses = zeros(3*len(self.atomlist))
        for i,(Zi,posi) in enumerate(self.atomlist):
            self.masses[3*i:3*(i+1)] = AtomicData.atom_masses[AtomicData.atom_names[Zi-1]]
    def setGeometry(self, atomlist):
        self.atomlist = atomlist
    def getGeometry(self):
        return self.atomlist
    def getCenterOfMass(self):
        pos = XYZ.atomlist2vector(self.atomlist)
        com = center_of_mass(self.masses, pos)
        print "Center of Mass"
        print "=============="
        print com
        return com
    def getInertialTensor(self):
        pos = XYZ.atomlist2vector(self.atomlist)
        return inertial_tensor(self.masses, pos)
    def getZmatrix(self):
        return cartesian2Zmatrix(self.atomlist)
    def getEulerAngles(self):
        pos = XYZ.atomlist2vector(self.atomlist)
        return euler_angles_inertia(self.masses, pos)
    def getInternalCoordinates(self):
        pos = XYZ.atomlist2vector(self.atomlist)
        return cartesian2internal(self.masses, pos, self.atomlist)
    def setInternalCoordinates(self, internal):
        pos = internal2cartesian(self.masses, internal)
        self.atomlist = XYZ.vector2atomlist(pos, self.atomlist)
    def getJacobian(self):
        pos = XYZ.atomlist2vector(self.atomlist)
#        pos = self.standard_orientation()
        jac, jac_inv = jacobian_cartesian2internal(self.masses, pos, self.atomlist)
        print "Jacobian  J_ij = d( q_j )/d( x_i )"
        print "====================================="
        print utils.annotated_matrix(jac,["x_%s" % i for i in range(0, len(jac[:,0]))], ["q_%s" % j for j in range(0, len(jac[0,:]))])

        print "Jacobian  J^(-1)_ij = d( x_j )/d( q_i )"
        print "====================================="
#        jac_inv = inv(jac)
        print utils.annotated_matrix(jac_inv,["q_%s" % i for i in range(0, len(jac_inv[:,0]))], ["x_%s" % j for j in range(0, len(jac_inv[0,:]))])


#        print "Jacobian  J^(-1)_ij = d( x_j )/d( q_i )"
#        print "====================================="
#        jac_inv = inv(jac)
#        print utils.annotated_matrix(jac_inv,["q_%s" % i for i in range(0, len(pos))], ["x_%s" % j for j in range(0, len(pos))])

        return jac, jac_inv
    def standard_orientation(self):
        pos = XYZ.atomlist2vector(self.atomlist)
        return molpro_standard_orientation(self.masses, pos)            

def standard_orientation(atomlist):
    """
    bring a molecular geometry into standard orientation
    """
    mol = MolecularCoords(atomlist)
    pos_std = mol.standard_orientation()
    atomlist_std = XYZ.vector2atomlist(pos_std, atomlist)
    return atomlist_std

def nact_cartesian2internal(atomlist, nact):
    mol = MolecularCoords(atomlist)
    jac, jac_inv = mol.getJacobian()
    for (I,J),nactIJ in nact.iteritems():
        #
        RxGrad = zeros(3)
        Grad = zeros(3)
        for i,(Zi, posi) in enumerate(atomlist):
            RxGrad += cross(array(posi), array(nactIJ[i][1]))
            Grad += array(nactIJ[i][1])
        print "Grad"
        print Grad
        print "RxGrad"
        print RxGrad
        #
        nactIJvec = XYZ.atomlist2vector(nactIJ)
        nactIJinternal = dot(jac_inv, nactIJvec)
        print "Non-adiabatic coupling vectors between state %s and %s" % (I,J)
        print "in cartesian coordinates"
        print nactIJvec
        print "in internal coordinates"
        print nactIJinternal

def grad_cartesian2internal(atomlist, gradients):
    mol = MolecularCoords(atomlist)
    jac, jac_inv = mol.getJacobian()
    for I,gradI in enumerate(gradients):
        gradIvec = XYZ.atomlist2vector(gradI)
        gradIinternal = dot(jac_inv, gradIvec)
        print "Gradient on state %s" % I
        print "in cartesian coordinates"
        print gradIvec
        print "in internal coordinates"
        print gradIinternal


