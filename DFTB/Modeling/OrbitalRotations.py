"""
When a molecule is rotated the orbitals have to be rotated too. This module
implements the rotation matrices for real s-,p- and d-orbitals.
"""

import numpy as np
from numpy import sin, cos
from itertools import groupby

from DFTB.Modeling import MolecularCoords as MolCo

def orbital_rotation_matrices(a,b,g):
    """
    Parameters:
    ===========
    a,b,g: Euler angles in z-y-z convention. 
           Warning: The Euler angles DFTB.MolecularCoords follow a different convention.

    Returns:
    ========
    rotation matrices for s-, p- and d-block

    The matrices were computed using the script 'DFTB/Mathematica/rotated_Yreal.nb'
    """
    rotS = np.array([1.0])

    sa, ca = np.sin(a), np.cos(a)
    sb, cb = np.sin(b), np.cos(b)
    sg, cg = np.sin(g), np.cos(g)
    
    rotP = np.array([
        [ca*cg-cb*sa*sg,    sb*sg,   cg*sa+ca*cb*sg  ],
        [sa*sb,             cb,          -ca*sb      ],
        [-cb*cg*sa-ca*sg,   cg*sb,   ca*cb*cg - sa*sg]])

    s2a,c2a = np.sin(2*a), np.cos(2*a)
    s2b,c2b = np.sin(2*b), np.cos(2*b)
    s2g,c2g = np.sin(2*g), np.cos(2*g)

    sb2,cb2 = np.sin(b/2.0), np.cos(b/2.0)

    sapg, capg = np.sin(a+g), np.cos(a+g)
    samg, camg = np.sin(a-g), np.cos(a-g)

    s2apg, c2apg = np.sin(2*a+g), np.cos(2*a+g)
    s2amg, c2amg = np.sin(2*a-g), np.cos(2*a-g)

    SQ3 = np.sqrt(3.0)
    
    rotD = np.array([
        [0.25*(4*c2a*cb*c2g-(3+c2b)*s2a*s2g),       sb*(ca*c2g-2*cb*cg*sa*sg),                     SQ3/2.0 * sb**2 * s2g,        0.25*(4*c2g*sa*sb + 2*ca*s2b*s2g),                    0.25*(4*cb*c2g*s2a+c2a*(3+c2b)*s2g)     ],
        [-2*cb2**3*c2apg*sb2-2*cb2*c2amg*sb2**3,    cb2**2*(-1+2*cb)*capg+(1+2*cb)*camg*sb2**2,    SQ3*cb*sb*sg,                 sb2**2*samg+2*cb*sb2**2*samg+cb2**2*(-1+2*cb)*sapg,   -sb2**2*sb*s2amg-2*cb2**3*sb2*s2apg     ],
        [-SQ3/2.0 * s2a*sb**2,                      SQ3/2.0*sa*s2b,                                0.25*(1+3*c2b),               -SQ3/2.0 * ca*s2b,                                    SQ3/2.0 * c2a*sb**2                     ],
        [-sb2**2*sb*s2amg+2*cb2**3*sb2*s2apg,       (1+2*cb)*sb2**2*samg-cb2**2*(-1+2*cb)*sapg,    SQ3*cb*cg*sb,                 cb2**2*(-1+2*cb)*capg-(1+2*cb)*camg*sb2**2,           -2*cb2**3*c2apg*sb2+2*cb2*c2amg*sb2**3  ],
        [-sb2**4*sin(2*(a-g))-cb2**4*sin(2*(a+g)),  sb2**2*sb*sin(a-2*g)-2*cb2**3*sb2*sin(a+2*g),  SQ3/2.0*c2g*sb**2,            2*cb2**3*cos(a+2*g)*sb2-2*cb2*cos(a-2*g)*sb2**3,      cb2**4*cos(2*(a+g))+cos(2*(a-g))*sb2**4]])
    # matrix for rotating coefficients of s-, p- and d-block
    orbital_rotations = {0: rotS.transpose(),
                         1: rotP.transpose(),
                         2: rotD.transpose()}
    return orbital_rotations
                        
def rotate_orbitals(atomlist, valorbs, orbs, euler_angles):
    """
    computes MO coefficients of a rotated orbital. The Euler angles
    follow the z-y-z convention.
    """
    a,b,g = euler_angles
    # find rotation matrices for s-,p- and d-orbitals
    orb_rot_byl = orbital_rotation_matrices(a,b,g)
    rotation_matrices = []
    # iterate over atoms
    for (Zi,posi) in atomlist:
        # iterate over all orbitals grouped by angular momentum shell
        for (n,l),lshell_orbs in groupby(valorbs[Zi], lambda nlm: (nlm[0], nlm[1])):
            assert len(list(lshell_orbs)) == 2*l+1
            rotation_matrices.append(orb_rot_byl[l])
    # rotate orbitals
    i = 0
    orbs_rotated = np.copy(orbs)
    for rot in rotation_matrices:
        dim = rot.shape[0]
        orbs_rotated[i:i+dim] = np.dot(rot, orbs[i:i+dim])
        i += dim

    return orbs_rotated


if __name__ == "__main__":
    a,b,g = 0.345345, 1.234, 0.56
    orb_rot = orbital_rotation_matrices(a,b,g)
    R = MolCo.EulerAngles2Rotation(a,b,g, convention="z-y-z")
    print "R"
    print R
#    print np.dot(R, R.transpose())
#    print "Rorb"
#    print orb_rot[1]
    mapping = {0: 1, 1: 2, 2: 0}
    Rmapped = np.zeros((3,3))
    for i in range(0, 3):
        for j in range(0, 3):
            Rmapped[mapping[i],mapping[j]] = orb_rot[1][i,j]
    print "Rmapped"
    print Rmapped
    print R-Rmapped
#    print np.dot(orb_rot[1], orb_rot[1].transpose())

    for l in range(0, 2+1):
        print "*** l=%d ***" % l
        Rorb = orb_rot[l]
        print "Rorb"
        print Rorb
        print np.dot(Rorb, Rorb.transpose())
        print np.round(np.dot(Rorb.transpose(), Rorb), 8)

        
