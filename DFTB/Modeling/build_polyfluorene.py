#!/usr/bin/env python
"""
assemble a polyfluorene oligomer of desired chain length
"""
import numpy as np
import numpy.linalg as la

from DFTB import XYZ, AtomicData
from DFTB.Modeling import MolecularCoords as MolCo

def shift_atoms(atomlist, vec):
    """
    shift all atoms by a displacement vector
    """
    atomlist_shifted = []
    for (Z,pos) in atomlist:
        pos = np.array(pos)
        atomlist_shifted.append( (Z, pos + vec) )
    return atomlist_shifted

def rotate_atoms(atomlist, axis, angle):
    """
    rotate all atoms around an axis
    """
    R = MolCo.rotation_matrix(axis, angle)
    atomlist_rot = []
    for (Z,pos) in atomlist:
        atomlist_rot.append( (Z, np.dot(R,pos)) )
    return atomlist_rot
    
def rotate_around_bond(fragment, atom1, atom2, angle):
    axis = np.array(atom2[1]) - np.array(atom1[1])
    axis /= la.norm(axis)
    # shift center of bond to the origin
    center = 0.5 * (np.array(atom2[1]) + np.array(atom1[1]))
    # rotate and shift back
    fragment_rot = shift_atoms(
        rotate_atoms(
            shift_atoms(fragment, -center),
            axis, angle),
        center)
    return fragment_rot
    
def build_polyfluoerene(distances, angles):
    F1 = XYZ.read_xyz("F1_S0opt.xyz")[0]
    plane = F1[0:21]  # list of atoms that lie in molecular plane fluoerene
    Cl = [F1[21], F1[22],F1[23]]   # left CH2 group 
    Cr = [F1[24], F1[25],F1[26] ]  # right CH2 group
    Hl =             [F1[28]]  # 3rd hydrogen in left CH3 group
    Hr =             [F1[27]]  # 3rd hydrogen in right CH3 group

    # lattice vector along x-axis
    avec = np.array([-1,0,0])
    # build oligomer by adding monomers from left to right
    oligomer = []
    # CH3 group at the left end
    oligomer += Cl + Hl
    dtot = 0.0
    for d,a in zip(distances, angles):
        monomer = plane + Cr
        monomer = rotate_around_bond(monomer, Cl[0], Cr[0], a)
        oligomer += shift_atoms(monomer, dtot*avec)
        dtot += d
    # complete CH3 group at the right end by adding a hydrogen
    oligomer += shift_atoms(Hr, (dtot-d)*avec)

    return oligomer

def scan_angle_F2():
    n = 2
    distances = 2.5 / AtomicData.bohr_to_angs * np.ones(n)
    scan_geometries = []
    for a in np.linspace(0.0, np.pi, 80.0):
        print (a*180.0/np.pi)
        angles = [a/2.0,-a/2.0]
        F2 = build_polyfluoerene(distances, angles)
        scan_geometries.append( F2 )
    XYZ.write_xyz("F2_angle_scan.xyz", scan_geometries)

def eclipsed_oligomers():    
    # build eclipsed oligomers
    for n in [2,3,4,5,6,7,8]:
        distances = 2.5 / AtomicData.bohr_to_angs * np.ones(n)
        angles = np.zeros(n)
        Fn = build_polyfluoerene(distances, angles)
        XYZ.write_xyz("F%d_eclipsed.xyz" % n, [Fn])

def alternately_rotated_oligomers():
    # neighbouring fluorenes are rotated by +-10.5 degrees
    for n in [2,3,4,5,6,7,8]:
        distances = 2.5 / AtomicData.bohr_to_angs * np.ones(n)
        # alternate fluorenes with + or - 10.5 degrees
        angles = 21.0/2.0 * np.pi/180.0 * pow(-1,np.arange(n))
        Fn = build_polyfluoerene(distances, angles)
        XYZ.write_xyz("F%d_rotated.xyz" % n, [Fn])
    
        
if __name__ == "__main__":
    #scan_angle_F2()
    eclipsed_oligomers()
    alternately_rotated_oligomers()
    
