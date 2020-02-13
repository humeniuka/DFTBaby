#!/usr/bin/env python
"""
remove all molecule outside sphere of radius R from xyz file
"""
import numpy as np
from numpy import linalg as la

from DFTB import AtomicData
from DFTB.Modeling import MolecularCoords

def connectivity(atomlist, Rcov_scale=1.3):
    """
    If atoms i and j are closer than the
    bond length for a Zi-Zj bond 
    (approximated as the sum of the covalent atomc radii x Rcov_scale)
    the connectivity matrx Con[i,j] = 1 otherwise 0
    """
    # find unique atom types
    atom_types = []
    for (Zi, posi) in atomlist:
        atom_types.append(Zi)
    atom_types = list(set(atom_types))
    # estimate bond lengths between all possible dimer combinations
    bond_lengths = {}
    print "Cut-offs for bonding"
    print "===================="
    for a,Za in enumerate(atom_types):
        name_a = AtomicData.atom_names[Za-1]
        Rcova = AtomicData.covalent_radii[name_a]
        for Zb in atom_types[a:]:
            name_b = AtomicData.atom_names[Zb-1]
            Rcovb = AtomicData.covalent_radii[name_b]
            bl = Rcov_scale * (Rcova + Rcovb)
            bond_lengths[(Za,Zb)] = bl
            bond_lengths[(Zb,Za)] = bl
            print "r(%s - %s)   <   %2.5f bohr" % (name_a.ljust(3), name_b.ljust(3), bl)

    print "build connectivity matrix ..."
    Nat = len(atomlist)
    Con = np.zeros((Nat, Nat), dtype=int)
    for i in range(0, Nat):
        Con[i,i] = 1 # every atom is connected to itself
        Zi, posi = atomlist[i]
        for j in range(i+1, Nat):
            Zj, posj = atomlist[j]
            Rij = la.norm(np.array(posi) - np.array(posj))
#            print "Rij = %s <? %s" % (Rij, bond_lengths[(Zi,Zj)]) 
            if Rij < bond_lengths[(Zi,Zj)]:
                Con[i,j] = 1
                Con[j,i] = 1
#    print "Connectivity matrix"
#    print Con
    return Con

def remove_element(elem, l):
    """
    remove element from list
    """
    pass

def find_connected(Con, i, sel=None):
    """
    recursively find all atoms connected to atom with index i
    """
    if sel == None:
        sel = set(range(0, Con.shape[0]))
    sel -= set([i])
    # find indeces of connected atoms
    directly_connected = []
    for j in sel:
        if Con[i,j] == 1:
            directly_connected.append(j)
    directly_connected = set(directly_connected)
    print "directly connected = %s" % directly_connected
    connected = []
    for at in directly_connected:
        sel_at = sel - set([at])
        connected += find_connected(Con, at, sel=sel_at)
        print "      connected = %s" % connected
    connected += directly_connected
    connected = set(connected)
    return connected

def cut_sphere(atomlist, R=20.0):
    """
    remove all atoms outside sphere

    Parameters:
    ===========
    atomlist
    R: radius of sphere in bohr

    Returns:
    ========
    atomlist with atoms inside sphere
    """
    # shift to center of mass
    masses = AtomicData.atomlist2masses(atomlist)
    pos = XYZ.atomlist2vector(atomlist)
    pos_shifted = MolecularCoords.shift_to_com(pos, masses)
    atomlist = XYZ.vector2atomlist(pos_shifted, atomlist)

    Con = connectivity(atomlist)
    print "recursively remove connected atoms..."
    removed = [] # list of removed indeces
    for i,(Zi,posi) in enumerate(atomlist):
        print "i = %s" % i
        if la.norm(posi) > R:
            if not (i in removed):
                print "remove %s%d" % (AtomicData.atom_names[Zi-1], i)
                removed.append(i)
                # remove connect atoms
                connected = find_connected(Con, i)
                print "and connected atoms %s" % connected
                removed += connected
    removed = set(removed)
    cut_atomlist = []
    for i,(Zi,posi) in enumerate(atomlist):
        if i in removed:
            pass
        else:
            cut_atomlist.append( (Zi, posi) )
    return cut_atomlist


if __name__ == "__main__":
    import sys
    from os.path import expandvars, expanduser

    from DFTB import XYZ

    usage = "python %s <in .xyz> <out .xyz> <radius in bohr>\n" % sys.argv[0]
    usage += "  removes all atoms outside sphere of radius R and all connected atoms\n"

    if len(sys.argv) < 4:
        print usage
        exit(-1)

    xyz_in = expandvars(expanduser(sys.argv[1]))
    xyz_out = expandvars(expanduser(sys.argv[2]))
    R = float(sys.argv[3])

    atomlist = XYZ.read_xyz(xyz_in)[0]
    cut_atomlist = cut_sphere(atomlist, R)
    print "Removed %d atoms" % (len(atomlist) - len(cut_atomlist))
    XYZ.write_xyz(xyz_out, [cut_atomlist])
