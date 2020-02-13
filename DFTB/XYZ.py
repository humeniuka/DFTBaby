"""\
 adapted from 

 Code for reading/writing Xmol XYZ files.
 
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""
from AtomicData import atom_names, atomic_number, bohr_to_angs, covalent_radii
from numpy import zeros
import numpy as np
from numpy import linalg as la
import os.path

def read_xyz_it(filename, units="Angstrom", fragment_id=atomic_number):
    """
    same as read_xyz, with the difference that an iterator to the geometries
    is returned instead of a list.
    For very large MD trajectories not all geometries can be kept in memory.  
    The iterator returns one geometry at a time which can be processed in a pipeline.

    Parameters:
    ===========
    filename: path to xyz-file
    units: specify the units of coordinates in xyz-file, "Angstrom" or "bohr"
    fragment_id: a function that takes the name of the atom and assigns a
       number to it, usually the atomic number Z.

    Returns:
    ========
    iterator to individual structures
    """
    assert units in ["bohr", "Angstrom", "hartree/bohr", ""]
    fh = open(filename)
    igeo = 1
    while 1:
        line = fh.readline()
        if not line:
            # end of file reached
            break
        words = line.strip().split()
        if words[0] == "Tv":
            # skip lattice vectors
            continue
        try:
            nat = int(words[0])
        except ValueError as e:
            print e
            raise Exception("Probably wrong number of atoms in xyz-file '%s'" % filename)
        # skip title
        title = fh.readline()
        # read coordinates of nat atoms
        atoms = []
        for i in xrange(nat):
            line = fh.readline()
            words = line.split()
            atno = fragment_id(words[0])
            x,y,z = map(float,words[1:4])
            if units == "Angstrom":
                x,y,z = map(lambda c: c/bohr_to_angs, [x,y,z])
            atoms.append((atno,(x,y,z)))
        igeo += 1
        yield atoms

def read_xyz(filename, units="Angstrom", fragment_id=atomic_number):
    """
    read geometries from xyz-file

    Parameters:
    ===========
    filename: path to xyz-file
    units: specify the units of coordinates in xyz-file, "Angstrom" or "bohr"
    fragment_id: a function that takes the name of the atom and assigns a
       number to it, usually the atomic number Z.
    
    Returns:
    ========
    list of structures, each structure is a list of atom numbers and positions 
       [(Z1, (x1,y1,z1)), (Z2,(x2,y2,z2)), ...]
    """
    geometries = []
    for atoms in read_xyz_it(filename, units=units, fragment_id=fragment_id):
        geometries.append(atoms)
    return geometries

def extract_keywords_xyz(filename):
    """
    The title line in an xyz-file may be used to convey additional information 
    (such as total electronic charge) in the form of key-value pairs.

    Example of xyz-file

    2
    charge=+1 
    H  0.0 0.0 0.0
    H  1.0 0.0 0.0

    Parameters:
    ===========
    filname: path to xyz-file

    Returns:
    ========
    dictionary with key-value pairs
    """
    fh = open(filename)
    # skip first line
    fh.readline()
    # read title
    title = fh.readline()
    # find key-value pairs in title line
    kwds = {}
    parts = title.split()
    for p in parts:
        if "=" in p:
            try:
                k,v = p.strip().split("=")
            except ValueError:
                raise Exception("Format of '%s' in title-line not understood. Note that key-value pairs in the title line of an xyz-file must not contain white spaces: for example, use 'charge=+1' instead of 'charge = +1'." % p)
            try:
                kwds[k] = eval(v)
            except (SyntaxError, NameError): # got a string
                kwds[k] = v
#    print "The following keywords were found in %s: %s" % (filename, kwds)
    return kwds

def write_xyz(filename,geometries,title=" ", units="Angstrom", mode='w'):
    """
    write geometries to xyz-file

    Parameters:
    ===========
    geometries: list of geometries, each is a list of tuples
      of type (Zi, (xi,yi,zi))

    Optional:
    =========
    title: string, if a list of strings is provided, the number
           titles should match the number of geometries
    """
    fh = open(filename,mode)
    txt = xyz2txt(geometries, title, units)
    fh.write(txt)
    fh.close()
    return

def xyz2txt(geometries, title="", units="Angstrom"):
    """
    write nuclear geometry to a string in xyz-format.

    Parameters:
    ===========
    geometries: list of geometries, each is a list of tuples
      of type (Zi, (xi,yi,zi))
    
    Optional:
    =========
    title: string, if a list of strings is provided, the number
           titles should match the number of geometries

    Returns:
    ========
    string 
    """
    txt = ""
    for i,atoms in enumerate(geometries):
        if type(title) == list:
            if i < len(title):
                current_title = title[i]
            else:
                current_title = " "
        else:
            current_title = title
        txt += _append_xyz2txt(atoms, current_title, units)
    return txt

def _append_xyz2txt(atoms,title="", units="Angstrom"):
    txt = "%d\n%s\n" % (len(atoms),title)
    for atom in atoms:
        atno,pos = atom
        x,y,z = pos[0], pos[1], pos[2]
        if units == "Angstrom":
            x,y,z = map(lambda c: c*bohr_to_angs, [x,y,z])
        try:
            atname = atom_names[atno-1]
        except TypeError:
            # leave it as is
            atname = atno
        txt += "%4s %10.15f %10.15f %10.15f\n" \
                   % (atname.capitalize(),x,y,z)
    return txt

def update_xyz(filename,atomlist,title="", units="Angstrom"):
    """
    For an existing file the coordinates are updates with the atomic positions given
    in ``atomlist``. The atomic symbols in the original file are retained. In this way
    different atoms of the same type can be distinguished (e.g. C1,C2) and the symbols
    are not overwritten.

    Parameters:
    ===========
    filename: path to existing xyz-file. If the file does not exist an exception is raised.
    atomlist: list of tuples (Zat,[xi,yi,zi]) with new positions in bohr
    """
    # open original file and read symbols
    fh = open(filename)
    nat = int(fh.readline())
    assert nat == len(atomlist), "number of atoms changed in %s (old: %d, new: %d)" % (filename, nat, len(atomlist))
    # read comment
    comment = fh.readline().strip()
    symbols = []
    for i in xrange(nat):
        line = fh.readline()
        words = line.split()
        symbols.append( words[0] )
    fh.close()
    #
    fh = open(filename, "w")
    print>>fh, "%d" % nat
    print>>fh, "%s" % comment
    for i,atom in enumerate(atomlist):
        atno,pos = atom
        x,y,z = pos[0], pos[1], pos[2]
        if units == "Angstrom":
            x,y,z = map(lambda c: c*bohr_to_angs, [x,y,z])
        print>>fh, "%s %10.15f %10.15f %10.15f" \
                   % (symbols[i],x,y,z)
    fh.close()
        
def read_initial_conditions(filename, units="Angstrom", fragment_id=atomic_number):
    """
    read initial positions and velocities from an structure_*.in file. The velocities
    for each atom follow directly after the geometry specification, e.g.

    2
    C   0.0 0.0 0.0
    N   1.5 0.0 0.0
    0.0 0.0 0.0
    -0.3 0.0 0.0
    """
    assert units in ["bohr", "Angstrom", "hartree/bohr", ""]
    fh = open(filename)
    line = fh.readline()
    try:
        nat = int(line.split()[0])
    except ValueError as e:
        print e
        raise Exception("Error while reading initial conditions from '%s'" % filename)
    # read coordinates
    coords = []
    for i in xrange(nat):
        line = fh.readline()
        words = line.split()
        atno = fragment_id(words[0])
        x,y,z = map(float,words[1:4])
        if units == "Angstrom":
            x,y,z = map(lambda c: c/bohr_to_angs, [x,y,z])
        coords.append((atno,(x,y,z)))
    # read velocities
    vels = []
    for i in xrange(nat):
        line = fh.readline()
        words = line.split()
        atno = coords[i][0]
        vx,vy,vz = map(float,words)
        if units == "Angstrom":
            vx,vy,vz = map(lambda c: c/bohr_to_angs, [vx,vy,vz])
        vels.append((atno,(vx,vy,vz)))
    return coords, vels

def atomlist2vector(atomlist):
    """
    convert a list of atom positions [(Z1,(x1,y1,z1)), (Z2,(x2,y2,z2)),...]
    to a numpy array that contains on the positions: [x1,y1,z1,x2,y2,z2,...]
    """
    vec = zeros(3*len(atomlist))
    for i,(Zi,posi) in enumerate(atomlist):
        vec[3*i+0] = posi[0]
        vec[3*i+1] = posi[1]
        vec[3*i+2] = posi[2]
    return vec

def vector2atomlist(vec, ref_atomlist):
    """
    convert a vector [x1,y1,z1,x2,y2,z2,...] to a list of atom positions 
x    with atom types [(Z1,(x1,y1,z1)), (Z2,(x2,y2,z2)),...].
    The atom types are assigned in the same order as in ref_atomlist.
    """
    atomlist = []
    for i,(Zi,posi) in enumerate(ref_atomlist):
        atomlist.append( (Zi, vec[3*i:3*i+3]) )
    return atomlist

def load_pointcharges(point_charges_xyz=None):
    """
    load position and charges of point charge field
    point_charges: path to xyz-file with point charges (use charges instead of atomic numbers in first column of xyz-file)
    """
    if point_charges_xyz != None:
        pc_file = os.path.expanduser(os.path.expandvars(point_charges_xyz))
        point_charges = read_xyz(pc_file)[0]
    else:
        point_charges = []
    return point_charges

def read_charges(chg_file, units="Angstrom"):
    """
    read geometry and transition charges from chg-file (as created by Multiwfn)

    The chg-file should have the following format:

      1st line: number of atoms Nat
      2nd line: comment
      line 3 to Nat+3: five columns with
         ELEMENT  X Y Z  CHARGE

    Parameters
    ----------
    chg_file    :  path to input file

    Optional
    --------
    units       :  units of positions in input, "Angstrom" or "bohr"

    Returns
    -------
    atomlist    : list of Nat tuples (Zat,[x,y,z])
    charges     : numpy array with Nat charges
    """
    assert units in ["Angstrom", "bohr"]
    fh = open(chg_file)
    # read number of atoms
    try:
        Nat = int(fh.readline())
    except ValueError as e:
        print e
        print "First two lines of chg-file should contain number of atoms and a comment!"
        exit(-1)
    # read comments
    comment = fh.readline()
    # read a line for each atom with atomic number Zat,
    # coordinates X Y Z (in Angstrom) and charges
    charges = np.zeros(Nat)
    atomlist = []
    for i in range(0, Nat):
        words = fh.readline().split()
        Zat = atomic_number(words[0])
        x,y,z = map(float, words[1:4])
        if units == "Angstrom":
            x,y,z = map(lambda c: c/bohr_to_angs, [x,y,z])
        q = float(words[4])

        atomlist.append( (Zat, (x,y,z)) )
        charges[i] = q

    return atomlist, charges

def write_charges(filename, atomlist, charges, title=" ", units="Angstrom", mode='w'):
    """
    write atom positions and partial charges to file

    The chg-file will have the following format:

      1st line: number of atoms Nat
      2nd line: comment
      line 3 to Nat+3: five columns with
         ELEMENT  X Y Z  CHARGE
    
    Parameters
    ==========
    filename   :  path to output file
    atomlist   :  list of tuple (Zat,[x,y,z]) with atomic numbers
                  and positions (in bohr)
    charges    :  list atomic monopole charges

    Optional
    ========
    title      :  2nd line
    units      :  output units, 'Angstrom' or 'bohr'
    mode       :  write file ('w) or append to file ('a')
    
    """
    assert units in ["Angstrom", "bohr"]
    assert len(atomlist) == len(charges)
    # 1st and 2nd line
    txt = "%d\n%s\n" % (len(atomlist),title)
    for i,(Zat,(x,y,z)) in enumerate(atomlist):
        if units == "Angstrom":
            x,y,z = map(lambda c: c*bohr_to_angs, [x,y,z])
            atname = atom_names[Zat-1]
        q = charges[i]
        txt += "%3s    %+10.8f  %+10.8f  %+10.8f    %+10.8f\n" \
               % (atname.capitalize(),x,y,z, q)

    fh = open(filename, mode)
    print>>fh, txt
    fh.close()
    
            
def connectivity_matrix(atomlist, search_neighbours=None, thresh=1.3, hydrogen_bonds=False, debug=0):
    """
    compute matrix that shows which atoms are connected by bonds

    C[i,j] = 0    means that atoms i and j are not connected
    C[i,j] = 1    means that they are connected, they are closer than the <thresh>*(bond length between i and j)

    Paramters:
    ==========
    atomlist: list of tuples (Zi,posi) for each atom

    Optional:
    =========
    search_neighbours: search for connected atoms among the 
         <search_neighbours> atoms that are listed right after the current atom. 
         If search_neighbours == None, all atoms are checked.
    thresh: bond lengths can be larger by this factor and are still recognized
    hydrogen_bonds: include hydrogen bonds, too. If hydrogen donor i and hydrogen acceptor
         j are connected through a hydrogen atom k, the elements C[i,k] and C[j,k] are set
         to 1.

    Returns:
    ========
    Con: 2d numpy array with adjacency matrix
    """
    Nat = len(atomlist)
    Con = np.zeros((Nat,Nat), dtype=int)
    if search_neighbours == None:
        search_neighbours = Nat
    for A in range(0, Nat):
        ZA,posA = atomlist[A]
        for B in range(A+1, min(A+search_neighbours, Nat)):
            ZB,posB = atomlist[B]
            RAB = la.norm(np.array(posB) - np.array(posA))
            # approximate bond length by the sum of the covalent radii (in bohr)
            bond_length = covalent_radii[atom_names[ZA-1]] + covalent_radii[atom_names[ZB-1]]
            if RAB < thresh * bond_length:
                Con[A,B] = 1
                Con[B,A] = 1
                
    if (hydrogen_bonds == True):
        # Hydrogen bonds should have
        # 1) donor-acceptor-distances <= 3.5 Angstrom
        # 2) hydrogen-donor-acceptor angles <= 30 degrees
        #
        max_donor_acceptor_distance = 3.5
        max_H_donor_acceptor_angle = 30.0
        # The following atoms can donate or accept a hydrogen bond: O, N
        donor_acceptor_atoms = [8,9] 
        for A in range(0, Nat):  # loop over possible donors or acceptors
            ZA,posA = atomlist[A]
            posA = np.array(posA)
            if not (ZA in donor_acceptor_atoms):  # cannot donate or accept a hydrogen bond
                continue
            for B in range(A+1, min(A+search_neighbours, Nat)): # loop over possible donors or acceptors
                ZB,posB = atomlist[B]
                posB = np.array(posB)
                if not (ZB in donor_acceptor_atoms):
                    continue
                # donor-acceptor-distance
                RAB = la.norm(posB - posA)
                if (RAB*bohr_to_angs > max_donor_acceptor_distance):
                    continue
                for C in range(0, Nat): # loop over hydrogens
                    ZC,posC = atomlist[C]
                    posC = np.array(posC)
                    if ZC != 1: # not a hydrogen atom
                        continue
                    # Which of the atoms A or B is the donor atom? The one that is closer to the hydrogen
                    RAC = la.norm(posC - posA)
                    RBC = la.norm(posC - posB)
                    if RAC < RBC:
                        # A atom is the hydrogen donor
                        r_donor_H = posC - posA
                        r_donor_acceptor = posB - posA
                    else:
                        # B atom is the hydrogen donor
                        r_donor_H = posC - posB
                        r_donor_acceptor = posA - posB
                    
                    # hydrogen-donor-acceptor angle
                    # angle angle(hydrogen --> donor --> acceptor)
                    angle = np.arccos( np.dot(r_donor_H,r_donor_acceptor)/(la.norm(r_donor_H)*la.norm(r_donor_acceptor)) )
                    if angle*180.0/np.pi < max_H_donor_acceptor_angle:
                        # hydrogen bond found
                        # donor/acceptor -- H
                        Con[A,C] = 1
                        Con[C,A] = 1
                        # H -- acceptor/donor
                        Con[B,C] = 1
                        Con[C,B] = 1

                        if debug > 0:
                            print "hydrogen bond %s(%2.d)--H(%2.d)--%s(%2.d)     distance= %8.4f Ang    angle= %8.4f degrees" \
                            % (atom_names[ZA-1].upper(), A+1, C+1, atom_names[ZB-1].upper(), B+1, RAB*bohr_to_angs, angle*180.0/np.pi)

        # 
    return Con

def hydrogens_to_end(atomlist):
    """
    The atomlist is reordered so that all hydrogen atoms come at the end
    """
    hydrogens = []
    others = []
    for (Zi, posi) in atomlist:
        if Zi == 1:
            hydrogens.append( (Zi, posi) )
        else:
            others.append( (Zi, posi) )
    return others+hydrogens
    
__all__ = ["read_xyz", "write_xyz"]

if __name__ == "__main__":
    atomlist = read_xyz("DFTB/test_structures/benzene.xyz")[0]
    vec = atomlist2vector(atomlist)
    atlist = vector2atomlist(vec, atomlist)
    print "atomlist:"
    print atomlist
    print "atlist  :"
    print atlist
