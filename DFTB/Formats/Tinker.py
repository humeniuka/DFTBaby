"""
read Tinker xyz files
"""
from CCS import AtomicData

def read_tinker_xyz_it(filename, units="Angstrom"):
    """
    Parameters:
    ===========
    filename: path to xyz-file
    units: specify the units of coordinates in xyz-file, "Angstrom" or "bohr"

    Returns:
    ========
    iterator to individual structures
    each structure is a tuple with
      - atomlist: list of (Zi,[xi,yi,zi])
      - atomtypes: list of integers encoding type
      - connectivity: list of index lists
           connectivity[2] = [3,4,5]  means that atom 2 is bonded to atoms 3,4 and 5
    """
    assert units in ["bohr", "Angstrom", "hartree/bohr"]
    fh = open(filename)
    # index of the current structure in the file
    igeo = 1
    while 1:
        line = fh.readline()
        if not line: 
            break
        try:
            # first line contains number of atoms in one structure
            nat = int(line.split()[0])
        except ValueError as e:
            print e
            raise Exception("Probably wrong number of atoms in xyz-file '%s'" % filename)
        atomlist = []
        atomtypes = []
        connectivity = []
        for i in xrange(nat):
            line = fh.readline()
#            print line
            words = line.split()
            atid = int(words[0])
            assert atid == i+1
            atno = AtomicData.atomic_number(words[1])
            # cartesian coordinates of the atom
            x,y,z = map(float,words[2:5])
            attype = int(words[5])
            # indeces of atoms connected to the current atom, starting from 0
            atcon = map(lambda i: int(i)-1, words[6:])
            # convert units
            if units == "Angstrom":
                x,y,z = map(lambda c: c/AtomicData.bohr_to_angs, [x,y,z])
            atomlist.append((atno,(x,y,z)))
            atomtypes.append(attype)
            connectivity.append(atcon)
        igeo += 1
        yield atomlist, atomtypes, connectivity

def read_tinker_xyz(filename, units="Angstrom"):
    """
    same as read_tinker_xyz
    """
    return list(read_tinker_xyz_it(filename, units=units))

def txyz2txt(atomlist, atomtypes, connectivities, title="", units="Angstrom"):
    txt = "         %d   (%s units=%s)\n" % (len(atomlist), title, units)
    for i,(Zi, posi) in enumerate(atomlist):
        x,y,z = posi[0], posi[1], posi[2]
        if units == "Angstrom":
            x,y,z = map(lambda c: c*AtomicData.bohr_to_angs, [x,y,z])
        try:
            atname = AtomicData.atom_names[Zi-1]
        except TypeError:
            # leave it as is
            atname = str(Zi)
        # list of connections
        constr = " ".join(map(lambda j: str(j+1), connectivities[i]))
        txt += "%d  %4s    %+10.15f  %+10.15f  %+10.15f   %4s   %s\n" \
                   % (i+1, atname.upper(),x,y,z, atomtypes[i], constr)
    return txt

def write_tinker_xyz(filename, geometries, atomtypes, connectivities, title="uff parameters", units="Angstrom", mode='w'):
    """
    write list of geometries to Tinker xyz-file
    """
    
    fh = open(filename,mode)
    for atomlist in geometries:
        txt = txyz2txt(atomlist, atomtypes, connectivities, title=title, units=units)
        fh.write(txt)
    fh.close()
    return


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print "Usage: python %s <tinker .xyz file>" % sys.argv[0]
        exit(-1)
    xyz_file = sys.argv[1]
    atomlist, atomtypes, connectivities = read_tinker_xyz(xyz_file)[0]
    print atomlist
    print atomtypes
    print connectivities

    write_tinker_xyz("/tmp/tinker.xyz", [atomlist], atomtypes, connectivities)
