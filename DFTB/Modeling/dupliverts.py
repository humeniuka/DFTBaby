#!/usr/bin/env python
"""
Create an arrangement of molecules
by placing fragments at positions specified in an xyz-file.
Example:

  2
  charge=+1
  C60   0.0   0.0   0.0
  C60   0.0   0.0  20.0

will create an xyz-file with the positions of all atoms in two bucky balls
whose center of masses are separated by 20 Ang.
"""
from numpy import array, sqrt, dot, sin, cos, pi, eye
from os.path import expandvars, expanduser, dirname, join, exists
from DFTB import XYZ

def load_fragments(frag_links_file):
    """
    Load xyz-files for each fragment. To map the fragment names to the fragment
    xyz-files a file with the links has to be provided.
    Example for a fragments.lnk file:

       C60  C60.xyz
       Et   ethylene.xyz

    Parameters:
    ===========
    frag_links_file: path that links fragment names with geometry definitions
    
    Returns:
    ========
    fragments: dictionary with e.g. d["C60"] = (geometry of C60)
    """
    frag_dir=dirname(frag_links_file)
    fh = open(frag_links_file)
    fragments = {}
    for l in fh:
        frag_name, frag_xyz_path = l.split()
        if not exists(frag_xyz_path):
            # not an absolute path, try relative path toe frag_dir
            frag_xyz_path=join(frag_dir, frag_xyz_path)
        fragments[frag_name] = XYZ.read_xyz(frag_xyz_path)[0]
    fh.close()
    return fragments

def rotation_matrix(axis,theta):
    """
    compute rotation matrix that rotates a 3D vector around
    an axis k=(kx,ky,kz) by an angle theta

    Parameters:
    ===========
    axis: numpy array [kx,ky,kz]
    theta: rotation angle in radians

    Returns:
    ========
    rotation matri R

    stolen from http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    """
    axis = axis/sqrt(dot(axis,axis))
    a = cos(theta/2)
    b,c,d = -axis*sin(theta/2)
    return array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def read_xyz_rotation_it(filename, units="Angstrom"):
    """
    Read coordinates and rotations of fragments from an xyz-file. The xyz-file may contain
    for each atom in addition to the positions also an axis of rotation and an angle.
    For example:

       2

       A 0.0 0.0 0.0
       B 0.0 0.0 0.0   1 0 0 90

    The second fragment is rotated by 90 deg around the x-axis.

    An iterator to the geometries is returned.
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
    assert units in ["bohr", "Angstrom", "hartree/bohr"]
    fh = open(filename)
    igeo = 1
    while 1:
        line = fh.readline()
        if not line: break
        try:
            nat = int(line.split()[0])
        except ValueError as e:
            print e
            raise Exception("Probably wrong number of atoms in xyz-file '%s'" % filename)
        title = fh.readline()
        fragments = []
        for i in xrange(nat):
            line = fh.readline()
            words = line.split()
            fragname = words[0]
            x,y,z = map(float,words[1:4])
            if units == "Angstrom":
                x,y,z = map(lambda c: c/XYZ.bohr_to_angs, [x,y,z])
            if len(words) == 8:
                # got rotation 
                kx,ky,kz, angle = map(float, words[4:9])
                # translate angle from degree to radians
                angle *= pi/180.0
                fragments.append((fragname,(x,y,z), ((kx,ky,kz),angle)))
            else:
                fragments.append((fragname,(x,y,z), None))
        igeo += 1
        yield fragments

def substitute_fragments(fraglist, fragnames, fragments):
    atomlist = []
    for (fragname, pos_frag, rotation) in fraglist:
        if (rotation != None):
            # rotation in the frame of the fragment
            (axis, angle) = rotation
            R = rotation_matrix(array(axis), angle)
        else:
            R = eye(3)
        # find fragment that is to be placed at pos_frag
        atomlist_frag = fragments[fragname]
        for (Zi, posi) in atomlist_frag:
            # shift rotated atom to center
            atomlist.append( (Zi, list(array(pos_frag) + dot(R, array(posi)))) )
    return atomlist

def dupliverts(xyz_file, frag_links_file, out_xyz_file):
    """
    
    """
    fragments = load_fragments(frag_links_file)
    fragnames = fragments.keys()
    # copy keywords like charge, etc.
    kwds = XYZ.extract_keywords_xyz(xyz_file)
    title=reduce(lambda a,b: a+" "+b, ["%s=%s " % (k,v) for (k,v) in kwds.iteritems()], "")
    # 
    for i,fraglist in enumerate(read_xyz_rotation_it(xyz_file)):
        if i == 0:
            mode='w'
        else:
            mode='a'
        XYZ.write_xyz(out_xyz_file, [substitute_fragments(fraglist, fragnames, fragments)], title=title, mode=mode)
    

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print "Usage: python %s <xyz-file with fragment arrangement> <file with links to fragment xyz-files> <xyz output file>" % sys.argv[0]
        exit(-1)
    dupliverts(sys.argv[1], sys.argv[2], sys.argv[3])
