"""
module for reading Gaussian 09 output files
- extract geometry
- forces
and writing Gaussian 09 command files
"""
from numpy import array
from DFTB.AtomicData import bohr_to_angs, atom_names

def read_geometry(log_file):
    """
    read cartesian positions of atoms from Gaussian output file

    Parameters:
    ===========
    log_file: Gaussian output file

    Returns:
    ========
    dictionary with atom positions
    atoms[i] = (Z, array([x,y,z])) contains the atom position of center i in bohr 
    """
    fh = open(log_file, "r")
    atoms = {}
    while True:
        l = fh.readline()
        if l == "":
            break
        if "Input orientation" in l:
            print "Found geometry"
            
            # skip the next three lines
            for i in range(0, 4):
                skip = fh.readline()
            while True:
                l = fh.readline()
                if "--" in l:
                    # geometry specification finished
                    break
                else:
                    center, atnum, attype, x,y,z = l.strip().split()
                    atoms[int(center)-1] = (int(atnum), array(map(float, [x,y,z]))/bohr_to_angs)
            break
    fh.close()
    return atoms

class FormatError(IOError):
    pass

def read_forces(log_file):
    """
    read cartesian forces in hartrees/bohr from Gaussian 09 output

    Parameters:
    ===========
    log_file: gaussian output file

    Returns:
    ========
    dictionary F
    F[i] = (Z, array([ForceX,ForcY,ForceZ])) 
    is the force acting on atomic center i, Z is the atomic number in hartrees/bohr
    """
    fh = open(log_file, "r")
    forces = {}
    while True:
        l = fh.readline()
        if l == "":
            break
        if "Forces (Hartrees/Bohr)" in l:
            print "Found forces"
            
            # skip the next lines
            for i in range(0, 2):
                skip = fh.readline()
                if i == 1:
                    assert "-----" in skip
            while True:
                l = fh.readline()
                if "--" in l:
                    # geometry specification finished
                    break
                else:
                    center, atnum, Fx,Fy,Fz = l.strip().split()
                    forces[int(center)-1] = (int(atnum), array(map(float, [Fx,Fy,Fz])))
            break
    if len(forces) == 0:
        raise FormatError("No force found in gaussian output file %s." % log_file)
    fh.close()
    return forces

def read_total_energy(log_file):
    """
    read total SCF energy from Gaussian 09 log file
    """
    fh = open(log_file, "r")
    lines = fh.readlines()
    fh.close()
    
    enTot = None
    for l in lines:
        if "SCF Done:  E(" in l:
            words = l.split()
            enTot = float(words[4])
    
    if enTot == None:
        raise FormatError("No SCF energy found in gaussian output file %s." % log_file)
    return enTot

def cpu_count():
    """
    find the number of CPU's available. If the environment variable OMP_NUM_THREADS is set, this
    number is used otherwise the total number of CPU's are used.
    """
    import multiprocessing
    import os
    nr_cpus = int(os.environ.get("OMP_NUM_THREADS", multiprocessing.cpu_count()))
    #print "%s CPU's are available" % nr_cpus
    return nr_cpus

def write_input(com_file, atomlist, chk_file="", route="", title="", charge=0, multiplicity=1):
    """
    produce Gaussian 09 command file.

    Parameters:
    ===========
    com_file: path to input file
    atomlist: list of atom types and positions in bohr (Zi,[xi,yi,zi])
    route: "# PBE/6-311G(3df,3pd)++ Force
    """
    com = ""
    com += "%Mem=100MW\n"
    com += "%%Nproc=%d\n" % cpu_count()
    if chk_file != "":
        com += "%%Chk=%s\n" % chk_file
    com += "%s\n\n" % route
    com += "%s\n\n" % title
    com += "%s %s\n" % (charge, multiplicity)
    for (Zi,posi) in atomlist:
        x,y,z = map(lambda pos: pos*bohr_to_angs, posi)
        com += "%s         %2.7f         %2.7f         %2.7f\n" % (atom_names[Zi-1], x,y,z)
    com += "\n"
    fh = open(com_file,"w")
    fh.write(com)
    fh.close()

if __name__ == "__main__":
    import sys
    log_file = sys.argv[1]
    print read_geometry(log_file)
    print read_forces(log_file)
