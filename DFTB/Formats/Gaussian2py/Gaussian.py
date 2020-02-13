"""
module for reading Gaussian 09 output files
- extract geometry
- forces
and writing Gaussian 09 command files
"""
import numpy as np
import numpy.linalg as la

from DFTB.AtomicData import bohr_to_angs, atom_names
from DFTB import XYZ

import Checkpoint

def read_geometry(log_file):
    """
    read cartesian positions of atoms from Gaussian output file

    Parameters:
    ===========
    log_file: Gaussian output file

    Returns:
    ========
    atomlist with tuples (Zi, array([x,y,z]))
    """
    fh = open(log_file, "r")
    atoms = {}
    while True:
        l = fh.readline()
        if l == "":
            break
        if ("Input orientation" in l) or ("Standard orientation" in l):
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
                    atoms[int(center)-1] = (int(atnum), np.array(map(float, [x,y,z]))/bohr_to_angs)
            break
    fh.close()
    # convert to atomlist
    atomlist = atoms.values()
    if len(atomlist) == 0:
        raise Exception("No coordinates found in %s" % log_file)
    return atomlist

def read_geometries_it(log_file):
    """
    read list of geometries from Gaussian output file

    The geometries are returned in the order in which they appear in
    the log file.

    Parameters:
    ===========
    log_file: Gaussian output file

    Returns:
    ========
    iterator:  each item represent one molecular geometry
               as a list of tuples (Zi, array([x,y,z]))
    """
    fh = open(log_file, "r")
    atoms = {}
    while True:
        l = fh.readline()
        if l == "":
            break
        if ("Input orientation" in l) or ("Standard orientation" in l):
            #print "Found geometry"
            
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
                    atoms[int(center)-1] = (int(atnum), np.array(map(float, [x,y,z]))/bohr_to_angs)
            # convert to atomlist
            atomlist = atoms.values()
            yield atomlist
            
    fh.close()

def read_scf_energies_it(log_file):
    """
    extract SCF energies belonging to each geometry.

    Parameters
    ----------
    log_file   :   path to Gaussian log-file

    Returns
    -------
    energy     :   iterator of floats
    """
    fh = open(log_file, "r")
    atoms = {}
    while True:
        l = fh.readline()
        if l == "":
            break
        if ("Input orientation" in l) or ("Standard orientation" in l):
            #print "Found geometry"

            # find energy belonging to this geometry
            energy=None
            while True:
                l = fh.readline()
                if "SCF Done:  E(" in l:
                    words = l.split()
                    energy = float(words[4])
                    break
                if l == "":
                    break
            if energy is None:
                break
            
            yield energy
            
    fh.close()
    
    
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
            #print "Found forces"
            
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
                    forces[int(center)-1] = (int(atnum), np.array(map(float, [Fx,Fy,Fz])))
            break
    if len(forces) == 0:
        raise FormatError("No force found in gaussian output file %s." % log_file)
    fh.close()
    return forces.values()

def read_symmetric_matrix(fh, N):
    Mat = np.zeros((N,N))
    def read_block(fh, jbl):
        l = fh.readline()
        cols = map(int, l.strip().split())
        for i in range(jbl,N):
            l = fh.readline().replace("D", "E")
            parts = l.strip().split()
            row = int(parts[0])
            mijs = map(float, parts[1:])
            for j,m in enumerate(mijs):
                Mat[i,jbl+j] = m
    for jbl in range(0, N, 5):
        print "Read block for columns %s-%s" % (jbl, jbl+5)
        read_block(fh, jbl)
    Mdiag = np.diag(Mat)
    Moff = Mat - np.diag(Mdiag)
    Mat = Moff + Moff.transpose() + np.diag(Mdiag)
    return Mat
        
def read_force_constants(log_file):
    """
    read hessian matrix from Gaussian frequency calculation

    Add iop(7/33=1) to the route section to print out the hessian matrix, e.g.:

      #P AM1 Opt Freq iop(7/33=1)
    """
    fh = open(log_file, "r")
    atomlist = read_geometry(log_file)
    Nat = len(atomlist)
    H = None
    while True:
        l = fh.readline()
        if l == "":
            break
        if "Force constants in Cartesian coordinates:" in l:
            print "Found force constants"
            H = read_symmetric_matrix(fh, 3*Nat)
    fh.close()
    if H == None:
        raise Exception("No force constants found in %s!" % log_file)
    return H

def read_spectrum(log_file):
    """
    read excitation energies from TD-DFT calculation
    """
    fh = open(log_file, "r")
    exc_energies = []
    spins = []
    irreps = []
    oszis = []
    while True:
        l = fh.readline()
        if l == "":
            break
        if "Excited State" in l:
            parts = l.split()
            spin,sym = parts[3].split("-")
            en_eV = float(parts[4]) # excitation energy in eV
            f = float(parts[8][2:]) # oscillator strength
            exc_energies.append( en_eV )
            irreps.append( sym )
            spins.append( spin )
            oszis.append(f)
    fh.close()
    return spins, irreps, exc_energies, oszis

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

def write_input(com_file, atomlist, chk_file="", route="", title="",
                charge=0, multiplicity=1, connectivity=None,
                verbose=0):
    """
    produce Gaussian 09 command file.

    Parameters:
    ===========
    com_file: path to input file
    atomlist: list of atom types and positions in bohr (Zi,[xi,yi,zi])
    route: route section, e.g. '# PBE/6-311G(3df,3pd)++ Force'
    connectivity: list of connected atoms, connectivity[i] contains the list of atoms which are bonded
       to atom i
    """
    com = ""
    com += "%Mem=100MW\n"
    com += "%%Nproc=%d\n" % cpu_count()
    if chk_file != "":
        com += "%%Chk=%s\n" % chk_file
    com += "%s\n\n" % route
    com += "%s\n\n" % title
    com += "%d %d\n" % (charge, multiplicity)
    for i,(Zi,posi) in enumerate(atomlist):
        x,y,z = map(lambda pos: pos*bohr_to_angs, posi)
        com += "%s         %2.7f         %2.7f         %2.7f\n" % (atom_names[Zi-1], x,y,z)
    com += "\n"
    # connectivity
    if connectivity != None:
        assert "connectivity" in route.lower(), \
            "connectivity matrix specified although the option Geom=Connectivity is missing in the route!"
        for i,(Zi,posi) in enumerate(atomlist):
            con_line = " %d  " % (i+1)
            for j in connectivity[i]:
                # atom j is connected to atom i
                con_line += " %d 1.0" % (j+1)
            com += con_line + "\n"
    fh = open(com_file,"w")
    if verbose > 1:
        print "Gaussian Input:"
        print com
    fh.write(com)
    fh.close()
                   
######### Handlers for running calculations #####
import os
from os.path import join
import uuid # used to create unique temporary file names

class UFF_handler:
    def __init__(self, atomlist, embedding="electrostatic", verbose=0, name="uff", unique_tmp=False):
        self.embedding = embedding
        self.verbose = verbose
        self.atomlist = atomlist
        # find the scratch directory where Gaussian writes to
        try:
            scratch_dir = os.environ["GAUSS_SCRDIR"]
        except KeyError as e:
            print "WARNING: Environment variable GAUSS_SCRDIR not set!"
            print "         Check that g09 is installed correctly!"
            #raise e
            scratch_dir="./"
        if unique_tmp == True:
            name += "-" + str(uuid.uuid4())    
        # write input to temporary file
        self.com_file = join(scratch_dir, "%s.gcom" % name)
        self.chk_file = join(scratch_dir, "%s.chk" % name)
        self.fchk_file = join(scratch_dir, "%s.fchk" % name)
        self.log_file = join(scratch_dir, "%s.log" % name)
        # determine the connectivity of the atoms. The connectivity should not change
        # during a dynamics simulation, since this would lead to kinks in the potential
        # energy
        ConMat = XYZ.connectivity_matrix(atomlist)
        self.connectivity = []  # connectivity[i] is a list with atoms bonded to i
        for i in range(0, len(atomlist)):
            self.connectivity.append( list(np.where(ConMat[i,:] != 0)[0]) )
        # first calculation is done only in order to get the
        # correct partial charges from Qeq
        """
        if self.embedding == "electrostatic":
            self.calc(atomlist, do_Qeq=True)
        """
        self.have_charges = False
    def calc(self, atomlist, do_Qeq=False, charge=0, multiplicity=1):
        """do the actual calculation and parse the results"""
        self.atomlist = atomlist
        if self.embedding == "electrostatic" and self.have_charges == False:
            do_Qeq = True
            self.have_charges = True
        if do_Qeq == True:
            # The charge equilibration with Qeq should be performed only once.
            # For repeated calculations (dynamics, optimization) the partial
            # charges are set to 0, and the electrostatic interaction is added
            # later.
            route = "#P UFF=Qeq Force NoSymm Geom=(Connectivity, NoCrowd)"
        else:
            route = "#P UFF Force NoSymm Geom=(Connectivity, NoCrowd)"
        write_input(self.com_file, self.atomlist, chk_file=self.chk_file, route=route,
                    connectivity=self.connectivity, 
                    charge=charge, multiplicity=multiplicity, title="compute energies and gradients")
        # execute g09
        if self.verbose > 0:
            print "Gaussian 09 ..."
            print "  route: %s" % route
        ret = os.system("g09 < %s 2>&1 > %s" % (self.com_file, self.log_file))
        error_msg = "G09 Calculation failed! Check the log-file %s" % self.log_file
        assert ret == 0, error_msg
        # create formatted checkpoint file
        ret = os.system("formchk %s 2>&1 > /dev/null" % self.chk_file)
        assert ret == 0, "ERROR: formchk failed!"
        # parse the checkpoint file and extract:
        #   - the total energy
        #   - the forces
        #   - the partial MM charges from the charge equilibration
        Data = Checkpoint.parseCheckpointFile(self.fchk_file)
        # save results for later use
        self.mm_energy = Data["_Total_Energy"]
        # The gradient does not contain the contribution from electrostatics
        self.mm_gradient = Data["_Cartesian_Gradient"]
        if do_Qeq == True:
            # At the first step, when Qeq is performed
            # gradient and energy contain electrostatic interactions, already
            self.mm_charges = Data["_MM_charges"]
        elif self.embedding == "electrostatic":
            # add the electrostatic energy
            enCoul, gradCoul = self._electrostatics()
            self.mm_energy += enCoul
            self.mm_gradient += gradCoul
    def clean_up(self):
        if self.verbose > 0:
            print "removing temporary Gaussian files"
        # remove temporary files
        for f in [self.com_file, self.chk_file, self.fchk_file, self.log_file]:
            os.remove(f)
    def _electrostatics(self):
        """
        computes the electrostatic interaction between the MM charges,
         enCoul =  - sum_i<j qi*qj/|Ri-Rj|
        and the gradient w/r/t to the positions of the individual atoms
        """
        Nat = len(self.atomlist)
        x = XYZ.atomlist2vector(self.atomlist)
        #
        enCoul = 0.0
        gradCoul = np.zeros(3*Nat)  # gradient
        for i in range(0, Nat):
            for j in range(i+1,Nat):
                Rij_vec = x[3*i:3*(i+1)] - x[3*j:3*(j+1)]
                Rij = la.norm(Rij_vec)
                # contribution to energy
                Vij = self.mm_charges[i]*self.mm_charges[j]/Rij
                enCoul += Vij
                # contribution to the gradient
                eij = Rij_vec/Rij # unit vector
                gij = Vij / Rij * eij
                gradCoul[3*i:3*(i+1)] -= gij
                gradCoul[3*j:3*(j+1)] += gij
                
        if self.verbose > 0:
            print "MM electrostatic energy: %e" % enCoul        
                
        return enCoul, gradCoul
    def get_MM_Energy(self):
        return self.mm_energy
    def get_MM_Gradient(self):
        return self.mm_gradient
    def get_MM_Charges(self):
        return self.mm_charges
        
if __name__ == "__main__":
    import sys
    """
    log_file = sys.argv[1]
    print read_geometry(log_file)
#    print read_forces(log_file)
#    print read_force_constants(log_file)
    """
    xyz_file = sys.argv[1]
    atomlist = XYZ.read_xyz(xyz_file)[0]
    uff_handler = UFF_handler(atomlist)
    print uff_handler.get_MM_Charges()
    uff_handler.calc(atomlist)
    print uff_handler.get_MM_Gradient()
    print uff_handler.get_MM_Energy()

    uff_handler.calc(atomlist)
    print uff_handler.get_MM_Gradient()
    print uff_handler.get_MM_Energy()
    
#    uff_handler.clean_up()
    
