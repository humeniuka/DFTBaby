# NOT YET PROPERLY TESTED!!!  BETTER TO USE SLATER ORBITALS INSTEAD OF GAUSSIAN ORBITALS
"""
export to Molden format to visualize orbitals.
To this end we convert Kohn-Sham atomic orbitals from DFTB calculation into STO-6G basis. Note
that the plotted orbitals are only a rough approximation as we replace exact atomic orbitals
by a minimal contraction of Gaussians.
"""
from DFTB import AtomicData, XYZ
from DFTB.AtomicData import atom_names
from STO3G_valence import STO3G_valence
from numpy import array, zeros, dot
import numpy as np
from itertools import groupby
#from scipy.linalg import block_diag
from utils import block_diag
from os.path import expandvars, expanduser

def reorder_orbitals(orbs, atomlist, valorbs):
    """
    orbitals in one l-shell are ordered by their m-values, -l,...,0,...,l
    but molden expects orbitals to be in cartesian order.
    The reordering rules are:
            s                 ->       s
        py, pz, px            ->   px, py, pz
        dxy,dyz,dz2,dzx,dx2y2 -> dz2,dzx,dyz,dx2y2,dxy

    Parameters:
    ===========
    orbs: MO coefficients, orbs[:,i] are the coefficients of the i-th orbital
    atomlist: list of tuples (Za,array([xa,ya,za])) with atomic number and atom positions 
    valorbs: dictionary that gives for each atomic number Za a list of the quantum numbers
      (nj,lj,mj) of each valence orbital for the respective atom type.
      e.g. valorbs[1] = [(1,0,0),(1,1,-1),(1,1,0),(1,1,1)] for 2s,2py,2pz,2px of carbon
    """
    def permutation_matrix(new_order):
        """
        find permutation matrix that reorders elements of a vector.
        
        Parameters:
        ===========
        new_order: dictionary with entries   (old_position in vector):(new position in vector)

        e.g.    {0:1, 1:0} produces the swap matrix [[0,1],[1,0]]

        Returns:
        ========
        permutation matrix: 2D numpy array
        """
        dim = len(new_order)
        P = zeros((dim,dim))
        for old,new in new_order.iteritems():
            P[new,old] = 1.0
        return P

    # permutation matrices for s,p and d blocks
    sph2cart_order = {0 : permutation_matrix({0:0}),
                      1 : permutation_matrix({0:1,1:2,2:0}),
                      2 : permutation_matrix({2:0,3:1,1:2,4:3,0:4})}
    permutation_matrices = []
    for (Zi,posi) in atomlist:
        # iterate over all orbitals grouped by angular momentum shell
        for (n,l),lshell_orbs in groupby(valorbs[Zi], lambda nlm: (nlm[0], nlm[1])):
            assert len(list(lshell_orbs)) == 2*l+1
            perm = sph2cart_order[l]
            permutation_matrices.append(sph2cart_order[l])            
##    P = block_diag(*permutation_matrices)
##    orbs_reordered = dot(P,orbs)
    # reorder orbitals
    i = 0
    orbs_reordered = np.copy(orbs)
    for perm in permutation_matrices:
        dim = perm.shape[0]
        orbs_reordered[i:i+dim] = dot(perm, orbs[i:i+dim])
        i += dim

    return orbs_reordered

class MoldenExporter:
    def __init__(self, dftb2, title="DFTB calculation"):
        """
        Parameters:
        ===========
        dftb2: converged instance of DFTB2
        title (optional): title to include in the header
        """
        self.dftb2 = dftb2
        self.txt = "[Molden Format]\n"
        self.txt += "[Title]\n"
        self.txt += "%s\n" % title
        if dftb2 != None:
            # by default Kohn-Sham orbitals are exported
            self.orbs = self.dftb2.getKSCoefficients()
            self.orbe = self.dftb2.getKSEnergies()
            self.f = self.dftb2.getOccupation()
            self.occs = np.where(self.f > 0.0)[0]
            self.virts = np.where(self.f == 0.0)[0]
    def setOrbitals(self, orbs, orbe, f):
        """
        provide different set of orbitals that should be exported

        Parameters
        ----------
        orbs : MO coefficients
        orbe : MO eigenvalues
        f    : occupation numbers
        """
        self.orbs = orbs
        self.orbe = orbe
        self.f = f
        self.occs = np.where(self.f > 0.0)[0]
        self.virts = np.where(self.f == 0.0)[0]
    def Atoms(self):
        """section with nuclear geometry"""
        self.txt += "[Atoms] AU\n"
        for i,(Zi,posi) in enumerate(self.dftb2.getGeometry()):
            self.txt += "%s     %s     %s       %2.7f %2.7f %2.7f\n" \
                % (atom_names[Zi-1], i+1, Zi, posi[0], posi[1], posi[2])
    def GTO(self):
        """section with contracted Gaussian type orbitals"""
        self.txt += "[GTO]\n"
        for i,(Zi,posi) in enumerate(self.dftb2.getGeometry()):
            self.txt += "%s 0\n" % (i+1)
            self.txt += STO3G_valence[atom_names[Zi-1]]
            # empty line
            self.txt += "\n"
    def MO(self, nr_active_occ=100, nr_active_virt=100):
        """Kohn-Sham orbital coefficients"""
        self.txt += "[5D]\n"
        self.txt += "[MO]\n"
        orbs = reorder_orbitals(self.orbs,\
                                self.dftb2.getGeometry(),\
                                self.dftb2.getValorbs())
        # select which orbitals should be exported around the HOMO-LUMO
        active = np.hstack((self.occs[-nr_active_occ:],
                            self.virts[:nr_active_virt]))
        for i in active:
            self.txt += " Sym=  %sa\n" % (i+1)
            self.txt += " Ene=  %2.7f\n" % self.orbe[i]
            self.txt += " Spin= Alpha\n" 
            self.txt += " Occup=%2.7f\n" % self.f[i]
#            j = 0
#            for (Zi,posi) in self.dftb2.getGeometry():
#                for (ni,li,mi) in 
            for j,aoj in enumerate(orbs[:,i]):
                self.txt += "    %s    %2.7f\n" % (j+1,aoj)
    def writeFile(self, filename):
        """write specifications to molden file"""
        fh = open(filename, "w")
        fh.write(self.txt)
        fh.close()
    def export(self, molden_file=None, molden_nr_occ=100, molden_nr_virt=100):
        """
        export geometry, basis and KS coefficients to molden format
        Molden-File.molden_file: Save molecular Kohn-Sham orbitals to this file in the Molden format
        Molden-File.molden_nr_occ: only these highest occupied orbitals will be exported to the Molden file
        Molden-File.molden_nr_virt: only these lowest unoccupied orbitals will be exported to the Molden file
        """
        if molden_file == None:
            return
        self.Atoms()
        self.GTO()
        self.MO(nr_active_occ=molden_nr_occ, nr_active_virt=molden_nr_virt)

        path = expandvars(expanduser(molden_file))
        self.writeFile(path)
        print "Saved molden file to %s" % path

# add sections (geometry. frequencies) individually
class MoldenExporterSectioned(MoldenExporter):
    def __init__(self, dftb2, title=""):
        self.dftb2 = dftb2
        self.txt = "[Molden Format]\n"
        self.txt += "[Title]\n"
        self.txt += "%s\n" % title
    def addVibrations(self, freq_atomlist, freqs, norm_coords, labels=None):
        """
        freqs should be in atomic units
        """
        self.freq_atomlist = freq_atomlist
        if not hasattr(self, "atomlist"):
            self.atomlist = freq_atomlist
        self.freqs = freqs
        self.norm_coords = norm_coords
        self.labels = labels
        # for chained methods
        return self
    #
    def Frequencies(self):
        self.txt += "[FREQ]\n"
        for f in self.freqs:
            # convert to cm^-1
            self.txt += "%s\n" % (f*AtomicData.hartree_to_wavenumbers)
        self.txt += "[FR-COORD]\n"
        for i,(Zi,posi) in enumerate(self.freq_atomlist):
            self.txt += "%s      %2.7f %2.7f %2.7f\n" \
                % (AtomicData.atom_names[Zi-1], posi[0], posi[1], posi[2])
        self.txt += "[FR-NORM-COORD]\n"
        if self.labels==None:
            comments = ["%d" % i for i in range(0, len(self.freqs))]
        else:
            comments = self.labels
        for ivib, vibdisp in enumerate(self.norm_coords):
            datomlist = XYZ.vector2atomlist(vibdisp, self.freq_atomlist)
            # vibrational displacements for each atom
            self.txt += "# %s\n" % comments[ivib]
            self.txt += " vibration %d\n" % ivib
            for i,(Zi,dposi) in enumerate(datomlist):
                self.txt += " %+2.7f %+2.7f %+2.7f\n" \
                    % (dposi[0], dposi[1], dposi[2])
        
    def export(self, molden_file="/tmp/molden.input"):
        self.Atoms()
        if hasattr(self, "freqs"):
            self.Frequencies()
        self.writeFile(molden_file)
        print "Saved molden file to %s" % molden_file
        return self
