import numpy as np
import numpy.linalg as la

from DFTB import AtomicData, XYZ
from DFTB import SphericalHarmonics
from DFTB.SlaterKoster.SKIntegrals import spline_wavefunction

############ Numerically Exact Basis Set ######

def import_pseudo_atom(Z):
    """
    load python module with data for free and confined
    pseudo atoms

    Parameters:
    ===========
    Z: atomic number

    Returns:
    ========
    confined_atom, free_atom
    """
    # load confined pseudo atom
    confined_atom = __import__("DFTB.SlaterKoster.confined_pseudo_atoms.%s" % AtomicData.atom_names[Z-1], fromlist=["orbital_occupation"])
    # load free pseudo atom to get the true onsite energies
    free_atom = __import__("DFTB.SlaterKoster.free_pseudo_atoms.%s" % AtomicData.atom_names[Z-1], fromlist=['nshell','angular_momenta','energies'])
    return confined_atom, free_atom

def load_pseudo_atoms_old(atomlist, confined=True):
    atomtypes = list(set([Zi for (Zi,posi) in atomlist]))
    atomtypes.sort()
    radial_val = {}
    valorbs = {}
    for Zi in atomtypes:
        # load radial orbitals of pseudo atoms
        confined_atom, free_atom = import_pseudo_atom(Zi)
        if confined == True:
            atom = confined_atom
        else:
            atom = free_atom
        # definition of valence shell
        valorbs[Zi] = []
        radial_val[Zi] = []
        for i in atom.valence_orbitals:
            R_spl = spline_wavefunction(atom.r, \
                                        atom.radial_wavefunctions[i])
            n, l = atom.nshell[i], atom.angular_momenta[i]
            for m in range(-l,l+1):
                valorbs[Zi].append((n-1,l,m))
                radial_val[Zi].append(R_spl)
    return valorbs, radial_val

def load_pseudo_atoms(atomlist, confined=True, orbital_set="valence"):
    """
    load radial parts of valence or core orbitals for each atom

    Optional
    ========
    confined    :  if True, the confined orbitals are loaded instead
                   of the free ones
    orbital_set :  load either 'core' or 'valence' orbitals

    Returns
    =======
    qnumbers    :  dictionary with list of quantum numbers (n-1,l,m)
                   of the loaded orbitals for each atom type
    radial_wfn  :  dictionary with list of splines evaluating the
                   radial wavefunctions for each atom type
    """
    atomtypes = list(set([Zi for (Zi,posi) in atomlist]))
    atomtypes.sort()
    # radial wavefunctions of orbitals for each atom
    radial_wfn = {}
    # quantum numbers (n-1,l,m) of orbitals for each atom
    qnumbers = {}
    for Zi in atomtypes:
        # load radial orbitals of pseudo atoms
        confined_atom, free_atom = import_pseudo_atom(Zi)
        if confined == True:
            atom = confined_atom
        else:
            atom = free_atom

        # Which orbitals (core or valence) should be loaded?
        if orbital_set == "core":
            # find indices of core orbitals
            core_orbitals = []
            for i,occ in enumerate(atom.orbital_occupation):
                if (occ > 0) and not (i in atom.valence_orbitals):
                    # Core orbitals are occupied orbitals that
                    # do not form part of the valence set.
                    core_orbitals.append(i)
                    
            orbital_indices = core_orbitals
        elif orbital_set == "valence":
            orbital_indices = atom.valence_orbitals
        else:
            # by default valence orbitals are loaded
            orbital_indices = atom.valence_orbitals
            
        # 
        qnumbers[Zi] = []
        radial_wfn[Zi] = []
        for i in orbital_indices:
            R_spl = spline_wavefunction(atom.r, \
                                        atom.radial_wavefunctions[i])
            n, l = atom.nshell[i], atom.angular_momenta[i]
            for m in range(-l,l+1):
                qnumbers[Zi].append((n-1,l,m))
                radial_wfn[Zi].append(R_spl)
    return qnumbers, radial_wfn

class AtomicBasisFunction:
    def __init__(self, Zi, center, n, l, m, radial_orbital, atom_index):
        self.Zi = Zi
        self.center = center
        self.R_spl = radial_orbital
        self.n = n
        self.l = l
        self.m = m
        # index of atom in the list of atoms to which this basis function belongs
        self.atom_index = atom_index
    def amp(self, x,y,z):
        x0,y0,z0 = self.center
        dx,dy,dz = x-x0, y-y0, z-z0
        # add some offset to avoid division by zero
        dx += 1.0e-15
        dy += 1.0e-15
        dz += 1.0e-15
        #
        r = np.sqrt(dx*dx + dy*dy + dz*dz)
        # The phase factor (-1)^l is chosen such that matrix elements (S and H0)
        # computed numerically between the atomic basis functions agree with
        # those from the Slater-Koster rules, which contain the factor
        #  orbital_parity = (-1)**(l1+l2). (see DFTB/SlaterKoster/SKIntegrals.py)
        A = (-1)**self.l * self.R_spl(r) * SphericalHarmonics.Yreal(self.l,self.m, (dx,dy,dz))
        return A
    def __call__(self, x,y,z):
        """
        basis function can be evaluated 
        as if it were a normal function f(x,y,z)
        """
        return self.amp(x,y,z)
    def radial_wavefunction(self, r):
        return self.R_spl(r)
    def __str__(self):
        x,y,z = self.center
        txt = "<AO atom=%s%d  center=(%.5f,%.5f,%.5f) n=%s l=%s m=%s />" % \
            (AtomicData.atom_names[self.Zi-1], self.atom_index, x,y,z, self.n, self.l, self.m)
        return txt

class AtomicBasisSet:
    def __init__(self, atomlist, confined=True, orbital_set="valence"):
        qnumbers, radial_wfn = load_pseudo_atoms(atomlist, confined=confined, orbital_set=orbital_set)
        self.bfs = []
        for i,(Zi,posi) in enumerate(atomlist):
            for indx,(ni,li,mi) in enumerate(qnumbers[Zi]):
                basis_function = AtomicBasisFunction(Zi, posi, ni,li,mi, radial_wfn[Zi][indx], i)
                self.bfs.append(basis_function)
        #
#        print "Basis Set"
#        print "========="
#        for bf in self.bfs:
#            print bf


############ LCAO orbitals ######

class MolecularBasisFunction:
    def __init__(self, bfsAO, coeffs):
        """
        Parameters
        ==========
        bfsAO       :  list of atomic basis functions
        coeffs      :  MO coefficients
        """
        self.bfsAO = bfsAO
        self.coeffs = coeffs
    def amp(self, x,y,z):
        """ evaluate MO on the grid (x,y,z) """
        mo = 0.0*x
        for c,bf in zip(self.coeffs, self.bfsAO):
            mo += c * bf.amp(x,y,z)
        return mo
    def __call__(self, x,y,z):
        return self.amp(x,y,z)

class LCAOBasisSet:
    def __init__(self, bfsAO, orbs):
        """
        molecular orbital basis is a linear combination of atomic orbitals

           MO(i) = sum_k orbs[k,i] * AO(k)

        Parameters
        ==========
        bfsAO       :  list of atomic basis functions
        orbs        :  MO coefficients, orbs[:,i] are the coefficients
                       of the i-th MO
        """
        nAO,nMO = orbs.shape
        # list of MOs
        self.bfs = []
        for i in range(0, nMO):
            basis_function = MolecularBasisFunction(bfsAO, orbs[:,i])
            self.bfs.append(basis_function)
        
        
        
############ Auxiliary Basis ####

# UGLY: same name as in GammaApproximation
class AuxiliaryBasisFunction:
    def __init__(self, Zi, center, sigma, powers=(0,0,0)):
        self.Zi = Zi
        self.center = center
        self.sigma = sigma
        self.powers = powers
    def amp(self, x,y,z):
        x0,y0,z0 = self.center
        dx,dy,dz = x-x0,y-y0,z-z0
        r2 = dx*dx+dy*dy+dz*dz
        dr = [dx,dy,dz]
        poly = 1.0
        for n,p in enumerate(self.powers):
            poly *= dr[n]**p
        
        if sum(self.powers) == 0:
            # s-like orbital
            N = 1.0/(2.0*np.pi*self.sigma**2)**(1.5)
        elif sum(self.powers) == 1:
            # p-like orbital
            N = 1.0/((2.0*np.pi*self.sigma**2)**(1.5) * self.sigma**2)
        F = N * poly * np.exp(-r2/(2.0*self.sigma**2))
        F[abs(poly) < 1.0e-12] = 0.0
        return F

class AuxiliaryBasisSet:
    def __init__(self, atomlist, hubbard_U):
        sigmas_s = {}
        sigmas_p = {}
        for Z,U in hubbard_U.iteritems():
            sigmas_s[Z] = 1.0/(np.sqrt(np.pi)*U)
            sigmas_p[Z] = 1.0/(6.0*np.sqrt(np.pi)*U**3) # sigmas_s[Z] #1.0/(6.0*np.sqrt(np.pi)*U**3)
        #
        self.bfs = []
        for i,(Zi, posi) in enumerate(atomlist):
            s  = AuxiliaryBasisFunction(Zi, posi, sigmas_s[Zi], powers=(0,0,0))
            px = AuxiliaryBasisFunction(Zi, posi, sigmas_p[Zi], powers=(1,0,0))
            py = AuxiliaryBasisFunction(Zi, posi, sigmas_p[Zi], powers=(0,1,0))
            pz = AuxiliaryBasisFunction(Zi, posi, sigmas_p[Zi], powers=(0,0,1))
            self.bfs += [s,px,py,pz]

############# Gaussian Basis Set #######
# contracted Gaussian type functions as used by Turbomole or Gaussian
# partly stolen from PyQuante

from scipy import special

class ContractedGaussianBasisFunction:
    def __init__(self, Zi, center, l, m, prims):
        """
        spherical Gaussian
          phi(r,th,phi) = sum_i c_i * r^l * exp(-a_i * r^2) * Yreal_lm(theta,phi)

        Parameters:
        ===========
        Zi: atomic number
        center: tuple (xi,yi,zi) in bohr
        l, m: angular quantum numbers
        prims: list of tuples (coefficient, exponent) for each primitve Gaussian
        """
        self.Zi = Zi
        self.center = center
        self.l = l
        self.m = m
        self.prims = prims
        # find normalization factor
        norm2 = 0.0
        p = 3.0/2.0 + self.l
        for ai,ci in self.prims:
            for aj,cj in self.prims:
                norm2 += ci*cj * pow(ai+aj,-p)
        norm2 *= 0.5 * special.gamma(p)
        self.norm = np.sqrt(norm2)
    def radial_wavefunction(self, r):
        """radial part of Gaussian orbitals"""
        R = 0.0*r
        for alpha,coef in self.prims:
            R += coef * np.exp(-alpha*r**2)
        R /= self.norm
        return R
    def amp(self, x,y,z):
        x0,y0,z0 = self.center
        dx,dy,dz = x-x0, y-y0, z-z0
        # add some offset to avoid division by zero
        dx += 1.0e-15
        dy += 1.0e-15
        dz += 1.0e-15
        #
        r = np.sqrt(dx*dx + dy*dy + dz*dz)
        A = self.radial_wavefunction(r) * SphericalHarmonics.Yreal(self.l,self.m, (dx,dy,dz))
        return A

TM_orbital_ordering = {0: [0],
                       1: [1,-1,0],
                       2: [0,1,-1,2,-2]}
    
class GaussianBasisSet:
    def __init__(self, atomlist, basis_data):
        self.atomlist = atomlist
        self.bfs = []
        sym2angmom = {"S": 0, "P": 1, "D": 2, "F": 3}
        for i,(Zi,posi) in enumerate(atomlist):
            bs = basis_data[Zi]
            for sym,prims in bs:
                l = sym2angmom[sym]
                #for m in range(-l,l+1):
                for m in TM_orbital_ordering.get(l,range(-l,l+1)):
                    bf = ContractedGaussianBasisFunction(Zi, posi, l, m, prims)
                    self.bfs.append(bf)
    def getTransformation_GTO_to_DFTB(self):
        """
        find the transformation matrix that maps the MO coefficients in a Gaussian basis
        such as used by Turbomole to the basis of numerically exact atomic orbitals used in DFTB

        Results:
        T: transformation matrix, such that MO(DFTB) = T.MO(GTO), the DFTB orbitals will not be normalized
             and have to be normalized afterwards.
        """
        # create DFTB basis set
        dftb_basis = AtomicBasisSet(self.atomlist)
        m = len(dftb_basis.bfs)  # number of numerical basis functions
        n = len(self.bfs)        # number of GTO basis functions
        # transformation matrix, MOs(dftb) = T.MOs(GTO)
        T = np.zeros((m,n))
        # overlaps between
        # radial grid for overlap calculations
        rs = np.linspace(0.0000000001, 30.0, 3000)
        dr = np.ediff1d(rs, to_end=rs[-1]-rs[-2])
        olaps = {}
        for i,bfi in enumerate(dftb_basis.bfs):
            for j,bfj in enumerate(self.bfs):
                # we only want to relate orbitals on the same atom
                Rij = la.norm( np.array(bfj.center) - np.array(bfi.center) )
                if Rij > 1.0e-10:
                    continue
                #
                assert bfi.Zi == bfj.Zi
                # same angular momentum
                if (bfi.l != bfj.l) or (bfi.m != bfj.m):
                    continue
                Z = bfi.Zi
                l = bfi.l
                # compute overlap
                if olaps.has_key((Z,l)):
                    olap = olaps[(Z,l)]
                else:
                    olap = np.sum( bfi.radial_wavefunction(rs) * bfj.radial_wavefunction(rs) * rs**2 * dr )
                    olaps[(Z,l)] = olap
                
                T[i,j] = olap
        return T

"""
def extract_turbomole_orbitals(tm_dir):
    from DFTB.Formats import Turbomole
    print "Reading geometry and basis from Turbomole directory %s" % tm_dir
    Data = Turbomole.parseTurbomole(tm_dir + "/coord")
    atomlist = Data["coord"]
    Data = Turbomole.parseTurbomole(tm_dir + "/basis")
    basis_data = Data["basis"]
    bs = GaussianBasisSet(atomlist, basis_data)
    T = bs.getTransformation_GTO_to_DFTB()
    # load Turbomole orbitals
    Data = Turbomole.parseTurbomole(tm_dir + "/mos")
    orbe,orbs = Data["scfmo"]
    # transform them to DFTB basis
    orbs_dftb = np.dot(T, orbs)
    # normalize orbitals
    
    
    from DFTB.Scattering.SlakoScattering import save_dyson_orbitals
    XYZ.write_xyz("/tmp/transformed.xyz", [atomlist])
    save_dyson_orbitals("/tmp/transformed.dyson", ["%d" % i for i in range(0, len(orbe))], orbe, orbs_dftb, mode="w")
        
if __name__ == "__main__":
    import sys
    extract_turbomole_orbitals(sys.argv[1])
"""
