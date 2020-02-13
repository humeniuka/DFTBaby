"""
primitive QM/MM with electrostatic embedding (but without polarization of QM part)

 - QM part with DFTB 
 - MM part with UFF as implemented in Gaussian 09
        or with rudimentary implementation of the DREIDING force field (for periodic boundary conditions)

The total system is partioned into an inner region I, which is treated quantum mechanically
and an outer region O, which is described using a force field. The partioning into inner 
and outer region should be done in such a way that both I and O contain closed shell molecules.
The total energy is computed using a subtractive scheme:

E = E_MM(I+O) + E_QM(I) - E_MM(I)                                 
     + sum_(i in O) sum_(i<j in O) (dQ^MM_i * dQ^MM_j)/|Ri-Rj|    outer MM electrostatic energy
     + sum_(i in I) sum_(j in O)   (dQ^MM_i * dQ^MM_j)/|Ri-Rj|    MM electrostatic interaction inner/outer
  = E_QMMM + E_QM(I)

E_MM includes all force field terms except for the non-bonded electrostatic term.
The MM charges are determined using the Qeq algorithm from UFF and remain constant during optimization
or molecular dynamics simulations. The UFF charges are NOT included in the one-particle hamiltonian
of the inner QM region as point charges. Therefore the the QM-part is not polarized by the external charges.
However, the electrostatic interaction between the MM charges in the inner and outer region is taken into account.
Using the DFTB Mulliken charges for the inner region would lead to the wrong analytical gradients, since
dq^QM({Ri}) depends implicitely on the geometry through the MO coefficients. Although probably the 
analytic gradients could be modified to fix this, this task is not an easy one.

"""
import numpy as np

from DFTB.Formats.Gaussian2py import Gaussian # UFF
from DFTB.ForceField.PeriodicForceField import read_force_field, PeriodicForceField  # periodic DREIDING 

class QMMM:
    def __init__(self, atomlist_full, inner_indeces, embedding="electrostatic", pff_file=None, verbose=0):
        """
        Parameters:
        ===========
        atomlist_full: list of tuples (Zi,[xi,yi,zi]) for all atoms (QM + MM)
        inner_indeces: list of the indeces which belong to the QM atoms
        """
        
        assert embedding in ["mechanical", "electrostatic"]
        if "-" in inner_indeces:
            # index list contains ranges such as "9-14"
            inner_indeces = parseAtomTags(inner_indeces)
        else:
            inner_indeces = list(eval(inner_indeces))
        inner_indeces = list(set(inner_indeces))  # remove duplicate indeces
        inner_indeces.sort()
        if verbose > 0:
            print "Indeces of QM atoms:"
            print inner_indeces
            print "number of QM atoms: %d" % len(inner_indeces)
        # counting in the partitioning file starts at 1
        inner_indeces = np.asarray(inner_indeces, dtype=int)-1

        Nat = len(atomlist_full)
        all_indeces = set(range(0, Nat))
        outer_indeces = list(all_indeces - set(inner_indeces))
        outer_indeces.sort()
        self.inner_indeces = inner_indeces
        self.outer_indeces = outer_indeces
        # sort atoms into inner and outer region
        self.atomlist_full = atomlist_full
        self.atomlist_inner = [self.atomlist_full[i] for i in self.inner_indeces]
        self.atomlist_outer = [self.atomlist_full[i] for i in self.outer_indeces]

        if (pff_file == None):
            # 
            # prepare the drivers for the UFF calculations
            # E^MM(I+O)
            self.FF_full = Gaussian.UFF_handler(self.atomlist_full,
                                                embedding=embedding, verbose=verbose, unique_tmp=True)
            # E^MM(I)
            self.FF_inner = Gaussian.UFF_handler(self.atomlist_inner,
                                                 embedding=embedding, verbose=verbose, unique_tmp=True)
        else:
            # load definitions for periodic force field
            if verbose > 0:
                print "periodic MM calculations with DREIDING"
            atomlist_full_ff, atomtypes_full, charges_full, lattice_vectors = read_force_field(pff_file)
            atomtypes_inner = [atomtypes_full[i] for i in self.inner_indeces]
            charges_inner = [charges_full[i] for i in self.inner_indeces]
            if (len(atomlist_full_ff) != Nat):
                raise ValueError("Wrong number of atoms in '%s'. Expected %d atoms but got %d !" \
                                 % (pff_file, Nat, len(atomlist)))
            # prepare drivers for DREIDING calculations
            # E^MM(I+O)
            self.FF_full =  PeriodicForceField(self.atomlist_full, atomtypes_full, charges_full,
                                               lattice_vectors, [], verbose=verbose)
            # E^MM(I)
            self.FF_inner = PeriodicForceField(self.atomlist_inner, atomtypes_inner, charges_inner,
                                               lattice_vectors, [], verbose=verbose)
            
        self.charge = 0
    def getGeometryFull(self):
        return self.atomlist_full
    def partitionGeometry(self, atomlist):
        """
        divides the geometry into an inner and an outer region

        Parameters:
        ===========
        atomlist: atomlist contains the inner and outer region.

        Returns:
        ========
        atomlist_inner: a list of the atoms that should belong to the inner region. 
           This is the atomlist that the DFTB modules will see and work with.
        """
        self.atomlist_full = atomlist

        # sort atoms into inner and outer region
        self.atomlist_inner = [self.atomlist_full[i] for i in self.inner_indeces]
        self.atomlist_outer = [self.atomlist_full[i] for i in self.outer_indeces]
        
        return self.atomlist_inner
    def setCharge(self, nelec, charge):
        if nelec % 2 == 0:
            self.multiplicity = 1
        else:
            self.multiplicity = 2
        self.charge = charge
    def getEnergy(self):
        """
        
        """
        self.FF_full.calc(self.atomlist_full, charge=self.charge, multiplicity=self.multiplicity)
        self.FF_inner.calc(self.atomlist_inner, charge=self.charge, multiplicity=self.multiplicity)
        #
        en_MM_full = self.FF_full.get_MM_Energy()
        en_MM_inner = self.FF_inner.get_MM_Energy()
        # total energy E^QMMM = E^MM(I+O) - E^MM(I)
        en_QMMM = en_MM_full - en_MM_inner
        return en_QMMM
    def getGradientFull(self, grad_QM_inner):
        """
        
        """
        grad_MM_full = self.FF_full.get_MM_Gradient()
        grad_MM_inner = self.FF_inner.get_MM_Gradient()

        grad = grad_MM_full
        # grad E = grad E^MM(I+O) + grad E^QM(I) - grad E^MM(I)
        for i,iI in enumerate(self.inner_indeces):
            grad[3*iI:3*(iI+1)] += grad_QM_inner[3*i:3*(i+1)] - grad_MM_inner[3*i:3*(i+1)]
        return grad

def parseAtomTags(index_str):
    """
    The 'Atom Group Editor' of GaussView produces a list of atom indeces which may contain
    atom ranges (e.g. 9-14) interspersed with a comma-separated list of indeces. This function
    evaluates all ranges and constructs a complete list.
    """
    indeces = []
    parts = index_str.split(",")
    for p in parts:
        if "-" in p:
            start, end = map(int, p.split("-"))
            indeces += range(start,end+1)
        else:
            indeces += [ int(p) ]
    return indeces

