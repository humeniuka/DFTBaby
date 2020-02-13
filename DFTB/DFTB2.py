#!/usr/bin/env python
"""
self-consistent-charge density-functional tight-binding method 
(derived from expansion to 2nd order of DFT energy around a reference density).

Approximations:
  - minimal basis set only includes valence orbitals of confined pseudo-atoms
  - three-center integrals are neglected
  - density profile of charge fluctuations on individual atoms
    is assumed to be Gaussian and spherically symmetric
"""
from numpy import array, zeros, sqrt, log, exp, dot, argsort, where, sum, vstack, pi, mean, trace, copy
import numpy as np
import numpy.linalg as la
from scipy.linalg import eig, eigh, norm, inv, sqrtm
from scipy.special import erf 
from scipy.optimize import fminbound, fmin
from scipy.misc import derivative
import itertools
import multiprocessing
from DFTB import utils
from DFTB.utils import annotated_matrix
from DFTB.SlaterKoster.SKIntegrals import SlakoTransformations, NoSKDipolesException
from DFTB.RepulsivePotential.RepulsivePotential import RepulsivePotential, read_repulsive_potential, load_repulsive_potentials
from DFTB.SlaterKoster import slako_tables
from DFTB.SKMatrixElements import H0andS, DipoleMatrix, OnSiteMatrix
from DFTB import GammaApproximation as GM
from DFTB.AtomicData import atom_names, kBoltzmann, hartree_to_eV, valence_electrons
from DFTB import Parameters
from DFTB.ConvAcceleration import DIIS_80, DIIS_82
from DFTB import XYZ
from DFTB import Periodic
from DFTB.Parallization import Jobs
from DFTB.Timer import GlobalTimer as T
from DFTB.Modeling import MolecularCoords
from DFTB import Symmetry
from DFTB import Mulliken
from DFTB.extensions import mulliken # for fast Mulliken charge calculation
from DFTB.BasisSets import import_pseudo_atom
from DFTB.DensityFitting import DensityFitting
from DFTB.Dispersion import DispersionCorrection
from DFTB.QMMM import QMMM
from DFTB.Cavity import Cavity
from DFTB.DiskMemory import GlobalDiskMemory as Mem

from DFTB import OrbitalLocalization
from DFTB.Molden import MoldenExporter


from os.path import expandvars, expanduser
from collections import defaultdict

# old long_range_radius=3.03

class SCFNotConverged(Exception):
    pass

class SCFChargeInconsistency(Exception):
    pass

class DFTB2(object):
    def __init__(self, atomlist, parameter_set="homegrown",
                 missing_reppots="error", 
                 reppot_paths="[]",
                 point_charges_xyz=None, initial_charge_guess=None,
                 save_converged_charges=None, verbose=1,
                 distance_cutoff=30.0, 
                 long_range_correction=1, long_range_radius=3.03, long_range_T=1.0, long_range_switching="erf", lc_implementation=None,
                 onsite_correction=0,
                 tune_range_radius=0, save_tuning_curve=None,
                 nr_unpaired_electrons=0,
                 use_symmetry=0,
                 hubbard_scaling=1.0,
                 hubbard_override={},
                 reppot_scaling=1.0,
                 fluctuation_functions='Gaussian',
                 mulliken_dipoles=0,
                 dispersion_correction=None,
                 qmmm_partitioning=None, qmmm_embedding="electrostatic", periodic_force_field=None,
                 cavity_radius=0.0, cavity_force_constant=0.0000005,
                 scratch_dir=None,
                 cpks_solver='direct',
                 **opts):
        """
        prepare everything for a DFTB calculation. The nuclear geometry
        specified by the list of atoms is not important and is not stored. 
        At the initialization stage only the atom types are needed to load 
        the required Slater-Koster files and determine the valence orbitals.

        Parameters:
        ===========
        atomlist:  list of tuples (Zi, [xi,yi,zi])
                   Zi is the atomic number of atom i
                   [xi,yi,zi] is the cartesian position vector of atom i
        Parametrization.parameter_set:  
                   "homegrown" - my own parametrization
                   "hotbit"  - hotbit parameters
                   "mio" - mio-0-1 DFTB+ parameters set
        Parametrization.reppot_paths: list of additional search paths for repulsive potentials. Folders are searched for the repulsive potential modules `a_b.py` in the order: 1) any folder specified in `reppot_paths` 2) the default location for repulsive potentials, i.e. 'DFTB/RepulsivePotentials/reppot_tables/'.
        Parametrization.missing_reppots: determines how missing reppots are treated,  
                   "error" - stop program
                   "dummy" - for the missing potential a dummy with Vrep=0 is used. Electronic spectra do not depend on the repulsive potentials, however optimizations and dynamics simulation will be nonsense with dummy potentials.
        Charges.point_charges_xyz: path to file with point charges which interact electrostatically with the molecular partial charges 
        Charges.initial_charge_guess: path to file with initial partial charges, on each line the partial charge for one atom should be given in the same order as in the xyz file
        Charges.save_converged_charges: path to a file where the partial charges of the converged SCF calculation are stored, on each line the partial charge for one atom is written in the same order as in the xyz file. A subsequent calculation can use this as an initial guess via the 'initial_charge_guess' option.
        Output.verbose: controls amount of output, 
             0 - silent
             1 - in each iteration write energies and Mulliken charges
            >1 - print matrix elements and MO coefficients
        Parametrization.distance_cutoff: orbitals whose centers are further away than this cutoff do not interact (in bohr)
        Long-Range-Correction.long_range_correction:
             0 - disabled
             1 - add exact exchange for large distances
        Long-Range-Correction.long_range_switching: "erf" or "erfgau", determines how the exact HF exchange is turned on as a function of the distance, \
"erf" => erf(R_AB/R_lr), \
"erfgau" => erf(R_AB/R_lr) - 2/(sqrt(pi)*R_lr)*R_AB*Exp(-1/3 (R_AB/R_lr)^2)
        Experimental (not tested properly).long_range_T: 
              parameter T in Fermi function, see option long_range_switching
        Long-Range-Correction.lc_implementation: The long-range correction can be implemented in slightly different ways. In order to reproduce the calculations from JCP 143 134120 (2015) set this option to "old"
        Long-Range-Correction.long_range_radius: 
              distance R_lr in bohr where the long range correction is switched on
        Experimental (not tested properly).onsite_correction: 
              add 1-center Coulomb integrals that are neglected in the Mulliken approximation
        Experimental (not tested properly).tune_range_radius: 
              optimize the long range radius such that the ionization potential equals the negative of the HOMO energy
        Experimental (not tested properly).save_tuning_curve: Save a chart of J(R_lr) plotted agains R_lr to this file
        Experimental (not tested properly).mulliken_dipoles: 1 - the electrostatic interaction between monopoles-dipoles and dipoles-dipoles is included in the coulomb part of the hamiltonian, 0 - only monopole-monopole interaction
        Experimental (not tested properly).dispersion_correction: 1 - the Grimme dispersion correction is added to the repulsive potential.
        nr_unpaired_electrons: number of unpaired electrons for open shell singlet states
        use_symmetry: (0 = off, >0 = on) Brings the molecule into standard orientation and tries to classify excited states by the irreducible representations of the symmetry group detected. The assignment only works for non-degenerate irreps. Not all point groups implemented. Note that symmetry adapted orbitals are not used, therefore this option does not speed up the calculation. The symmetry is only determined for the N=use_symmetry lowest excited states.
        Parametrization.hubbard_scaling: 
              scale all Hubbard parameters U by a common factor s, i.e. all U are replaced by s*U
        Parametrization.hubbard_override: 
              override default Hubbard parameters using a dictionary with atom names (lower case) as keys and associated Hubbard parameters (in Hartree) as values, e.g. "dict(n=0.20, o=0.30)" to set the Hubbard parameters for nitrogen and oxygen
        Parametrization.reppot_scaling: 
              scale all repulsive potentials by a common factor f (f < 1.0 makes them less repulsive, f > 1.0 more repulsive)
        Parametrization.fluctuation_functions: Functional form of the spherical charge fluctuations dn(r) on which the gamma-matrix is based, can be 'Gaussian', 'Slater' or 'Numerical' (the gamma-matrix is interpolated using values tabulated for each atom pair).
        QM/MM.qmmm_partitioning: The system is partitioned into an inner region that is treated with DFTB (QM-part) and an outer region treated with the UFF force field (MM-part). This option should contain the indeces of the atoms in the xyz-file that should belong to the QM part. The counting should start at 1. Indeces can be a list, e.g. "[1,2,3]", or a range, e.g. "range(1,4)", or a combination, e.g. "[1,2,3]+range(4,10)". Since the MM atom types are assigned automatically you need to make sure that the geometry in the xyz-file is optimized and reasonable for a force field. The UFF calculations require a working installation of Gaussian 09
        QM/MM.qmmm_embedding: 'mechanical' or 'electrostatic'. The 'mechanical' embedding neglects all electrostatic interactions between MM atoms. In the case of 'electrostatic' embedding, the MM charges are not included into the DFTB-hamiltonian and therefore do not polarize the orbitals, however, the electrostatic interaction between MM charges in the inner and outer region is included. 
        QM/MM.periodic_force_field: If periodic boundary conditions are require, a rudimentary implementation of the DREIDING force field is used instead of UFF. This parameter should contain the path to a file with the force field definition. The file should specify atom position, atom types and lattice vectors. Examples are given in DFTB/ForceField/examples/.
        Cavity.cavity_radius: Sets the cavity radius r0 in bohr. A spherical confining cavity with radius r0 holds the geometry together. The walls push atoms back towards the center. The potential is Vcav(r) = 1/2 * k * mass * (r-r0)^2   for r > r0. If r0==0, the confinement is turned off. The potential is proportional to the mass, so that all particles experience the same acceleration at the boundary of the cavity.
        Cavity.cavity_force_constant: Sets the force constant k [Hartree bohr^-2 mass_e^-1].
        Experimental (not tested properly).scratch_dir: Large arrays can be mapped to scratch files to economize memory usage. These files will be created temporarily in the directory provided in 'scratch_dir'. Functions using this mapping will run much slower
        CPKS.cpks_solver: Determines how the Coupled-perturbed Kohn-Sham (CPKS) equations should be solved - iteratively with a Krylov subspace method ('iterative') or directly ('direct'). The direct solution requires an extremely large amount of memory, whereas the iterative method takes a long time. The CPKS equations are only solved, if gradients of the MO coefficients are required. Gradients of the ground and excited state are computed using Furche's auxiliary function method, so the CPKS equations are not needed in this case.
        """
        self.verbose = verbose
        self.distance_cutoff = distance_cutoff
        self.long_range_correction = long_range_correction
        assert lc_implementation in ["old", None]
        self.lc_implementation = lc_implementation
        self.long_range_radius = long_range_radius
        self.long_range_T = long_range_T
        self.long_range_switching = long_range_switching
        self.onsite_correction = onsite_correction
        self.tune_range_radius = tune_range_radius
        self.save_tuning_curve = save_tuning_curve
        self.save_converged_charges = save_converged_charges
        self.nr_unpaired_electrons = nr_unpaired_electrons
        self.use_symmetry = use_symmetry
        self.hubbard_scaling = hubbard_scaling
        self.hubbard_override = hubbard_override
        self.fluctuation_functions = fluctuation_functions
        self.mulliken_dipoles = mulliken_dipoles
        self.dispersion_correction = dispersion_correction
        if qmmm_partitioning != None:
            self.qmmm = QMMM(atomlist, qmmm_partitioning, embedding=qmmm_embedding,
                             pff_file=periodic_force_field, verbose=verbose)
            # filter the atomlist so that the rest of the DFTB
            # module only sees the QM part
            atomlist = self.qmmm.partitionGeometry(atomlist)
            # now atomlist contains only the QM part
        else:
            self.qmmm = None
        if cavity_radius > 0.0:
            self.cavity = Cavity(r0=cavity_radius, k=cavity_force_constant)
        else:
            self.cavity = None
        # If scratch_dir is different from None, some functions will
        # wrapped so that they use memmap arrays instead of usual
        # numpy-arrays
        global Mem
        Mem.setScratchDirectory(scratch_dir)
        # CPKS
        self.cpks_solver=cpks_solver
        #
        # find unique atom types
        atomtypes = list(set([Zi for (Zi,posi) in atomlist]))
        atomtypes.sort()
        # dictionaries with one entry per unique atom type
        # find quantum numbers (n,l,m) of valence orbitals
        self.valorbs = {}
        # number of valence electrons
        self.ne_val = {}
        # occupation of valence orbitals
        self.valorbs_occupation = {}
        # hubbard parameters
        self.hubbard_U = {}
        # energies of free orbitals by (n,l) quantum numbers
        self.orbital_energies = {}
        # load default U parameters, will be overridden for hotbit
        hubbard_U_byZ = Parameters.get_hubbard_parameters(parameter_set)
        self.parameter_set = parameter_set
        # Additional paths where to look for repulsive potentials
        self.reppot_paths = reppot_paths
        # How should missing repulsive potentials be treated
        self.missing_reppots = missing_reppots
        for Zi in atomtypes:
            if self.verbose > 0:
                print "Loading pseudoatom for element %s" % atom_names[Zi-1].capitalize()
            if parameter_set in ["homegrown", "mio"]:
                self.hubbard_U[Zi] = hubbard_U_byZ[Zi]
                try:
                    atom, free_atom = import_pseudo_atom(Zi)
                except ImportError as e:
                    print e
                    raise Exception("Could not import pseudo atom %s" % atom_names[Zi-1])
            elif parameter_set == "hotbit":
                # load hotbit element description
                from DFTB.SlaterKoster.hotbit_format import load_hotbit_pseudoatom
                import os.path
                dftb_dir = os.path.dirname(os.path.abspath(__file__))
                atom = load_hotbit_pseudoatom("%s/hotbit_parameters/%s.elm" % (dftb_dir, atom_names[Zi-1].capitalize()))
                free_atom = atom
                self.hubbard_U[Zi] = atom.hubbard_U
            else:
                raise Exception("Unknown parameter set %s" % parameter_set)
            
            # scale all Hubbard parameters by a common factor
            self.hubbard_U[Zi] *= self.hubbard_scaling
            # override default Hubbard parameters by those provided on the command line
            # or in the config file
            atname = atom_names[Zi-1]
            if self.hubbard_override.has_key(atname):
                if self.verbose > 0:
                    print "overriding Hubbard parameter for element %s" % atname.capitalize()
                self.hubbard_U[Zi] = self.hubbard_override[atname]
            
            self.valorbs[Zi] = []
            self.valorbs_occupation[Zi] = []
            for i in atom.valence_orbitals:
                n, l = atom.nshell[i], atom.angular_momenta[i]
                for m in range(-l,l+1):
                    self.valorbs[Zi].append((n-1,l,m))
                    self.valorbs_occupation[Zi].append(atom.orbital_occupation[i]/float(2*l+1))
            self.ne_val[Zi] = sum([atom.orbital_occupation[i]  for i in atom.valence_orbitals])
            self.orbital_energies[Zi] = dict([((n-1,l),en) for (n,l,en) in zip(free_atom.nshell,free_atom.angular_momenta, free_atom.energies)])
        # find unique atom pairs and initialize Slater-Koster tables
        atompairs = list(utils.combinations_with_replacement(atomtypes, 2))
        if self.verbose > 0:
            print "unique atom pairs = %s" % [attup for attup in atompairs]
        if parameter_set == "homegrown":
            try:
                slako_tables.load_atompairs(atompairs)
            except ImportError as e:
                print "\nERROR: Could not load Slater-Koster tables for all atom pairs!\n"
                raise e
            try:
                reppot_tables = load_repulsive_potentials(atompairs, missing_reppots=self.missing_reppots, reppot_paths=self.reppot_paths)
                
            except ImportError as e:
                print "\nERROR: Could not load repulsive potentials for all atom pairs!\n"
                print "       If you are only interested in electronic spectra you can"
                print "       disable this error message via the option"
                print "           missing_reppots=dummy                "
                print "       which sets all missing repulsive potentials to Vrep=0.\n"
                raise e
                
        self.SKT = {}
        self.VREP = {}
        for (Zi,Zj) in atompairs:
            assert Zi <= Zj
            # The repulsive potential for the atom pair (Z1,Z2)
            # can be made more or less repulsive by scaling the distance-axis
            # by a factor f=scaling_factors[(Z1,Z2)]
            # f < 1.0  -> less repulsive
            # f > 1.0  -> more repulsive
            scaling_factors = defaultdict(lambda: reppot_scaling, {})
            """
            # scale each potential by a different factor
            scaling_factors = defaultdict(
                lambda: 0.96, # default factor for atom pairs not specified below
                { # make core potentials less repulsive
                    (1,1) : 0.95, # H-H
                    (1,6) : 0.80, # H-C
                    (1,7) : 0.95, # H-N
                    (1,8) : 0.93, # H-O
                    (6,6) : 0.98, # C-C
                    (6,7) : 0.88, # C-N
                    (6,8) : 0.96, # C-O
                })
            """
            if self.verbose > 0:
                print "load parametrization for atom pair %s-%s" % (atom_names[Zi-1], atom_names[Zj-1])
            try:
                if parameter_set == "homegrown":
                    # load precalculated slako table
                    slako_module = slako_tables.atompairs[(Zi,Zj)]
                    # load repulsive potential
                    reppot_module = reppot_tables.atompairs[(Zi,Zj)]
                    # directly load hotbit repulsive potential
                    #import os.path
                    #dftb_dir = os.path.dirname(os.path.abspath(__file__))
                    #parfile = "%s/hotbit_parameters/%s_%s.par" % (dftb_dir,atom_names[Zi-1].capitalize(), atom_names[Zj-1].capitalize())
                    #reppot_module = read_repulsive_potential(parfile)
                elif parameter_set in ["hotbit", "mio"]:
                    from DFTB.SlaterKoster.SKIntegrals import read_slakoformat
                    import os.path
                    dftb_dir = os.path.dirname(os.path.abspath(__file__))
                    if parameter_set == "mio":
                        # or load mio-0-1 slako tables
                        skffile = "%s/dftb-plus_parameters/mio-0-1/%s-%s.skf" % (dftb_dir,atom_names[Zi-1].capitalize(), atom_names[Zj-1].capitalize())
                        slako_module = read_slakoformat(skffile)
                        reppot_module = read_repulsive_potential(skffile)
                    elif parameter_set == "hotbit":
                        # or load hotbit parameters
                        parfile = "%s/hotbit_parameters/%s_%s.par" % (dftb_dir,atom_names[Zi-1].capitalize(), atom_names[Zj-1].capitalize())
                        slako_module = read_slakoformat(parfile)
                        # load repulsive potential from the same file
                        reppot_module = read_repulsive_potential(parfile)
                else:
                    raise Exception("Unknown parameter set %s" % parameter_set)
            except KeyError as e:
                print e
                raise Exception("Could not import Slater-Koster table for atom pair %s-%s" % \
                                    (atom_names[Zi-1], atom_names[Zj-1]))
            self.SKT[(Zi,Zj)] = SlakoTransformations(slako_module)
            self.VREP[(Zi,Zj)] = RepulsivePotential(reppot_module,scaling_factors)
        if self.verbose > 0:
            print "  Unique Atoms"
            print "  ============"
            print "Atom   |   nr. of valence electrons | nr. of valence orbitals |    hubbard U    |  orbital energies"
            print "----------------------------------------------------------------------------------------------------"
            for Zi in atomtypes:
                orbital_energies = list(set([self.orbital_energies[Zi][(ni,li)] for (ni,li,mi) in self.valorbs[Zi]]))
                orbital_energies.sort()
                print atom_names[Zi-1].ljust(3) + 16*" " + str(self.ne_val[Zi]).rjust(3) \
                                   + 25*" " + str(len(self.valorbs[Zi])).rjust(3) \
                                   + 15*" " + ("%+10.8f" % self.hubbard_U[Zi]) + \
                                   + 5*" " + str(orbital_energies) + "\n" 
        # INITIALIZE GAMMA MATRIX
        if self.fluctuation_functions == "Gaussian":
            sigmas = GM.gaussian_decay(self.hubbard_U, self.valorbs, self.long_range_radius, self.long_range_correction, verbose=self.verbose)
            gf = GM.GammaGaussian(sigmas)
        elif self.fluctuation_functions == "Slater":
            sigmas = GM.slater_decay(self.hubbard_U, self.valorbs, self.long_range_radius, self.long_range_correction, verbose=self.verbose)
            gf = GM.GammaSlater(sigmas)
        elif self.fluctuation_functions == "Numerical":
            gf = GM.GammaNumerical(atompairs)
        else:
            raise ValueError("fluctuation_functions can be 'Gaussian', 'Slater' or 'Numerical' but not %s" % self.fluctuation_functions)

        if self.long_range_correction == 1:
            if self.lc_implementation == "old" or self.fluctuation_functions == "Numerical":
                # The screen Coulomb integral is approximated by a switching function
                # times the unscreened integral. This leads to the limit R->0 gamma_lc(R) = 0
                if self.long_range_switching == "erf":
                    sw = GM.ErfSwitching(self.long_range_radius)
                elif self.dftb.long_range_switching == "erfgau":
                    sw = GM.ErfGauSwitching(self.long_range_radius)
                else:
                    raise ValueError("long_range_switching can be 'erf' or 'erfgau' but not %s" % self.long_range_switching)
                gf_lc = GM.GammaLC_approx(gf, sw)
            else:
                # use the analytic integrals for the gamma_lc from
                # Lutsker, Aradi, Niehaus J. Chem. Phys. 143, 184107 (2015)
                # In particular the limit R->0 gamma_lc(R) is not zero
                if self.fluctuation_functions == "Gaussian":
                    gf_lc = GM.GammaGaussianLC(sigmas, self.long_range_radius) 
                elif self.fluctuation_functions == "Slater":
                    gf_lc = GM.GammaSlaterLC(sigmas, self.long_range_radius) 
        else:
            sw = GM.NoSwitching()
            gf_lc = GM.GammaLC_approx(gf, sw)
        # Coulomb gamma matrix
        self.gm = GM.GammaMatrix(gf)
        ###### with multipoles #############
        if self.mulliken_dipoles == 1:
            gf_multipoles = GM.GammaGaussianMultipole(self.hubbard_U)
            self.gm_multipoles = GM.GammaMatrixMultipoles(gf_multipoles)
        ###### density fitting ####
#        self.density_fitter = DensityFitting(atomlist, self.hubbard_U)
        ####################################
        # lc exchange gamma matrix
        self.gm_lc = GM.GammaMatrix(gf_lc)
        #
        ###### dispersion correction ###########
        if self.dispersion_correction == 1:
            print "dispersion correction turned on"
            self.dispersion = DispersionCorrection(atomlist)
        ###
        self.constraints = []
        self.lagrange = []
        self.Vc = None
#        print "OPTS"
#        print opts

        # point charges
        if point_charges_xyz != None: 
            # point charges, in the same format as atomlist, [(Zi,(xi,yi,zi)),...]
            point_charges = XYZ.load_pointcharges(expandvars(expanduser(point_charges_xyz)))
        else:
            point_charges = []
        self.setPointCharges(point_charges)
        # initial charge guess
        if initial_charge_guess != None:
            charge_guess = np.loadtxt(expandvars(expanduser(initial_charge_guess)))
        else:
            charge_guess = None
        self.setChargeGuess(charge_guess)
            
        if (opts.get("constraints_script",None) != None):
            self.loadConstraints(opts["constraints_script"])
        # for printing annotated matrices
        from DFTB.SlaterKoster.slako_transformations import orbital_symbols
        self.orbital_names = []
        self.atom_labels = []
        for i,(Zi,posi) in enumerate(atomlist):
            self.atom_labels.append(atom_names[Zi-1]+"-"+str(i))
            for (nm1i,li,mi) in self.valorbs[Zi]:
                self.orbital_names.append(atom_names[Zi-1] + "-" + str(i) + " " + str(nm1i+1) + orbital_symbols[(li,mi)])
    def count_orbitals(self, atomlist):
        """
        count the number of orbitals on all atoms in atomlist
        """
        dim = 0
        for i,(Zi,posi) in enumerate(atomlist):            
            dim += len(self.valorbs[Zi])
        return dim
    def setCharge(self, charge):
        self.charge = charge
        # total number of valence electrons
        self.Nelec_val = sum(self.q0) - self.charge
        # Check that there is at least one valence electron.
        if self.Nelec_val < 0.5:
            print "\nWARNING: There are no valence electrons in this molecule, \n         charge = %+d,  nr. valence electrons = %s !\n" % (self.charge, self.Nelec_val)
            
    def getSymmetryGroup(self):
        if self.use_symmetry > 0:
            return self.symmetry_group
        else:
            print "Without option use_symmetry > 0, no point group is associated to the current molecular geometry."
            return None
    def setGeometry(self, atomlist, charge=0.0, keep_charge=False, update_symmetry=True):
        """
        specify the nuclear geometry and total charge for which the electronic structure should
        be calculated.

        Parameters:
        ===========
        atomlist:  list of tuples (Zi, [xi,yi,zi])
                   Zi is the atomic number of atom i
                   [xi,yi,zi] is the cartesian position vector of atom i
        charge: total charge of molecule.
        keep_charge: if True ignore the parameter charge and keep the old one if it exists
        update_symmetry: determine molecular symmetry
        """
        #
        if self.qmmm != None:
            atomlist = self.qmmm.partitionGeometry(atomlist)
        if self.use_symmetry > 0 and update_symmetry == True:
            atomlist = MolecularCoords.standard_orientation(atomlist)
            self.symmetry_group = Symmetry.detect_symmetry(atomlist)
            
            # If the molecule is reoriented, vectorial properties such as NACVs are defined relative
            # to the new orientation. When plotting these vectors the standard orientation has to be used,
            # so we need to save it.
            XYZ.write_xyz("standard_orientation.xyz", [atomlist], title=" symmetry=%s" % self.symmetry_group.name())
            print "Molecular geometry in standard orientation has been saved to 'standard_orientation.xyz'"
            
        self.atomlist = atomlist
        self.orbitals_on_atom = [[] for i in self.atomlist] # indeces of orbitals sitting on atom i
        # number of orbitals per atom
        self.orbsPerAtom = [len(self.valorbs[Z]) for (Z,pos) in atomlist]

        self.q0 = [self.ne_val[Zi] for (Zi,posi) in self.atomlist] # number of valence electrons on neutral atom i
        # determine total number of orbitals
        self.dim = 0
        for i,(Zi,posi) in enumerate(self.atomlist):            
            self.orbitals_on_atom[i] = range(self.dim, self.dim + len(self.valorbs[Zi]))
            self.dim += len(self.valorbs[Zi])

        if keep_charge == True and hasattr(self, "charge"):
            pass
        else:
            self.setCharge(charge)
            if self.qmmm != None:
                self.qmmm.setCharge(self.Nelec_val, charge)

        if hasattr(self, "solvent_cavity") and self.solvent_cavity.implicit_solvent != 0:
            self.solvent_cavity.constructSAS(self.atomlist)
            self.solvent_cavity.constructCOSMO()

        if self.verbose > 0:
            import string
            print "  Total number of valence electrons: %s" % self.Nelec_val
            print "  Geometry"
            print "  ========"
            print string.ljust("Atom",15) + string.center("X/bohr",10) + string.center("Y/bohr",10) + string.center("Z/bohr",10)
            print "-"*45
            for (Zi, posi) in self.atomlist:
                print string.ljust(atom_names[Zi-1], 15) \
                    + string.rjust("%.5f" % posi[0],10) \
                    + string.rjust("%.5f" % posi[1],10) \
                    + string.rjust("%.5f" % posi[2],10)
            print ""
    def setPointCharges(self, point_charges):
        """
        point_charges: list of tuples (Zi, [xi,yi,zi]) with point charges
                   which interact electrostatically with the molecular partial charges 
        """
        self.point_charges = point_charges
        if self.verbose > 0 and len(self.point_charges) > 0:
            import string
            print "  Point Charges"
            print "  ============="
            print string.ljust("Charge",15) + string.center("X/bohr",10) + string.center("Y/bohr",10) + string.center("Z/bohr",10)
            print "-"*45
            for (Zi, posi) in point_charges:
                print string.ljust("%.3f" % Zi, 10) \
                    + string.rjust("%.5f" % posi[0],10) \
                    + string.rjust("%.5f" % posi[1],10) \
                    + string.rjust("%.5f" % posi[2],10)
            print ""
    def setChargeGuess(self, charge_guess):
        self.charge_guess = charge_guess
    def _constructH0andS_AB(self, dimA, dimB, atomlistA, atomlistB, AisB=False):
        """
        compute Hamiltonian and overlap matrix elements between two sets of atoms. If the sets
        A and B contain exactly the same structure AisB should be set to True to ensure that
        the diagonal elements of the Hamiltonian are replaced by the correct on-site energies.

        Parameters:
        ===========
        dimA: number of orbitals of all the atoms in set A
        dimB:  ''                                 in set B
        atomlistA, atomlistB: list of (Zi,(xi,yi,zi)) for each atom
        """
        H0 = zeros((dimA,dimB))
        S = zeros((dimA,dimB))
        # iterate over atoms
        mu = 0
        for i,(Zi,posi) in enumerate(atomlistA):
            # iterate over orbitals on center i
            for (ni,li,mi) in self.valorbs[Zi]:
                # iterate over atoms
                nu = 0
                for j,(Zj,posj) in enumerate(atomlistB):
                    # iterate over orbitals on center j
                    for (nj,lj,mj) in self.valorbs[Zj]:
                        if mu == nu and AisB == True:
                            assert Zi == Zj
                            H0[mu,nu] = self.orbital_energies[Zi][(ni,li)] # use the true single particle orbitals energies
                            S[mu,nu] = 1.0 # orbitals are normalized to 1
                        else:
                            if Zi <= Zj:
                                # the first atom given to getHamiltonian() or getOverlap()
                                # has to be always the one with lower atomic number
                                if i == j and AisB == True:
                                    # different orbitals on the same atom should be orthogonal
                                    assert mu != nu
                                    s = 0.0
                                    h0 = 0.0
                                else:
                                    s  = self.SKT[(Zi,Zj)].getOverlap(li,mi,posi, lj,mj,posj)
                                    h0 = self.SKT[(Zi,Zj)].getHamiltonian0(li,mi,posi, lj,mj,posj)
                            else:
                                # swap atoms if Zj > Zi
                                s  = self.SKT[(Zj,Zi)].getOverlap(lj,mj,posj, li,mi,posi)
                                h0 = self.SKT[(Zj,Zi)].getHamiltonian0(lj,mj,posj, li,mi,posi)
                            H0[mu,nu] = h0
                            S[mu,nu] = s
                        nu += 1
                mu += 1
        return S, H0
    def _distance_matrix(self):
        """
        calculate the distances between atoms and the 
        unit vectors pointing from one atom to another
        """
        atomlist = self.atomlist
        Nat = len(atomlist)
        distances = np.zeros((Nat,Nat))
        directions = np.zeros((Nat,Nat,3))
        for i,(Zi,posi) in enumerate(atomlist):
            for j,(Zj,posj) in enumerate(atomlist):
                if i == j:
                    distances[i,j] = 0.0
                    directions[i,j,:] = np.zeros(3)
                elif i < j:
                    R = np.array(posi) - np.array(posj)
                    Rij = la.norm(R)
                    eij = R/Rij
                    distances[i,j] = Rij
                    directions[i,j,:] = eij
                else:
                    distances[i,j] = distances[j,i]
                    directions[i,j,:] = -directions[j,i]
        self.distances = distances
        self.directions = directions
        return distances, directions
    def _proximity_matrix(self):
        """
        orbitals on atoms that are separated by a large distance will have zero
        overlap and zero hamiltonian matrix elements. This function predicts which
        matrix elements should be 0.
        """
        # Mproximity[i,j] == 1 if  |Ri-Rj| < cutoff
        Nat = len(self.atomlist)
        Mproximity = zeros((Nat, Nat), dtype=int)
        for i in range(0, Nat):
            Zi, posi = self.atomlist[i]
            Mproximity[i,i] = 1
            for j in range(i+1, Nat):
                Zj, posj = self.atomlist[j]
                Rij = self.distances[i,j] 
                if Rij < self.distance_cutoff:
                    Mproximity[i,j] = 1
                    Mproximity[j,i] = 1
        self.Mproximity = Mproximity
        return Mproximity
    @T.timer
    def _constructH0andS(self):
        """
        construct the matrix elements of H0 and the overlaps S for a fixed
        nuclear geometry using Slater-Koster transformations. 
        This has to be performed only once, since during the DFTB iterations 
        only the mo coefficients and the partial charges change.
        """
        if self.verbose > 0:
            print "constructing H0 and S from Slater-Koster tables"
        S, H0 = H0andS(self.atomlist, self.valorbs, self.SKT, self.orbital_energies, self.Mproximity)
        # write out matrix elements
        if self.verbose > 1:
            print "Overlaps S of atomic orbitals"
            print "============================="
            print annotated_matrix(S, self.orbital_names, self.orbital_names, colwidth=15)
            print "Matrix elements of H0 between atomic orbitals"
            print "============================================="
            print annotated_matrix(H0, self.orbital_names, self.orbital_names, colwidth=15)

        return S, H0
    @T.timer
    def _constructDipoleMatrix(self):
        Dipole = DipoleMatrix(self.atomlist, self.valorbs, self.SKT, self.Mproximity, self.S)
        """
        # check that dipole matrix is really symmetric
        assert np.sum(abs(Dipole[:,:,0] - Dipole[:,:,0].transpose())) < 1.0e-10
        assert np.sum(abs(Dipole[:,:,1] - Dipole[:,:,1].transpose())) < 1.0e-10
        assert np.sum(abs(Dipole[:,:,2] - Dipole[:,:,2].transpose())) < 1.0e-10
        """
        # write out matrix elements
        from DFTB.SlaterKoster.slako_transformations import orbital_symbols
        self.orbital_names = []
        for i,(Zi,posi) in enumerate(self.atomlist):
            for (nm1i,li,mi) in self.valorbs[Zi]:
                self.orbital_names.append(atom_names[Zi-1] + "-" + str(i) + " " + str(nm1i+1) + orbital_symbols[(li,mi)])
        if self.verbose > 1:
            print "dipole matrix elements of atomic orbitals"
            print "========================================="
            print "X component"
            print annotated_matrix(Dipole[:,:,0], self.orbital_names, self.orbital_names, colwidth=15)
            print "Y component"
            print annotated_matrix(Dipole[:,:,1], self.orbital_names, self.orbital_names, colwidth=15)
            print "Z component"
            print annotated_matrix(Dipole[:,:,2], self.orbital_names, self.orbital_names, colwidth=15)
        return Dipole
    #############################################################
    def _construct_h1(self, gamma, dq):
        """
        construct electrostatic potential due to charge 
        fluctuations
        
        h1_(mu,nu) = 0.5*(e_i + e_j)
        e_i = sum_k gamma_ik(R_ik) dq_k

        NOTE: Since h1_(mu,nu) does not depend on the orbital but 
        only on the atom where the orbital sits, it has block form.
        """
        estatpot = dot(gamma,dq) # electrostatic potential on each atom
                                 # due to partial charging
        self.h1 = zeros((self.dim,self.dim))
        # iterate over atoms
        mu = 0
        for i,(Zi,posi) in enumerate(self.atomlist):
            # iterate over orbitals on center i
            for (ni,li,mi) in self.valorbs[Zi]:
                # iterate over atoms
                nu = 0
                for j,(Zj,posj) in enumerate(self.atomlist):
                    # iterate over orbitals on center j
                    for (nj,lj,mj) in self.valorbs[Zj]:
                        self.h1[mu,nu] = 0.5*(estatpot[i] + estatpot[j])
                        nu += 1
                mu += 1
        #
        if self.verbose > 1:
            print "Electrostatic potential matrix h1"
            print "================================="
            print annotated_matrix(self.h1, self.orbital_names, self.orbital_names, colwidth=15)
        #
        return self.h1
    def _coulomb_hamiltonian(self):
        """
        constructs the part Hcoul of the ground state Kohn-Sham hamiltonian
        due to the electrostatic interaction between partial charges and
        partial dipoles
        """
        h1 = self._construct_h1(self.gamma, self.dq)
        Hcoul = h1*self.S
        """
        ### 
        Hcoul_test = self._coulomb_hamiltonian_alternative()
        err = np.sqrt(np.sum(abs(Hcoul - Hcoul_test)**2))
        assert err < 1.0e-10, "err = %s" % err
        ###
        """
        return Hcoul
    def _coulomb_hamiltonian_multipoles(self):
        """
        constructs the part Hcoul of the ground state Kohn-Sham hamiltonian
        due to the electrostatic interaction between partial charges and
        partial dipoles
        """
        Hcoul = np.zeros(self.S.shape)
        # electrostatic potential on each atom
        # due to partial multipoles
#        estatpot = 0.5*np.dot(self.gamma_multipoles + self.gamma_multipoles.transpose(), self.dqm)
        estatpot = np.dot(self.gamma_multipoles, self.dqm)
        #
        """
        print "estatpot"
        print estatpot
        print "estatpot sym"
        print 0.5*np.dot(self.gamma_multipoles + self.gamma_multipoles.transpose(), self.dqm)
        estatpot_monopoles = dot(self.gamma,self.dq) 
        print "gamma multipoles"
        print 0.5*(self.gamma_multipoles+self.gamma_multipoles.transpose())
        print "gamma monopoles"
        print self.gamma
        print "estatpot multipoles"
        print estatpot
        print "estatpot monopoles"
        print estatpot_monopoles
        err = np.sum(abs(estatpot[::4] - estatpot_monopoles))
        assert err < 1.0e-10, "err = %s" % err
        """
        # iterate over atoms
        mu = 0
        for i,(Zi,posi) in enumerate(self.atomlist):
            # iterate over orbitals on center i
            for (ni,li,mi) in self.valorbs[Zi]:
                # iterate over atoms
                nu = 0
                for j,(Zj,posj) in enumerate(self.atomlist):
                    # iterate over orbitals on center j
                    for (nj,lj,mj) in self.valorbs[Zj]:
                        # monopole
                        Hcoul[mu,nu] += 0.5*self.S[mu,nu]*(estatpot[4*i] + estatpot[4*j])
                        # dipole
                        Hcoul[mu,nu] += 0.5*np.dot(self.D[mu,nu,:] - self.S[mu,nu]*self.Rcc, \
                                               (estatpot[4*i+1:4*(i+1)] + estatpot[4*j+1:4*(j+1)]))
                        nu += 1
                mu += 1
        #
        return Hcoul

    def _construct_h_point_charges(self):
        """
        interaction between external point charges and partial charges leads to
        an additional term  
           hpc_mn * Smn
        in the Hamiltonian
        """
        
        # Vcoul_i = sum_I ZI/|Ri-RI|
        self.Vcoul_pc = np.zeros(len(self.atomlist))
        for i,(Zi,posi) in enumerate(self.atomlist):
            for I,(ZpcI,RI) in enumerate(self.point_charges):
                self.Vcoul_pc[i] += (-ZpcI) / la.norm(np.array(posi)-np.array(RI))
        self.hpc = np.zeros((self.dim,self.dim))
        # hpc_mn = sum_i Vcoul_i * (delta(m in I) + delta(n in I))
        # iterate over atoms
        mu = 0
        for i,(Zi,posi) in enumerate(self.atomlist):
            # iterate over orbitals on center i
            for (ni,li,mi) in self.valorbs[Zi]:
                self.hpc[mu,:] += self.Vcoul_pc[i]
                self.hpc[:,mu] += self.Vcoul_pc[i]
                mu += 1
        #
        if self.verbose > 1:
            print "Electrostatic point charges hpc  "
            print "================================="
            print annotated_matrix(self.hpc, self.orbital_names, self.orbital_names, colwidth=15)
        return self.hpc
    def _density_matrix(self, orbs, f, check=False):
        occ_indx = where(f > 0.0)[0]
        occ_orbs = orbs[:,occ_indx]
        P = dot(f[occ_indx]*occ_orbs, occ_orbs.transpose())
        if check == True:
            err = abs(sum(P * self.S) - self.Nelec_val)
            assert err < 1.0e-10, "err = %s" % err
        return P
    def density_matrix_ref(self):
        """construct reference density matrix"""
        # all atoms should be neutral
        P0diag = []
        # iterate over atoms
        for i,(Zi,posi) in enumerate(self.atomlist):
            # iterate over orbitals on center i
            for iv,(ni,li,mi) in enumerate(self.valorbs[Zi]):
                # how many valence electrons are put into the nl-shell?
                P0diag.append( self.valorbs_occupation[Zi][iv] )
                #        print "reference density matrix"
                #        print P0diag
        self.P0 = np.diag(P0diag)
        return self.P0
    def _constructDensityMatrix(self):
        """
        form density matrix P_mn = sum_a f_a C_ma* C_na
        """
        import warnings
        with warnings.catch_warnings(): # ignore warnings about overflow in exp((en-mu)/(kB*T))
            warnings.simplefilter("ignore")
            mu, self.f = fermi_occupation(self.orbe, self.Nelec_val-self.nr_unpaired_electrons, Nelec_unpaired=self.nr_unpaired_electrons, T=self.temperature)
#        assert sum(self.f) == self.Nelec_val
        self.P = self._density_matrix(self.orbs, self.f)
        return self.P
    def _constructEnDensityMatrix(self):
        """
        form energy weighted density matrix eP_mn = sum_a f_a en_a C_ma* C_na
        """
        occ_indx = where(self.f > 0.0)[0]
        occ_orbs = self.orbs[:,occ_indx]
        # WARNING: if nocc == nvirt, this line might cause a bug because it assumes that
        # orbe is multiplied with the second axis of occ_orbs!!!
        enP = dot(self.f[occ_indx]*self.orbe[occ_indx]*occ_orbs,occ_orbs.transpose())
        return enP
    def _balance_partial_charges(self):
        """
        ensure that partial charges add to the net charge of the molecule
        by distribution the charge mismatch over all atoms 
        """
        partial_net_charge = sum(self.dq)
        # internally charge is measured in units of e-, but self.charge contains the charge in units of (-e)
        charge_mismatch = partial_net_charge + self.charge
        # change partial charges slightly to ensure they sum to the total charge
        self.dq -= charge_mismatch/float(len(self.dq))
    def _center_of_nuclear_charge(self):
        """
        compute the center of the nuclear point charges

           Rcc = (sum_i Zi*Ri)/(sum_i Zi)
        """
        self.Rcc = np.zeros(3)
        nuclear_charge = 0.0
        for i,(Zi,posi) in enumerate(self.atomlist):
            self.Rcc += Zi*np.array(posi)
            nuclear_charge += Zi
        self.Rcc /= nuclear_charge

    @T.timer
    def _MullikenAnalysis(self):
        """
        perform Mulliken population analysis
        """      
        # MONOPOLES
        # faster fortran code
        self.q, self.dq = mulliken.mulliken.monopoles(self.P, self.P0, self.S, self.orbsPerAtom)
        """
        # slow python code
        q_py, dq_py = Mulliken.monopoles(self.atomlist, self.P, self.P0, self.S, self.valorbs)
        assert np.sum(abs(self.q-q_py)) < 1.0e-10
        assert np.sum(abs(self.dq-dq_py)) < 1.0e-10
        """
        #
        if abs( np.sum(self.dq) + self.charge ) > 1.0e-5:
            raise SCFChargeInconsistency("Mulliken charges do not add up to total charge, sum dq = %s != total charge = %s" % (-np.sum(self.dq), self.charge))
        
        if self.mulliken_dipoles == 1:
            # DIPOLES
            self.dip, self.ddip = Mulliken.dipoles(self.atomlist, self.P, self.P0, self.S, self.D, self.valorbs,
                                                   self.Rcc)
            # MULTIPOLES
            self.qm, self.dqm = Mulliken.multipoles(self.atomlist, self.P, self.P0, self.S, self.D, self.valorbs,
                                                    self.Rcc)
            #
            # MULTIPOLES from fitting
#            # overwrite Mulliken dipoles
#            self.dqm = self.density_fitter.fit(self.P-self.P0, np.sum(self.dq), np.sum(self.ddip, axis=0))
#            self.dq = self.dqm[::4]
#            self.ddip[:,0] = self.dqm[1::4]
#            self.ddip[:,1] = self.dqm[2::4]
#            self.ddip[:,2] = self.dqm[3::4]            
            #            # save Mulliken charges for each SCF step
#            Mulliken.save_partial_dipoles("/tmp/mulliken_dipoles_%.4d.dat" % self.i, self.atomlist, self.ddip)

        """
        self.q = zeros(len(self.atomlist)) # Mulliken charges
        q_mu = sum(self.P*self.S, axis=1) # charge in orbital mu
        ##
#        print "Q_MU"
#        for i,(Zi,posi) in enumerate(self.atomlist):
#            print "%s-%s = %s" % (Zi, i, q_mu[self.orbitals_on_atom[i]])
#        print "sum(q_mu) = %s" % sum(q_mu)
        ##
        # sum over all orbitals mu belonging to atom I
        for i in range(0, len(self.atomlist)):
            self.q[i] = sum(q_mu[self.orbitals_on_atom[i]])
            self.dq[i] = self.q[i] - self.q0[i] # excess charge
        assert abs( sum(self.dq) - (sum(q_mu) - sum(self.q0)) ) < 1.0e-7
#        print "sum of partial charges = %s e-" % sum(self.dq)
#        self._balance_partial_charges() # WHY IS THIS NEEDED?
        if self.i > 0:
            assert abs(sum(self.dq) + self.charge) < 1.0e-4, "sum(dq) = %s   charge = %s" % (sum(self.dq), self.charge)
        """
    def runNonSCC(self, temperature=0.0):
        self.temperature = temperature
        self._distance_matrix()
        self._proximity_matrix()
        self.S, self.H0 = self._constructH0andS()
        if len(self.point_charges) > 0:
            hpc = self._construct_h_point_charges()
            H = self.H0 + hpc*self.S
        else:
            H = self.H0
        self.orbe, self.orbs = eigh(H,self.S)
        self._constructDensityMatrix()
        self.HLgap = self.getHOMO_LUMO_gap()
        # dq should be zero anyway 
        self.gamma = self.gm.gamma_atomwise(self.atomlist, self.distances)[0]
        
        self.q = self.q0
        self.dq = zeros(len(self.atomlist))
        self.getEnergies()
        # write_iteration expects these member variables to be set
        self.i = 0
        self.relative_change = 0.0
        self.writeIteration()
    def runSCC(self, maxiter=1000, \
                   scf_conv=1.0e-10, \
                   start_from_previous=1, \
                   mixing_threshold=1.0e-3, \
#                   density_mixer=None, fock_interpolator=DIIS_82(5),
#                   density_mixer=None, fock_interpolator=None, \
#                   mixing_threshold=0.2, \
                   density_mixer=DIIS_80(5),  fock_interpolator=None, \
                   level_shift=0.1, HOMO_LUMO_tol=0.05, \
                   linear_mixing_coefficient=0.33):
#                   level_shift=0.5, HOMO_LUMO_tol=0.5, linear_mixing_coefficient=0.5):
        """
        run a self-consistent-charge calculation

        Parameters:
        ===========
        SCF-Convergence.maxiter: stop SCF calculation after maxiter iterations
        SCF-Convergence.scf_conv: convergence threshold for relative change in SCF-calculation
        SCF-Convergence.density_mixer: instance of Mixer object, that determines how the next density
          matrix is obtained from previous guesses, if None a linear mixer is used
        SCF-Convergence.mixing_threshold: if the relative change drops below this value density mixing is used
        SCF-Convergence.fock_interpolator: the KS hamiltonian is interpolated based on previous iteration steps
        SCF-Convergence.level_shift: shift virtual orbitals up in energy, this shift parameter is gradually reduced to zero as the density matrix converges
        SCF-Convergence.HOMO_LUMO_tol: level shifting is turned on, as soon as the HOMO-LUMO gap drops below this value
        SCF-Convergence.linear_mixing_coefficient: is no density mixer object is used (density_mixer=None) the next guess for the density matrix is constructed as P_next = a*P + (1-a)*P_last
        SCF-Convergence.start_from_previous: use density matrix and partial charges from a previous calculation as initial guess if available (1) or start from reference density and zero partial charges (0)
        """
        if type(density_mixer) == type(""):
            # passed from command line, could be a security problem
            density_mixer = eval(density_mixer)
        if type(fock_interpolator) == type(""):
            # passed from command line
            fock_interpolator = eval(fock_interpolator)
        self.maxiter = maxiter
        self.temperature = 0.0 # occupation of orbitals is smeared out by Fermi distribution with temperature T in Kelvin
        # guess initial partial charges {dq_j}
        if self.charge_guess != None:
            Nat = len(self.atomlist)
            self.dq = self.charge_guess
            assert len(self.dq) == Nat # need partial charge for each atom
            # dipoles
            self.ddip = np.zeros((Nat,3))
            # multipoles (monopoles+dipoles)
            self.dqm = np.zeros(4*Nat)
            self.dqm[::4] = self.dq
            if self.verbose > 0:
                print "sum of partial charges = %s" % sum(self.dq)
            self._balance_partial_charges()
            assert abs(sum(self.dq) - self.charge) < 1.0e-10 # partial charges must sum to total charge
        else:
            if not hasattr(self, "dq") or start_from_previous == 0:
                Nat = len(self.atomlist)
                self.dq = np.zeros(Nat)
                # dipoles
                self.ddip = np.zeros((Nat,3))
                # multipoles (monopoles+dipoles)
                self.dqm = np.zeros(4*Nat)
            else:
                # otherwise use previous charges
                if self.verbose > 0:
                    print "keep partial charges from last calculation as initial guess"
        if density_mixer != None:
            density_mixer.reset()
        if fock_interpolator != None:
            fock_interpolator.reset()
            # compute A = S^(-1/2)
            # 1. diagonalize S
            W,V = eigh(self.S)
            # 2. compute inverse square root of the eigenvalues
            W12 = np.diag(pow(W, -0.5))
            # 3. and transform back
            A = dot(V, dot(W12, V.transpose()))
            print "S*S^(-1)" 
            print dot(dot(A,A), self.S)
        if not hasattr(self,"P") or start_from_previous == 0:
            # initialize density matrix
            self.P = self.density_matrix_ref()
        else:
            if self.verbose > 0:
                print "keep density matrix from last calculation as initial guess"
        P_last = self.P
        converged = False
        shift_flag = False
        mixing_flag = False
        for self.i in range(0, self.maxiter):
            # build hamiltonian 
            if self.mulliken_dipoles == 1:
                Hcoul = self._coulomb_hamiltonian_multipoles()
            else:
                Hcoul = self._coulomb_hamiltonian()
            # point charges
            if len(self.point_charges) > 0:
                hpc = self._construct_h_point_charges()
                Hcoul += hpc*self.S
#            print "h1 = %s" % h1
#            print "dq = %s" % self.dq
            H = self.H0 + Hcoul
            if self.onsite_correction == 1:
                H_onsite = self._onsite_hamiltonian()
                H += H_onsite
            # long range Hartree-Fock exchange
            if self.long_range_correction == 1:
                Hx = self.lc_exact_exchange()
                H += Hx
            # solve generalized eigenvalue problem
            #   H C = e S C
            # to get Kohn-Sham energies orbe[i] and 
            # corresponding orbital coefficients orbs[:,i]
            ##############
            # level shifting
            b = level_shift
            if hasattr(self, 'orbe'):
                if self.HLgap < HOMO_LUMO_tol and mixing_flag == False:
                    if self.verbose > 0:
                        print "gap between HOMO and LUMO fell below tolerance => turn on level shifting"
                    shift_flag = True
            if shift_flag == True and b > 0.0:
                # adjust level shift
                if self.relative_change < 1.0e5 * scf_conv:
                    b *= 0.5 # slowly turn level shift off
                if self.relative_change < 1.0e3 * scf_conv:
                    b = 0.0  # turn level shifting off completely
                if self.verbose > 0:
                    print "level shift b = %s" % b
            if shift_flag == True and b > 0.0:
                norb = self.orbs.shape[0]
                sort_indx = argsort(self.orbe)
                # transform KS matrix into the basis of the last
                # MO coefficients
                B = self.orbs.copy()
                H = dot(B.transpose(), dot(H, B))
                # shift diagonals of the virtual-virtual block in H up by b
                vvBlock = zeros(H.shape)
                virts = sort_indx[where(self.f == 0.0)[0]]
                vvBlock[(virts,virts)] = 1.0
                Hshift = H + b * vvBlock
                # 
                if fock_interpolator != None and self.i > 0 and converged == False:
                    Hshift = fock_interpolator.interpolate_hamiltonian(Hshift)
                self.orbe,orbs = eigh(Hshift)
                # transform orbitals from orthonormal basis back into ao basis
                self.orbs = dot(B, orbs)
            else:
                if fock_interpolator != None and self.i > 0 and converged == False:
                    H = fock_interpolator.interpolate_hamiltonian(H)
                try:
                    self.orbe, self.orbs = eigh(H,self.S)
                except la.linalg.LinAlgError as e:
                    print "WARNING: %s" % e
                    print "Trying to solve eigenvalue problem H'.C' = C'.e after Loewdin orthogonalization."
                    # There seems to be a bug in scipy.linalg.eigh
                    
                    # convert generalized eigenvalue problem H.C = S.C.e into eigenvalue problem H'.C' = C'.e
                    # by Loewding orthogonalization, H' = X^T.H.X, where X = S^(-1/2)
                    X = la.inv(sqrtm(self.S))     # X = S^(-1/2)
                    # H' = X^t.H.X
                    Hp = np.dot(X.conjugate().transpose(), np.dot(H, X))
                    self.orbe, Cp = la.eig(H)
                    # C = X.C'
                    self.orbs = np.dot(X, Cp)
        #
##########################
#            print self.orbe
            self._constructDensityMatrix()
            # use DIIS or other mixing scheme to speed up convergence
            if density_mixer != None and converged == False and mixing_flag == True:
                next_P = density_mixer.next_approximation(self.P)
                self.P = next_P
            else:
                if hasattr(self, "P") and converged == False:
                    # linear damping
                    alpha = linear_mixing_coefficient
                    next_P = alpha * self.P + (1.0-alpha)*P_last
                    # normalize P so that Tr(R*S) = Nelec
                    next_P *= sum(self.P * self.S) / sum(next_P * self.S) 
                    self.P = next_P

            # update partial charges using Mulliken analysis
            self._MullikenAnalysis()
            # compute energy
            self.getEnergies()
            # check for convergence
            # Does the density matrix commute with the KS Hamiltonian?
            # err = H*D*S - S*D*H
            if fock_interpolator != None and self.i > 0 and converged == False:
                comm_err = dot(H, dot(self.P, self.S)) - dot(dot(self.S, self.P), H)
                # transform error vector to orthogonal basis
                comm_err = dot(A.transpose(), dot(comm_err, A))
                fock_interpolator.set_error_vector(comm_err)
                print "Estimate SCF convergence from the commutator |H*P*S - S*P*H| = %s" % norm(comm_err)
            #
            self.HLgap = self.getHOMO_LUMO_gap()
            if converged == False:
                if density_mixer != None and mixing_flag == True:
                    self.relative_change = density_mixer.relative_change()
                else:
                    self.relative_change = norm(P_last - self.P)/norm(self.P)
            if self.verbose > 0:
                self.writeIteration()
            if converged == True:
                if self.verbose > 0:
                    print "!!! CONVERGED after %s iterations (relative change = %.2e < %.2e = threshold)!!!" \
                    % (self.i+1, self.relative_change, scf_conv)
                break

            if self.relative_change < scf_conv:
                # perform a last iteration without damping and extrapolation
                # so as to obtain a variational ground state
                converged = True
            if self.relative_change < mixing_threshold and mixing_flag == False:
                if density_mixer != None:
                    # Turn on density mixing for slightly converged solution
                    mixing_flag = True
                    if self.verbose > 0:
                        print "density mixing is turned on"
            P_last = self.P
        else:
            raise SCFNotConverged("SCF calculation did not converge after %s iterations!!! You should try to play around with the option --linear_mixing_coefficient" % (self.i+1))
        if self.save_converged_charges != None:
            # overwrite initial charge guess
            np.savetxt(expandvars(expanduser(self.save_converged_charges)), self.dq)

#######################
# On-Site Correction
#######################
    def _load_onsite_integrals(self):
        self.G_onsite = OnSiteMatrix(self.atomlist, self.valorbs)
        # write out matrix elements
        if self.verbose > 1:
            print "1-center Coulomb integrals (mn|mn)"
            print "=================================="
            print annotated_matrix(self.G_onsite, self.orbital_names, self.orbital_names, colwidth=15)

    def _onsite_hamiltonian(self):
        """ addition to Hamiltonian due to on-site correction """
        dP = self.P - self.P0
        H_onsite = dP * self.G_onsite
        return H_onsite
    def _onsite_energy(self):
        """ additional energy due to on-site correction """
        dP = self.P - self.P0
        E_onsite = 0.5 * np.sum(dP*self.G_onsite*dP)
        return E_onsite
        
#######################
# Long range correction
#######################
    def _switching_function(self):
        if self.long_range_switching == "erf":
            R_lr = self.long_range_radius
            def switch(Rij):
                return erf(Rij/R_lr)
            def switch_deriv(Rij):
                return 2.0/np.sqrt(np.pi) * (np.exp(-pow(Rij/Rlr,2)) / Rlr)
        elif self.long_range_switching == "fermi":
            R0 = self.long_range_radius
            T  = self.long_range_T
            def switch(Rij):
                return 1.0/(np.exp(-(Rij-R0)/T)+1.0)
            def switch_deriv(Rij):
                return np.exp((R0-Rij)/T)/(T*pow(1.0+exp((R0-Rij)/T),2))
        elif self.long_range_switching == "erfgau":
            # see Phys. Rev. A (2004), 70, 062505
            # Toulouse, Colonna, Savin, "Long-range short-range separation of the electron-electron interaction in density functional theory"
            R0 = self.long_range_radius
            def switch(Rij):
                if Rij < 1.0e-8:
                    return 0.0
                sw = erf(Rij/R0)/Rij - 2.0/(np.sqrt(np.pi)*R0)*np.exp(-1.0/3.0 * (Rij/R0)**2)
                return sw
            def switch_deriv(Rij):
                r2 = pow(Rij/R0,2)
                sw_deriv = 4.0/(3.0*np.sqrt(np.pi)*R0**3) * np.exp(-r2/3.0) \
                    +2.0/np.sqrt(np.pi*R0) * np.exp(-r2)/Rij \
                    -erf(Rij/R0)/Rij**2
                return sw_deriv
        else:
            raise ValueError("--long_range_switching can be 'erf', 'fermi' or 'erfgau' but not '%s'" % self.long_range_switching)
        return switch, switch_deriv

    def _coulomb_hamiltonian_alternative(self):
        """
        compute the Hamiltonian Hcoul as

          Hcoul_mn = 1/2 S_mn * sum_k ( G_mk + Gnk ) sum_l (P - P0)_kl Skl
        
        This function returns the same matrix as `_coulomb_hamiltonian()`
        """
        ret = self.gm.gamma_AOwise(self.atomlist, self.valorbs, self.distances, self.directions)
        self.G = ret[2]
        
        dP = self.P - self.P0
        dPS = np.sum(dP*self.S, axis=1)
        GdPS = np.dot(self.G, dPS)
        Hcoul = 0.5 * self.S * np.add.outer(GdPS, GdPS)
        return Hcoul

    def getGamma_lr(self):
        return self.gamma_lr
    def lc_exact_exchange(self):
        """
        construct part of the Hamiltonian matrix corresponding to long range 
        Hartree-Fock exchange
          H^x_mn = -1/2 sum_ab (P_ab-P0_ab) (ma|bn)_lr
        The Coulomb potential in the electron integral is replaced by 
            1/r ----> erf(r/R_lr)/r
        """
        if self.verbose > 1:
            print "Compute long range HF exchange"
        Hx = zeros((self.dim, self.dim))
        Nat = len(self.atomlist)
        if self.lc_implementation == "old":
            # The exact exchange is added for the full density matrix
            # including the reference density
            dP = self.P
        else:
            # The exact exchange is added only for the difference density
            # as done in J. Chem. Phys. 143, 184107 (2015), where it is assumed
            # that the exchange energy due to the reference density P0
            # is already contained in the parametrization.
            dP = self.P-self.P0
        Hx += dot(self.G_lr * np.dot(self.S, dP), self.S)
        Hx += self.G_lr * dot(dot(self.S, dP), self.S)
        Hx += dot(dot(self.S, dP*self.G_lr), self.S)
        Hx += dot(self.S, dot(dP, self.S) * self.G_lr)
        Hx *= -1.0/8.0
        return Hx
    def lc_exchange_energy(self):
        """
        compute the long range exchange contribution to the total
        ground state energy
        """
        if self.lc_implementation == "old":
            # The exact exchange is added for the full density matrix
            # including the reference density
            dP = self.P
        else:
            # The exact exchange is added only for the difference density
            # as done in J. Chem. Phys. 143, 184107 (2015), where it is assumed
            # that the exchange energy due to the reference density P0
            # is already contained in the parametrization.
            dP = self.P-self.P0
        E_HF_x  = np.sum(np.dot(self.S, np.dot(dP, self.S)) * dP * self.G_lr)
        E_HF_x += np.sum(np.dot(self.S, dP) * np.dot(dP, self.S) * self.G_lr)
        E_HF_x *= -1.0/8.0
        return E_HF_x
    # for LAMBDA-diagnostic
    def _construct_gaussian_overlap(self):
        """
        the charge distribution on atom A is modeled by a spherical Gaussian

         F_A(r) = 1/(2*pi*sigma_A^2)^(3/2) * exp(-1/2 r^2/(2*sigma_A^2))
         
        where the width sigma_A is inversely proportional to the Hubbard parameters U_A

         sigma_A = 1.329/sqrt(8 * ln(2))

        This function computes Omega_AB = int F_A(r-RA) F_B(r-RB) d^3r
        """
        Natoms = len(self.atomlist)
        gaussOmega = zeros((Natoms, Natoms))
        # compute width sigma_i
        sigma = zeros(Natoms)
        for i,(Zi,posi) in enumerate(self.atomlist):
            Ui = self.hubbard_U[Zi]
            sigma[i] = 1.329/sqrt(8.0*log(2.0)) * 1.0/Ui

        for i,(Zi,posi) in enumerate(self.atomlist):
            Ri = array(posi)
            for j,(Zj,posj) in enumerate(self.atomlist):
                Rj = array(posj)
                Rij = Ri-Rj
                sigAB = sigma[i]**2 + sigma[j]**2
                gaussOmega[i,j] = 1.0/pow(2.0*pi*sigAB,3.0/2.0) \
                    * exp(-0.5 * 1.0/sigAB * dot(Rij,Rij))
#                gaussOmega[i,j] = exp(-0.5 * 1.0/sigAB * dot(Rij,Rij))
        self.gaussOmega = gaussOmega
        return gaussOmega
    def getGaussianOverlap(self):
        if hasattr(self, "gaussOmega"):
            return self.gaussOmega
        else:
            return self._construct_gaussian_overlap()

    #############################################
    # optimally tuned range separation parameter
    # see "Excitation Gaps of Finite-Sized Systems from Optimally Tuned Range-Separated Hybrid Functionals"
    # by Leeor Kronik et.al. J. Chem. Theory Comput. 2012, 8, 1515-1531
    #############################################
    def tune_long_range_radius(self):
        """
        minimize  J^2(R_lr) = (e_HOMO(R_lr) + IE(R_lr))^2
                            = (e_HOMO(R_lr) + E_gs(N-1,R_lr) - E_gs(N,R_lr))^2
        """
        if self.tune_range_radius == 0:
            return
        print "*******************************"
        print "* Tuning of long range radius *"
        print "*******************************"
        charge = self.charge
        self.dqN = None
        self.dqNm1 = None

        def Jfunc(g):
            self.long_range_radius = g
            ret = self.gm_lc.gamma_AOwise(self.atomlist, self.valorbs, self.distances, self.directions)
            self.gamma_lr, self.G_lr = ret[0], ret[2]
            # E_gs(N,g)
            if self.dqN != None:
                self.dq = self.dqN
            self.setCharge(charge)
            self.runSCC()
            self.getEnergies()
            EN = self.E_tot
            HOMO,LUMO = self.getFrontierOrbitals()
            eHOMO = self.orbe[HOMO]
            self.dqN = np.copy(self.dq)
            # E_gs(N-1,g)
            if self.dqNm1 != None:
                self.dq = self.dqNm1
            self.setCharge(charge-1)
            self.runSCC()
            self.getEnergies()
            ENm1 = self.E_tot
            self.dqNm1 = np.copy(self.dq)
            # reset charge
            self.setCharge(charge)
            IE = ENm1 - EN
            J2 = pow(eHOMO + IE,2)
            J = np.sqrt(J2)
            print "R_lr = %4.7f bohr" % g
            print "=========================="
            print "  IE     = %4.7f Hartree   %4.7f eV" % (IE, IE*hartree_to_eV)
            print "  -eHOMO = %4.7f Hartree   %4.7f eV" % (-eHOMO,-eHOMO*hartree_to_eV)
            print "  J      = %4.7f Hartree   %4.7f eV" % (J, J*hartree_to_eV)
            return IE, -eHOMO, J
        if self.save_tuning_curve != None:
            print "Calculate tuning curve"
            gs = np.linspace(0.5, 15.0, 100)
            lr = self.long_range_radius # save original value as it is modified 
                                        # in map
            IEhomoJ = map(Jfunc, gs)
            self.long_range_radius = lr # restore old value
            tuning_curve = np.c_[gs, np.array(IEhomoJ)]
            # tuning curve contains the following columns
            #
            # R_lr/bohr  IE / hartree  -Energy(HOMO) / hartree  J=|IE + eHOMO| / hartree"
            #
            np.savetxt(expandvars(expanduser(self.save_tuning_curve)), tuning_curve)

        res = fmin(lambda g: Jfunc(g)[-1], np.array([self.long_range_radius]))[0]
#        res = fminbound(J2func, 0.1, 15.0) 
        self.long_range_radius = res
        print "Optimal long range radius: %.7f bohr" % self.long_range_radius
        print "***********************"
        print "*   Finished Tuning   *"
        print "***********************"
    ########################################################################
    # Periodic DFTB (experimental)
    # This code is rather old and primitive. The DIIS doesn't seem to work
    # for charge-consistent periodic calculations.
    ########################################################################
    def gamma_function(self, Cij,Rij, deriv=0):
        """
        gamma function gamma_ij(R_ij) = erf(C_ji*R_ij)/R_ij describes the Coulomb
        energy of two spherically symmetric Gaussian charge distributions with
        widths FWHMi and FWHMj, where C_ij = sqrt(4*ln(2)/(FWHMi^2 + FWHMj^2))
        see Koskinen/Maekinen eq. (27)
        """
        assert Rij > 0.0  # otherwise we would need the limit erf(x)/x -- x-> 0 --> 2/sqrt(pi)
        if deriv == 0:
            g = erf(Cij*Rij)/Rij
            #
            """
            from matplotlib import pyplot as plt
            from numpy import linspace
            x = linspace(0.000000001, 5.0, 100)
            plt.plot(x, erf(Cij*x)/x)
            plt.show()
            """
            #
        elif deriv == 1:
            def erf_deriv(x):
                return 2.0/sqrt(pi) * exp(-x*x)
            g = (erf_deriv(Cij*Rij)*(Cij*Rij) - erf(Cij*Rij))/pow(Rij,2)
        else:
            raise Exception("Higher derivatives of gamma(Cij*Rij) not implemented")
        return g
    def _construct_gamma(self, deriv=0):
        """
        build the matrix gamma_ij (and its derivatives) as in eq.(31) of Koskinen/Maekinen
        """
        assert deriv in [0,1]
        Natoms = len(self.atomlist)
        gamma = zeros((Natoms,Natoms))
        for i,(Zi,posi) in enumerate(self.atomlist):
            for j,(Zj,posj) in enumerate(self.atomlist):
                Ui = self.hubbard_U[Zi]
                Uj = self.hubbard_U[Zj]
                if i == j:
                    if deriv == 0:
                        gamma[i,i] = Ui
                    else:
                        gamma[i,i] = 0.0
                elif i < j:
                    FWHMi = 1.329/Ui
                    FWHMj = 1.329/Uj
                    Cij = sqrt(4.0*log(2.0)/(pow(FWHMi,2)+pow(FWHMj,2)))
                    R = array(posi) - array(posj) 
                    Rij = norm(R)
                    gamma[i,j] = self.gamma_function(Cij,Rij,deriv=deriv)
                else:
                    gamma[i,j] = gamma[j,i]
        return gamma
    def _construct_periodic_gamma(self, lat, nmax=6, deriv=0):
        # TODO: Ewald summation
        T = lat.translation_vectors(nmax=nmax, space="real", include_origin=False)
        # T contains a list of lattice vectors excluding the vector t=(0,0,0)!
        assert deriv in [0,1]
        Natoms = len(self.atomlist)
        gamma = zeros((Natoms,Natoms))
        for i,(Zi,posi) in enumerate(self.atomlist):
            for j,(Zj,posj) in enumerate(self.atomlist):
                Ui = self.hubbard_U[Zi]
                Uj = self.hubbard_U[Zj]
                FWHMi = 1.329/Ui
                FWHMj = 1.329/Uj
                Cij = sqrt(4.0*log(2.0)/(pow(FWHMi,2)+pow(FWHMj,2)))
                R = array(posi) - array(posj) 
                for translation in T:
                    Rt = R - translation
                    Rij = norm(Rt)
                    gamma[i,j] += self.gamma_function(Cij,Rij,deriv=deriv)
        # add electrostatic interaction of charges within the unitcell
        gamma += self._construct_gamma(deriv=deriv)
        return gamma
    def _constructH0andS_translated(self, translation=zeros(3)):
        """
        construct the matrix elements of H0 and the overlaps S for a fixed
        nuclear geometry using Slater-Koster transformations. 
        This has to be performed only once, since during the DFTB iterations 
        only the mo coefficients and the partial charges change.

        Parameters:
        ===========
        translation: shift the position of the basis functions nu in the matrix element <mu|S or H|nu> by
           the vector tx,ty,tz, so calculate Smn(T) = <mu|nu(r-T)> and <mu|H|nu(r-T)> 
           This option is needed for periodic DFTB calculations
        """
        H0 = zeros((self.dim,self.dim))
        S = zeros((self.dim,self.dim))
        # iterate over atoms
        mu = 0
        zero_translation = False
        if la.norm(translation) < 1.0e-15:
            zero_translation = True
        for i,(Zi,posi) in enumerate(self.atomlist):
            # iterate over orbitals on center i
            for (ni,li,mi) in self.valorbs[Zi]:
                # iterate over atoms
                nu = 0
                for j,(Zj,posj0) in enumerate(self.atomlist):
                    posj = posj0 - translation
                    # iterate over orbitals on center j
                    for (nj,lj,mj) in self.valorbs[Zj]:
                        if mu == nu and zero_translation == True:
                            assert Zi == Zj
                            H0[mu,nu] = self.orbital_energies[Zi][(ni,li)] # use the true single particle orbitals energies
                            S[mu,nu] = 1.0 # orbitals are normalized to 1
                        else:
                            if Zi <= Zj:
                                # the first atom given to getHamiltonian() or getOverlap()
                                # has to be always the one with lower atomic number
                                if i == j and zero_translation == True:
                                    # different orbitals on the same atom should be orthogonal
                                    assert mu != nu
                                    s = 0.0
                                    h0 = 0.0
                                else:
                                    s  = self.SKT[(Zi,Zj)].getOverlap(li,mi,posi, lj,mj,posj)
                                    h0 = self.SKT[(Zi,Zj)].getHamiltonian0(li,mi,posi, lj,mj,posj)
                            else:
                                # swap atoms if Zj > Zi
                                s  = self.SKT[(Zj,Zi)].getOverlap(lj,mj,posj, li,mi,posi)
                                h0 = self.SKT[(Zj,Zi)].getHamiltonian0(lj,mj,posj, li,mi,posi)
                            H0[mu,nu] = h0
                            S[mu,nu] = s
                        nu += 1
                mu += 1
        return S, H0
    def nearest_neighbour_approximation(self, M, translation=zeros(3), cutoff=3.0):
        """
        Set those matrix elements in M to zero, that would not contribute
        in a nearest neighbour calculation
        """
        Mnn = copy(M)
        # iterate over atoms
        mu = 0
        for i,(Zi,posi) in enumerate(self.atomlist):
            # iterate over orbitals on center i
            for (ni,li,mi) in self.valorbs[Zi]:
                # iterate over atoms
                nu = 0
                for j,(Zj,posj0) in enumerate(self.atomlist):
                    posj = posj0 - translation
                    # iterate over orbitals on center j
                    for (nj,lj,mj) in self.valorbs[Zj]:
                        if norm(posj - posi) > cutoff:
                            Mnn[mu,nu] = 0.0
                        nu += 1
                mu += 1        
        return Mnn
    def displaced_matrix_elements(self, T, cutoff=500.0):        
        print "compute overlap matrices Smn(T) and Hmn(T) for all translation vectors ..."
        Smn = []
        H0mn = []

        for t in T:
            # TODO: parallelize over t
            if la.norm(t) > cutoff:
                SmnT = None
                H0mnT = None
            else:
                SmnT, H0mnT = self._constructH0andS_translated(translation=t)
            #
#            SmnT = self.nearest_neighbour_approximation(SmnT, translation=t, cutoff=4.0)
#            H0mnT = self.nearest_neighbour_approximation(H0mnT, translation=t, cutoff=4.0)
            #
            Smn.append(SmnT)
            H0mn.append(H0mnT)

        """
        J = Jobs()
        for t in T:
            J.job(lambda: self._constructH0andS_translated(translation=t))
        for SmnT, H0mnT in J.run_parallel():
            Smn.append(SmnT)
            H0mn.append(H0mnT)
        """
        return Smn, H0mn
    def solve_KS_at_kpoints(self, Smn, H0mn, h1, T, K):
        # solve the KS equations for each k-vector
        #  H0(k) C = S(k) C
        # where S(k) = sum_T e^(i k*T) S(T)
        # and H0(k) = sum_T e^(i k*T) H0(T)
        for ik,k in enumerate(K):
            # TODO: parallelize over k
            Sk = zeros((self.dim,self.dim), dtype=complex)
            H0k = zeros((self.dim,self.dim), dtype=complex)
            # build periodic hamiltonian and overlap
            for t,SmnT,H0mnT in itertools.izip(T,Smn,H0mn):
                if SmnT != None:
                    Sk += exp(1.0j*dot(k,t))*SmnT
                    H0k += exp(1.0j*dot(k,t))*H0mnT
            # Why are Sk and H0k not hermitian already??? -> because t is only positive
            Sk = 0.5 * (Sk + Sk.conjugate().transpose())
            H0k = 0.5 * (H0k + H0k.conjugate().transpose())
            Hk = H0k + h1 * Sk
            Hk = 0.5 * (Hk + Hk.conjugate().transpose())
            #
            orbe_k, orbs_k = eigh(Hk,Sk)
            # sort orbitals by energy
            orbe_k = orbe_k.real
            sort_indx = argsort(orbe_k)
            orbe_k = orbe_k[sort_indx]
            orbs_k = orbs_k[:,sort_indx]
            # construct density matrix for this k-point
            mu, fk = fermi_occupation(orbe_k, self.Nelec_val, T=0.0)
            occ_indx = where(fk > 0.0)[0]
            occ_orbs = orbs_k[:,occ_indx]
#            Pk = dot(fk[occ_indx]*occ_orbs,occ_orbs.conjugate().transpose())
            Pk = dot(fk[occ_indx]*occ_orbs.conjugate(),occ_orbs.transpose())
            # Mulliken charges
            qk = zeros(len(self.atomlist), dtype=complex)
            Qk = Pk * Sk
            q_mu = sum(0.5*(Qk + Qk.transpose().conjugate()), axis=1)
#            q_mu = sum(Qk,axis=1)
            # sum over all orbitals mu belonging to atom I
            for i in range(0, len(self.atomlist)):
                qk[i] = sum(q_mu[self.orbitals_on_atom[i]])
            yield (k, orbe_k, orbs_k, fk, qk, Sk, H0k)
    def runPeriodicSCC(self, lattice_vectors, ks, nmax=(3,3,3), maxiter=1000, scc_conv=1.0e-11, DIIS_memory=1, linear_mixing_coefficient=0.1):
        self.maxiter = maxiter
        lat = Periodic.CrystalSystem()
        lat.setLatticeVectors(lattice_vectors)
        T = lat.translation_vectors(nmax=nmax, space="real")
        # construct lattice to check if the symmetry is correct
        lattice = []
        atomlist = self.getGeometry()
        for t in T:
            for (Zi,posi) in atomlist:
                lattice.append( (Zi, posi + t) )
        #
        Smn, H0mn = self.displaced_matrix_elements(T)
        self.gamma = self._construct_periodic_gamma(lat, nmax=nmax)

        kpoints = lat.sample_BZ(nmax=nmax)

        diis = DIIS_80(DIIS_memory)
        # initial charge guess
        self.dq = zeros(len(self.atomlist))
        q_last = copy(self.q0)
        # repulsive energy
        self.E_rep = 0.0
        for t in T:
            self.E_rep += self.getRepulsiveEnergy(translation=t)
        self.E_rep /= float(len(T))
        self.E_nuc = self.E_rep
        #
        for self.i in range(0, self.maxiter):
            h1 = self._construct_h1(self.gamma, self.dq)
            # new Mulliken charges
            self.q = zeros(len(self.atomlist), dtype=complex) # Mulliken charges
            # band structure energy
            self.E_bs = 0.0 + 0.0j
            for (k, orbe_k, orbs_k, fk, qk, Sk, H0k) in self.solve_KS_at_kpoints(Smn, H0mn, h1, T, kpoints):
                self.q += qk
                Pk = self._density_matrix(orbs_k, fk, check=False)
                self.E_bs += sum(Pk * H0k)
            self.E_bs /= float(len(T))
            self.q /= float(len(kpoints))
            if DIIS_memory > 1:
                self.q = diis.next_approximation(self.q.real)
            else:
                alpha = linear_mixing_coefficient
                self.q = alpha*self.q + (1.0-alpha)*q_last
                pass
            # Mulliken excess charges
            for i in range(0, len(self.atomlist)):
                self.dq[i] = self.q[i] - self.q0[i] # excess charge
            # Coulomb energy
            self.E_coulomb = 0.5 * dot(self.dq, dot(self.gamma, self.dq))
            self.E_coulomb /= float(len(T))
            self.E_elec = self.E_bs + self.E_coulomb
            self.E_tot = self.E_elec + self.E_nuc
            if DIIS_memory > 1:
                self.relative_change = diis.relative_change()
            else:
                self.relative_change = norm(q_last - self.q)/norm(self.q)    
            #
            self.writeIteration(info=["convergence", "energies", "charges"])
            # 
            if self.relative_change < scc_conv:
                if self.verbose > 0:
                    print "!!! CONVERGED after %s iterations (relative change = %.2e < %.2e = threshold)!!!" \
                    % (self.i+1, self.relative_change, scc_conv)
                break

            q_last = copy(self.q)
        else:
            raise Exception("SCC Calculation did not converge after %s iterations!!!" % (self.i+1))
        # compute band structure on finer grid
        orbe = []
        print "Compute band structure on a finer grid with %s points" % len(ks)
        for (k, orbe_k, orbs_k, fk, qk, Sk, H0k) in self.solve_KS_at_kpoints(Smn, H0mn, h1, T, ks):
            orbe.append(orbe_k)
            if hasattr(self, "special_kpoint") == False:
                self.special_kpoint = zeros(3) # Gamma point
            if la.norm(k-self.special_kpoint) < 1.0e-10:
                self.orbs = orbs_k
                self.orbe = orbe_k
                self.f = fk
                # crystal orbital at a special point
                print "K-Vectors: %s" % k
                self.write_Crystal_orbital_coefficients(orbe_k, orbs_k, fk)
        return lattice, orbe, self.dq
    def setSpecialKPoint(self, special_kpoint):
        """
        define the k-point in the Brillouin zone for which crystal orbitals
        and KS orbital energies can be exported to Molden. 
        """
        self.special_kpoint = special_kpoint
    def write_Mulliken_charges(self):
        import string
        txt = "  Mulliken Charges (Excess Charges)\n"
        txt += "  =================================\n"
        for i,(Zi,posi) in enumerate(self.atomlist):
            txt += "%s: %.3f\t(%.3f)\n" % (string.center("%s-%s" % (atom_names[Zi-1], i),12), self.q[i].real, self.dq[i].real)
        print txt
        print "sum q0 = %s" % sum(self.q0)
        print "sum q  = %s" % sum(self.q)
        print "sum dq = %s" % sum(self.dq)        
    def write_Crystal_orbital_coefficients(self, orbe_k, orbs_k, fk):
        txt = ""
        txt += "  CO Coefficients (Real part)\n"
        txt += "  ===========================\n"
        txt += annotated_matrix(vstack((orbe_k.real, fk, orbs_k.real)),\
             ["en. (a.u.)", "occ.", "-"] + self.orbital_names,\
             ["CO%s" % i for i in range(0, len(orbe_k))], colwidth=10)
        txt += "  CO Coefficients (Imaginary part)\n"
        txt += "  ================================\n"
        txt += annotated_matrix(vstack((orbe_k.real, fk, orbs_k.imag)),\
             ["en. (a.u.)", "occ.", "-"] + self.orbital_names,\
             ["CO%s" % i for i in range(0, len(orbe_k))], colwidth=10)
        print txt
    def runPeriodicNonSCC(self, lattice_vectors, ks, dq=None, nmax=(10,10,10)):
        """
        Parameters:
        ===========
        lattice_vectors: list of 3d vectors, translations by integer multiples 
          of these vectors are symmetry operations
        ks: list of 3D k-vectors for which the band structure should be calculated
        nmax: number of translations to include along the lattice vectors, e.g. (40,40,40)
           for a 3D system, or (40,1,1) for a 1D system
        Returns:
        ========
        band_structure: list of numpy arrays with KS energies levels for each k-vector
        """
        assert len(lattice_vectors) <= 3
        # create all translations that are compatible with the symmetry
        T = []
        klm = list(itertools.product(range(-nmax[0], nmax[0]+1), range(-nmax[1],nmax[1]+1), range(-nmax[2],nmax[2]+1)))
        # (-k,-l,-m), ...(0,0,0), (0,0,1), ..., (k,l,m) for 3D grid
        for int_disp in klm:
            t = zeros(lattice_vectors[0].shape)
            for i in range(0, len(lattice_vectors)):
                t += int_disp[i] * lattice_vectors[i]   # t = k*a + l*b + m*c
            T.append(t)
        # construct lattice to check if the symmetry is correct
        lattice = []
        atomlist = self.getGeometry()
        for t in T:
            for (Zi,posi) in atomlist:
                lattice.append( (Zi, posi + t) )
        Smn, H0mn = self.displaced_matrix_elements(T)
        if dq != None:
            gamma = self.construct_periodic_gamma(lat)
            h1 = self._construct_h1(gamma, dq)
        else:
            h1 = zeros((self.dim,self.dim))

        # solve the KS equations for each k-vector
        #  H0(k) C = S(k) C
        # where S(k) = sum_T e^(i k*T) S(T)
        # and H0(k) = sum_T e^(i k*T) H0(T)
        orbe = []

        from numpy import identity

        self.q = zeros(len(self.atomlist), dtype=complex) # Mulliken charges
        self.dq = zeros(len(self.atomlist), dtype=complex) # Mulliken excess charges
        print "perform DFTB calculation for each k-vector ..."
        for k in ks:
            # TODO: parallelize over k
#            print "K-Vector = %s" % k
            Sk = zeros((self.dim,self.dim), dtype=complex)
            H0k = zeros((self.dim,self.dim), dtype=complex)
            # build periodic hamiltonian and overlap
            for t,SmnT,H0mnT in itertools.izip(T,Smn,H0mn):
                if SmnT != None:
                    Sk += exp(1.0j*dot(k,t))*SmnT
                    H0k += exp(1.0j*dot(k,t))*H0mnT
            Sk = 0.5 * (Sk + Sk.conjugate().transpose())
            H0k = 0.5 * (H0k + H0k.conjugate().transpose())
            Hk = H0k + h1*Sk
            Hk = 0.5 * (Hk + Hk.conjugate().transpose())
            orbe_k, orbs_k = eig(Hk,Sk)
            # sort orbitals by energy
            orbe_k = orbe_k.real
            sort_indx = argsort(orbe_k)
            orbe_k = orbe_k[sort_indx]
            orbs_k = orbs_k[:,sort_indx]
            orbe.append(orbe_k.real)
            # construct density matrix for this k-point
            mu, fk = fermi_occupation(orbe_k, self.Nelec_val, T=0.0)
            self.f = fk # UGLY, select GAMMA point
            occ_indx = where(fk > 0.0)[0]
            occ_orbs = orbs_k[:,occ_indx]
            Pk = dot(fk[occ_indx]*occ_orbs,occ_orbs.transpose())
            # Mulliken charges
            qk = zeros(len(self.atomlist), dtype=complex)
            dqk = zeros(len(self.atomlist), dtype=complex)
            Qk = Pk * Sk
            q_mu = sum(0.5*(Qk + Qk.transpose().conjugate()), axis=1)
            # sum over all orbitals mu belonging to atom I
            for i in range(0, len(self.atomlist)):
                qk[i] = sum(q_mu[self.orbitals_on_atom[i]])
                dqk[i] = qk[i] - self.q0[i] # excess charge
#            print "sum qk = %s" % sum(qk)
#            print "sum dqk= %s" % sum(dqk)
            #
            self.q += qk
            self.dq += dqk
            if la.norm(k) < 1.0e-10:
                print "K-Vector: %s (Gamma)" % k
                self.write_Crystal_orbital_coefficients(orbe_k, orbs_k, fk)
        Omega_BZ = la.norm(np.cross(lattice_vectors[0], lattice_vectors[1]))
        print "Omega_BZ = %s" % Omega_BZ
        self.dq *= Omega_BZ / float(len(ks))
        self.q  *= Omega_BZ / float(len(ks))
        import string
        txt = "  Mulliken Charges (Excess Charges)\n"
        txt += "  =================================\n"
        for i,(Zi,posi) in enumerate(self.atomlist):
            txt += "%s: %.3f\t(%.3f)\n" % (string.center("%s-%s" % (atom_names[Zi-1], i),12), self.q[i].real, self.dq[i].real)
        print txt
        print "sum q0 = %s" % sum(self.q0)
        print "sum q  = %s" % sum(self.q)
        print "sum dq = %s" % sum(self.dq)
        return lattice, orbe
        
###############################
# Charge Constrained DFTB
#
# Constraints allow to find the density of the lowest state that satiffies certain conditions
# on the distribution of charge over molecular fragments. 
# If the charge distribution is chosen cleverly the constraint solution 
# might even approximate an excited charge transfer (CT) state.
# By enumerating all possible combinations of constraints 
# (e.g. a charge sitting on any of the monomers in a conducting polymer such as polypyrrole ) 
# we can get a coarse-grained set of potential energy surfaces (one state per constraint). Non-adiabatic dynamics within these states should be able
# to simulate charge transport.
###############################
    def loadConstraints(self, constraints_script=None):
        """
        read constraints from a python file. In the main body of the 
        script the variable "constraints" has to be defined as a list
        of tuples with atomic indeces and excess charges. For example
          constraints.py:

            constraints = [([0,1,2],-1.0), ([10,11,12],+1.0)]

        This fixed the excess charge on the molecular fragment consisting of
        atoms 0,1,2 to -1 and the excess charge on the fragment consisting
        of atoms 10,11,12 to +1, so that a constrained DFTB calculation
        finds the lowest state satisfying this condition.
        Note that the first atom has index 0 and not 1!

        !!!So far only the charge on a single fragment can be fixed!!!

        Parameters:
        ===========
        Experimental (not tested properly).constraints_script: path to script with constraints
        """
        if constraints_script == None:
            return
        dic = {}
        try:
            execfile(expandvars(expanduser(constraints_script)), dic)
        except IOError as e:
            print "Cannot load constraints"
            print e
        constraints = dic.get("constraints", [])
        self._setConstraints(constraints)
    def _setConstraints(self, constraints):
        self.constraints = constraints
        assert len(self.constraints) == 1 # only one constraint so far
    def _constraining_potentials(self):
        """
        calculate the constraining potential matrices by which the constraint hamiltonian
        is modified to enforce the constraint for the current values of the Lagrange multipliers.

        Returns:
        ========
        Vconstr: vector of matrices [V^1_{mu,nu}, V^2_{mu,nu},...,V^N_{mu,nu}]
        """
        if len(self.constraints) == 0:
            return []
        Vconstr = [zeros((self.dim,self.dim)) for F in range(0, len(self.constraints))]
        # iterate over atoms
        mu = 0
        for i,(Zi,posi) in enumerate(self.atomlist):
            # iterate over orbitals on center i
            for (ni,li,mi) in self.valorbs[Zi]:
                # iterate over atoms
                nu = 0
                for j,(Zj,posj) in enumerate(self.atomlist):
                    # iterate over orbitals on center j
                    for (nj,lj,mj) in self.valorbs[Zj]:
                        # iterate over constraints and build constraining potential matrices
                        for F,(fragment_indecesF, NF) in enumerate(self.constraints):
#                            if i in fragment_indecesF:
#                                Vconstr[F][mu,nu] = self.S[mu,nu]
                            #
                            if i in fragment_indecesF: 
                                Vconstr[F][mu,nu] += 0.5*self.S[mu,nu]
                            if j in fragment_indecesF:
                                Vconstr[F][mu,nu] += 0.5*self.S[mu,nu]
                        #
                        nu += 1
                mu += 1
        return Vconstr
    def _constraints_values(self, Vconstr):
        """
        compute the value of the constraints C_F for the current coefficients
        (which in turn depend on the current values of the Lagrange multipliers)

        Parameters:
        ===========
        Vconstr: vector of constraining potentials as calculated by _constraining_potentials()

        Returns:
        ========
        Constraints: vector with constraint values, [C1(lambdas),C2(lambdas),...,CN(lambdas)]
        """
        Constr = zeros(len(self.constraints))
        for F,(fragment_indecesF, NF) in enumerate(self.constraints):
            # CF = sum_(mu,nu) P_(mu,nu) V^F_(mu,nu) - ( sum_(I in F) q_0^I  + NF )
            Constr[F] = sum(self.P * Vconstr[F]) # sum over mu and nu
            for i in range(0, len(self.atomlist)):
                if i in fragment_indecesF:
                    Constr[F] -= self.q0[i]
            Constr[F] -= (-NF) # since -NF is the number of additional electron
                               # NF=+1 means we have 1 electron less
        return Constr
    def runConstrainedSCC(self, maxiter=1000, scc_conv=1.0e-6, constr_conv=1.0e-12, temperature=0.0, DIIS_memory=4, search_bracket=(-3.0, 3.0), **opts): # convergence 10e-12
        """
        run a constrained self-consistent-charge calculation 
        where the Lagrange multiplier is determined in an inner loop
        by solving the non-SCC Kohn-Sham equations for fixed Mulliken charges.
        In an outer loop the Mulliken charges are determined self-consistently.

        Parameters:
        ===========
        constr_conv: convergence threshold for constraints
        temperature (in Kelvin): occupation of orbitals is smeared out 
          by Fermi distribution
        search_bracket: search this interval (a,b) for sign changes when adjusting
          the lagrange parameter so as to satisfy the constraint
        """
        self.maxiter = maxiter
        self.temperature = temperature
        # outer loop
        # constraining potential depends on overlap matrix which does not change
        Vconstr = self._constraining_potentials()
        # guess initial partial charges {dq_j}
        self.dq = zeros(len(self.atomlist))
        # initial guess for Lagrange multipliers
        lagrange = zeros(len(self.constraints))
        def solveKSnonSCC(lagrange):
            """find solution to KS equation for fixed values 
            of Lagrange multipliers"""
            self.Vc = zeros(self.S.shape)
            for F in range(0, len(self.constraints)):
                self.Vc += lagrange[F]*Vconstr[F]
            # build hamiltonian 
            h1 = self._construct_h1(self.gamma, 0.0*dq)
            # point charges
            if len(self.point_charges) > 0:
                hpc = self._construct_h_point_charges()
                h1 += hpc
            H = self.H0 + h1*self.S - self.Vc 
            # solve generalized eigenvalue problem
            #   H C = e S C
            # to get Kohn-Sham energies orbe[i] and 
            # corresponding orbital coefficients orbs[:,i]
            self.orbe, self.orbs = eigh(H,self.S)
            # update partial charges using Mulliken analysis
            self._constructDensityMatrix()
            self._MullikenAnalysis()
            # compute energies
            self.getEnergies()
        def find_sign_change(f, x1, x2):
            xs = [x1,x2]
            fs = []
            for x in xs:
                fs.append(f(x))
            for k in range(0, 10):
                for l in range(0, len(xs)-1):
                    if fs[l]*fs[l+1] < 0.0:
                        bracket = (xs[l], xs[l+1])
                        return bracket
                else:
                    # subdivide xs
                    xs_subdiv = []
                    fs_subdiv = []
                    for l in range(0, len(xs)-1):
                        xs_subdiv.append(xs[l])
                        fs_subdiv.append(fs[l])
                        # add midpoint
                        xnew = 0.5*(xs[l]+xs[l+1])
                        fnew = f(xnew)
                        xs_subdiv.append(xnew)
                        fs_subdiv.append(fnew)
                        xs_subdiv.append(xs[l+1])
                        fs_subdiv.append(xs[l+1])
                    xs = xs_subdiv
                    fs = fs_subdiv
            else:
                print "xs = %s" % xs
                print "fs = %s" % fs
                raise Exception("Cannot find any sign change by subdiving interval!")
        def bisect(f, x1, x2):
            #between x1 and x2 f(x) changes sign.
            f1 = f(x1)
            f2 = f(x2)
            if f1*f2 > 0.0:
                raise ValueError("f does not change sign between x1 = %s and x2 = %s, f(x1) = %s and f(x2) = %s" % (x1, x2, f1,f2))
                #(x1,x2) = find_sign_change(f, x1, x2)
            if f1 < 0.0:
                a = x1
                b = x2
            else:
                a = x2
                b = x1
            while True:
                midpoint = (a+b)/2.0
                fm = f(midpoint)
                if fm < 0.0:
                    a = midpoint
                else:
                    b = midpoint
                if abs(a-b) < 1.0e-20:
                    raise Exception("Interval in bisection has shrunk to a point")
                yield fm,midpoint,(a,b)
        def f(x):
            solveKSnonSCC([x])
            Constr = self._constraints_values(Vconstr)
#            print "Constr = %s" % Constr
            return Constr[0]

        self.relative_change = 0.0
        diis_dq = DIIS_80(DIIS_memory) # does not really work
        dq = zeros(len(self.atomlist))
        cur_constr = 1.0e10
        for j in range(0, maxiter):
            if self.verbose > 0:
                print "*** SCC iteration: %s ***" % j
            lagrange_arr = []
            constr_arr = []
            cur_constr = 1.0e10
            self.i = 0
            for fm,midpoint,(a,b) in bisect(f,search_bracket[0], search_bracket[1]):
                cur_constr = abs(fm)
                lagrange_arr.append( midpoint )
                constr_arr.append(fm)
                if self.verbose > 0:
                    print "*** adjusting Lagrange multiplier to satisfy constraint ***"
                    print "    deviation of constraint: %.7f (tolerance = %.7f)" % (abs(fm), constr_conv)
                    print "    bracket a = %s b = %s fm=%s cur_constr=%s" % (a,b, fm, cur_constr)

                    self.writeIteration()
                if cur_constr < constr_conv: #or abs(a-b) < 1.0e-14: # NOOOO
#                    diis_dq.reset()
                    break
                if self.i > maxiter:
                    from numpy import linspace, argmin
                    lagrange_arr = array(lagrange_arr)
                    constr_arr = array(constr_arr)
                    lagrange_arr2 = linspace(search_bracket[0], search_bracket[1], 1000)
                    constr_arr2 = array([f(lag) for lag in lagrange_arr2])
                    imin = argmin(abs(constr_arr2))
                    print "best lagrange multiplier = %s (constraint = %s)" % (lagrange_arr2[imin], constr_arr2[imin])
                    from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, savefig
                    xlabel("Lagrange multiplier", fontsize=15)
                    ylabel("deviation from constraint", fontsize=15)
                    plot(lagrange_arr, constr_arr, "o", lw=2, label="bisection points")
                    plot(lagrange_arr2, constr_arr2, lw=2)
                    legend()
                    savefig("/tmp/bisection_failed.png")
                    show()
                    raise Exception("Constraint could not be satisfied within %d iterations!" % maxiter)
                self.i += 1
            if DIIS_memory > 1 and j > DIIS_memory:
                    #and diis_dq.nr_trial_vectors() > DIIS_memory:
                dq = diis_dq.next_approximation(self.dq)
                rel_change = diis_dq.relative_change()
            else:
                #dq = 0.4*dq + 0.6*self.dq
                dq = 0.7*dq + 0.3*self.dq
                #dq = 0.9*dq + 0.1*self.dq
                rel_change = sum(abs(dq - self.dq))/sum(abs(self.dq))
                #
#                print "nr trial vectors = %s" % len(diis_dq.trial_vectors)
                if DIIS_memory > 1:
                    diis_dq.next_approximation(dq) # add trial vectors
                #
            if self.verbose > 0:
                print "relative change in Mulliken charges: %0.7f (tolerance = %.7f)" % (rel_change, scc_conv)
                print "constraint: %0.7f (tolerance = %.7f)" % (cur_constr, constr_conv)
            if rel_change < scc_conv and cur_constr < 1.0e-5:#constr_conv:
                print "CONVERGED    cur_constr = %.9f (constr_conv = %.9f)" % (cur_constr, constr_conv)
                self.lagrange = [midpoint]
                if self.verbose > 0:
                    print "LAGRANGE MULTIPLIER = %s" % self.lagrange
                break
        else:
            raise Exception("Not converged in %d iterations!" % maxiter)
        """
        # find root the hard way
        from numpy import linspace
#        lagrange_arr = linspace(-0.802, -0.8015, 100)
#        lagrange_arr = linspace(-0.7210, -0.7205, 1000)
#        lagrange_arr = linspace(0.3915, 0.392, 100)
        lagrange_arr = linspace(-10.0, 10.0, 200)
        constr_arr = []
        entot_arr = []
        for lagrange in lagrange_arr:
            print "LAGRANGE = %s" % lagrange
            solveKSnonSCC([lagrange])
            self.getEnergies()
            # evaluate constraints
            Constr = self._constraints_values(Vconstr)
            constr_arr.append(Constr[0])
            entot_arr.append(self.E_tot)
        constr_arr = array(constr_arr)
        from matplotlib.pyplot import plot, show, legend, xlabel, ylabel
        xlabel("Lagrange multiplier", fontsize=15)
        ylabel("Constraint $\int w_F \rho - N_F$")
        plot(lagrange_arr, constr_arr, lw=2)
#        plot(lagrange_arr, entot_arr, label="total energy")
        legend()
        show()
        """

###############################
    def writeIteration(self, info=["convergence", "energies", "orbitals", "charges"]):
        """
        print out information about SCC iteration: 
          convergence
          total energies
          Kohn-Sham orbital energies
          partial charges
          MO coefficients
        """
        import string
        txt  = "      *******************\n"
        txt += "      * Iteration: %.2d *\n" % self.i
        txt += "      *******************\n"
        if "convergence" in info:
            txt += "  Convergence\n"
            txt += "  ===========\n"
            txt += "relative change: %.6e\n" % self.relative_change
        if "energies" in info:
            txt += "  Electronic Energies\n"
            txt += "  ===================\n"
            txt += "band structure energy E_bs: %s %s\n" % (string.rjust("%.7f hartree" % self.E_bs,20), string.rjust("%.7f eV" % (self.E_bs*hartree_to_eV),20))
            txt += "Coulomb energy E_coulomb  : %s %s\n" % (string.rjust("%.7f hartree" % self.E_coulomb,20), string.rjust("%.7f eV" % (self.E_coulomb*hartree_to_eV),20))
            if hasattr(self, "E_multipoles"):
                txt += "            E_multipoles  : %s %s\n" % (string.rjust("%.7f hartree" % self.E_multipoles,20), string.rjust("%.7f eV" % (self.E_multipoles*hartree_to_eV),20))

            if hasattr(self, "Vcoul_pc"):
                txt += "point charges E_pc        : %s %s\n" % (string.rjust("%.7f hartree" % self.E_point_charges,20), \
                                                                    string.rjust("%.7f eV" % (self.E_point_charges*hartree_to_eV),20))
            if hasattr(self, "E_onsite"):
                txt += "Coulomb on-site E_onsite  : %s %s\n" % (string.rjust("%.7f hartree" % self.E_onsite,20), \
                                                                    string.rjust("%.7f eV" % (self.E_onsite*hartree_to_eV),20))
            if hasattr(self, "E_HF_x"):
                txt += "long range HF-x E_x       : %s %s\n" % (string.rjust("%.7f hartree" % self.E_HF_x,20), \
                                                                    string.rjust("%.7f eV" % (self.E_HF_x*hartree_to_eV),20))
                
            txt += "total electronic energy   : %s %s\n" % (string.rjust("%.7f hartree" % self.E_elec,20), string.rjust("%.7f eV" % (self.E_elec*hartree_to_eV),20))
            txt += "repulsion energy          : %s %s\n" % (string.rjust("%.7f hartree" % self.E_rep,20), string.rjust("%.7f eV" % (self.E_rep*hartree_to_eV),20))
            if hasattr(self, "dispersion"):
                txt += "dispersion energy         : %s %s\n" % (string.rjust("%.7f hartree" % self.E_disp,20), string.rjust("%.7f eV" % (self.E_disp*hartree_to_eV),20))
            if self.qmmm != None:
                txt += "QM/MM energy              : %s %s\n" % (string.rjust("%.7f hartree" % self.E_qmmm,20), string.rjust("%.7f eV" % (self.E_qmmm*hartree_to_eV),20))
            if self.cavity != None:
                txt += "Cavity energy             : %s %s\n" % (string.rjust("%.7f hartree" % self.E_cav,20), string.rjust("%.7f eV" % (self.E_cav*hartree_to_eV),20))
            if hasattr(self, "solvent_cavity") and self.solvent_cavity.implicit_solvent != 0:
                txt += "solvent screening energy  : %s %s\n" % (string.rjust("%.7f hartree" % self.E_screening,20), string.rjust("%.7f eV" % (self.E_screening*hartree_to_eV),20))
            txt += "total energy              : %s %s\n" % (string.rjust("%.7f hartree" % self.E_tot,20), string.rjust("%.7f eV" % (self.E_tot*hartree_to_eV),20))
        if "orbitals" in info:
            txt += "HOMO-LUMO gap             : %s %s\n" % (string.rjust("%.7f hartree" % self.HLgap,20), string.rjust("%.7f eV" % (self.HLgap*hartree_to_eV),20))
            txt += "  Orbital Energies\n"
            txt += "  ================\n"
            sort_indx = argsort(self.orbe) # calculated redundantly in fermi_occupation !!
            for i,a in enumerate(sort_indx):
                txt += "%s: %s %s  (%.3f e)\n" % (string.rjust(str(i+1),6), \
                                              string.rjust("%.7f hartree" % self.orbe[a], 20), \
                                              string.rjust("%.7f eV" % (self.orbe[a]*hartree_to_eV), 20), \
                                              self.f[a])
            if self.verbose > 1:
                txt += "  MO Coefficients\n"
                txt += "  ===============\n"
                txt += annotated_matrix(vstack((self.orbe[sort_indx], self.f[sort_indx], self.orbs[:,sort_indx])),\
                                     ["en. (a.u.)", "occ.", "-"] + self.orbital_names,\
                                     ["MO%s" % i for i in range(0, len(self.orbe))], colwidth=10)
        if "charges" in info:
            txt += "  Mulliken Charges   (Partial Charges)"
            if self.mulliken_dipoles == 1:
                txt += "[Partial Mulliken Dipoles]\n"
            else:
                txt += "\n"
            txt += "  NOTE: charges are in units of e-!\n"
            txt += "  =========================================================================\n"
            for i,(Zi,posi) in enumerate(self.atomlist):
                txt += "%s: %.3f\t(%+.3f)" \
                    % (string.center("%s-%s" % (atom_names[Zi-1], i),12), \
                       self.q[i], self.dq[i])
                if self.mulliken_dipoles == 1:
                    txt += "\t[%+.3f %+.3f %+.3f]\n" \
                        %(self.ddip[i,0], self.ddip[i,1], self.ddip[i,2])
                else:
                    txt += "\n"
            txt += "  =========================================================================\n"
            txt += "  sum of charges        : %+.3f\n" % np.sum(self.q)
            txt += "  sum of partial charges: %+.3f\n" % np.sum(self.dq)
            if self.mulliken_dipoles == 1:
                txt += "  Dipoles of electronic density\n"
                txt += "  =============================\n"
                txt += "  sum of dipoles        : [%+.3f %+.3f %+.3f]\n" % tuple(np.sum(self.dip, axis=0))
                txt += "  sum of partial dipoles: [%+.3f %+.3f %+.3f]\n" % tuple(np.sum(self.ddip, axis=0))
                txt += "  point charge approx.  : [%+.3f %+.3f %+.3f]\n" % tuple( self.getMullikenDipoleMoment() )
        print txt
    def getEnergies(self):
        """
        compute electronic energies

        Returns:
        ========
        band structure energy E_bs, Coulomb energy E_coulomb, total electronic energy, repulsive energy and total energy
        """
        self.E_bs = sum(self.P*self.H0) # band structure energy
                                 # E_bs = sum_a f_a <a|H0|a>
        # Coulomb energy from monopoles
        self.E_coulomb = 0.5*np.dot(self.dq, np.dot(self.gamma, self.dq))
        if self.mulliken_dipoles == 1:
            # Coulomb energy from monopoles + dipoles
            self.E_multipoles = 0.5*np.dot(self.dqm, np.dot(self.gamma_multipoles, self.dqm))
        self.E_elec = self.E_bs + self.E_coulomb
        if self.long_range_correction == 1:
            # compute long range Hartree-Fock exchange
            self.E_HF_x = self.lc_exchange_energy()
            self.E_elec += self.E_HF_x
        if self.onsite_correction == 1:
            self.E_onsite = self._onsite_energy()
            self.E_elec += self.E_onsite
        if hasattr(self, "solvent_cavity") and self.solvent_cavity.implicit_solvent == 1:
            # solvent energy is included in electrostatic energy
            self.E_screening = self.solvent_cavity.getScreeningEnergy(self.dq)
        if hasattr(self, "Vcoul_pc"):
            # point charges
            self.E_point_charges = dot(self.Vcoul_pc, self.dq)
            self.E_elec += self.E_point_charges
        if self.i == 0:
            # the repulsive potential, the dispersion correction and
            # the QM/MM energy only depend on the nuclear geometry and do not change during
            # the SCF cycle, therefore they are calculated only once at
            # the start
            self.E_rep = self.getRepulsiveEnergy()
            self.E_disp = self.getDispersionEnergy()
            self.E_qmmm = self.getQMMM_Energy()
            self.E_cav = self.getCavity_Energy()
            self.E_nuc = self.E_rep + self.E_disp + self.E_qmmm + self.E_cav
        self.E_tot = self.E_elec + self.E_nuc
        return self.E_bs, self.E_coulomb, self.E_elec, self.E_nuc, self.E_tot
    @T.timer
    def getEnergy(self, **opts):
        """
        perform a charge consistent calculation and return total energy.
        This function together with .getGradient() is the interface exhibited
        by DFTB2 to the dynamics module.
        """        
        # use partial charges from previous dynamics step as initial guess for next step
        # This leads to discontinuous energy (Why?)
#        self.runSCC(charge_guess=getattr(self, "dq", None), **opts) # WARNING: dq might not sum exactly to zero
        #
        self._distance_matrix()
        self._proximity_matrix()
        self._center_of_nuclear_charge()
        # overlap between charge fluctuation functions
        self._construct_gaussian_overlap()
        # construct H0 and S matrices
        self.S, self.H0 = self._constructH0andS()
        # construct dipole matrix elements
        try:
            self.D = self._constructDipoleMatrix()
        except NoSKDipolesException as e:
            # with other parametrizations dipole matrix elements are not available
            print e
            Norb,Norb = self.S.shape
            self.D = np.zeros((Norb,Norb,3))
        #
        self.gamma = self.gm.gamma_atomwise(self.atomlist, self.distances)[0]
        ##### solvent cavity ##############
        if hasattr(self, "solvent_cavity") and self.solvent_cavity.implicit_solvent == 1:
            # The screening effect due to a solvent cavity can be
            # incorporated by a modified gamma-matrix
            gamma_solvent = self.solvent_cavity.constructCOSMO()
            # The modified gamma-matrix accounts for a reduced
            # interaction between the Mulliken charges.
            self.gamma += gamma_solvent
            
        # multipoles
        if self.mulliken_dipoles == 1:
            self.gamma_multipoles = self.gm_multipoles.gamma_atomwise(self.atomlist, self.distances, self.directions)[0]
        ##### long-range correction #############
        if self.long_range_correction == 1:
            if self.tune_range_radius == 1:
                self.tune_long_range_radius()
            else:
                ret = self.gm_lc.gamma_AOwise(self.atomlist, self.valorbs, self.distances, self.directions)
                self.gamma_lr, self.G_lr = ret[0], ret[2]
        ###### on-site correction ##############
        if self.onsite_correction == 1:
            self._load_onsite_integrals()
        ########################################

        if len(self.constraints) == 0:
#            print "SCC"
            self.runSCC(**opts)
#            self.runNonSCC()  
        else:
            # constrained SCC calculation
            print "Constrained SCC calculation"
            self.runConstrainedSCC(**opts)
        Etot = self.getEnergies()[-1]
        return Etot
    def getRepulsiveEnergy(self, translation=None):
        """
        compute energy due to core electrons and nuclear repulsion
        """
        # nuclear energy because of repulsive potential
        E_nuc = 0.0
        for i,(Zi,posi) in enumerate(self.atomlist[1:]):
            for (Zj,posj) in self.atomlist[:i+1]:
                if Zi > Zj:
                    Z1 = Zj
                    Z2 = Zi
                else:
                    Z1 = Zi
                    Z2 = Zj
                R = array(posi) - array(posj)
                if translation != None:
                    R -= translation
                Rij = norm(R)
                # nucleus-nucleus and core electron repulsion
                E_nuc += self.VREP[(Z1,Z2)].getVrep(Rij)
        return E_nuc
    def getDispersionEnergy(self):
        """
        compute dispersion correction for energy
        """
        if hasattr(self, "dispersion"):
            return self.dispersion.getEnergy(self.atomlist)
        else:
            return 0.0
    def getQMMM_Energy(self):
        """
        compute
         E^QMMM = E^MM(I+O) - E^MM(I) 
                 + sum_(i in O) sum_(i<j in O) (dQ^MM_i *dQ^MM_j)/|Ri-Rj|
                 + sum_(i in I) sum(j in O)    (dQ^MM_i *dQ^MM_j)/|Ri-Rj|
        """
        if self.qmmm != None:
            return self.qmmm.getEnergy()
        else:
            return 0.0
    def getCavity_Energy(self):
        if self.cavity != None:
            if self.qmmm != None:
                atomlist = self.qmmm.getGeometryFull()
            else:
                atomlist = self.atomlist
            return self.cavity.getEnergy(atomlist)
        else:
            return 0.0
    def _energy_func(self, x, atomlist, which_en):
        """helper function for numerical differentiation"""
        atpos = XYZ.vector2atomlist(x, atomlist)
        self.setGeometry(atpos)
        self.getEnergy()
        # UGLY: getEnergies() is called twice, first inside getEnergy() and then here again
        if which_en == "Ebs":
            en = self.getEnergies()[0] # band structure gradient
        elif which_en == "Ecoulomb":
            en = self.getEnergies()[1]
        elif which_en == "Eelec":
            en = self.getEnergies()[2]
        elif which_en == "Enuc":
            en = self.getEnergies()[3]
        elif which_en == "Etot": 
            en = self.getEnergies()[4]
        return en
    # LR_TDDFTB accesses ground state quantities through the following functions
    # These functions should only be called after a converged ground state calculation
    def getOverlapMatrix(self):
        """overlap matrix between atomic orbitals"""
        if not hasattr(self, "S"):
            self._distance_matrix()
            self._proximity_matrix()
            # construct H0 and S matrices
            self.S, self.H0 = self._constructH0andS()
        return self.S
    def getKSCoefficients(self):
        """molecular orbital coefficients of Kohn-Sham orbitals"""
        return self.orbs
    def getKSEnergies(self):
        """Kohn-Sham energies of molecular orbitals"""
        return self.orbe
    def getOccupation(self):
        """occupation of KS orbitals"""
        return self.f
    def getFrontierOrbitals(self):
        """find indeces of HOMO and LUMO orbitals (starting from 0)"""
        HOMO = where(self.f > 0.0)[0][-1]
        LUMO = where(self.f == 0.0)[0][0]
        return HOMO,LUMO
    def getHOMO_LUMO_gap(self):
        """compute HOMO-LUMO gap in Hartree"""
        HOMO,LUMO = self.getFrontierOrbitals()
        gap = self.orbe[LUMO] - self.orbe[HOMO]
        return gap
    def getGeometry(self):
        return self.atomlist
    def getValorbs(self):
        """list of valence orbitals (nA,lA,mA) for each atom type ZA"""
        return self.valorbs
    def getGamma(self):
        """
        two electron integral matrix gamma_IJ
        """
        return self.gamma
    def getDensityMatrix(self):
        return self.P
    def getDipoleMatrix(self):
        try:
            DipMat = self._constructDipoleMatrix()
        except NoSKDipolesException as e:
            print e
            DipMat = None        
        return DipMat
    def getPartialCharges(self):
        return self.dq
    # 
    def getMullikenDipoleMoment(self):
        """
        calculate dipole moment from Mulliken charges, assuming each atom carries its
        Mulliken charge:
         D = sum_i Ri qi
        
        Returns:
        ========
        dipole moment as array([Dx,Dy,Dz])
        """
        dipole = zeros(3)
        for i,(Zi,posi) in enumerate(self.atomlist):
            dipole += (array(posi)-self.Rcc) * self.q[i]
        return dipole
    def getDipoleMoment(self):
        """
        expectation value of dipole operator in the ground state.
        requires a converged SCC calculation

        Returns:
        ========
        dipole moment as array([Dx,Dy,Dz])
        """
        DipMat = self.getDipoleMatrix()
        # Di = sum_(mu,nu) P_(mu,nu) Di_(mu,nu)   i=x,y,z
        print "density matrix"
        print self.P
        dipole = array([sum(self.P*DipMat[:,:,0]), \
                        sum(self.P*DipMat[:,:,1]), \
                        sum(self.P*DipMat[:,:,2])])
        return dipole

    def setSolventCavity(self, solvent_cavity):
        self.solvent_cavity = solvent_cavity

    def saveLocalizedOrbitals(self, localize_orbitals=None):
        """
        Localized Orbitals.localize_orbitals: Molecular orbitals are localized to disconnected fragments and are saved in the Molden format to the file 'localized_orbitals.molden'. The only localization method is Pipek-Mezey ('PM')
        """        
        if localize_orbitals == "PM":
            print "Pipek-Mezey localization"
            orbs_loc, orbe_loc, Uloc, frags = OrbitalLocalization.localize_pipek_mezey(self.atomlist, self.orbs, self.orbe, self.f, self.S, self.valorbs)
        else:
            if not (localize_orbitals is None):
                raise ValueError("Acceptable values for option 'localized_orbitals' are: 'PM' (for Pipek-Mezey localization), but got '%s'" % localize_orbitals)
            return
        
        # save orbitals in molden format
        nocc = len(self.f[self.f > 0.0])
        nvirt = len(self.f[self.f == 0.0])
        
        molden_file = "localized_orbitals.molden"

        molden = MoldenExporter(self)

        molden.setOrbitals(orbs_loc, orbe_loc, self.f)
        molden.export(molden_file=molden_file, molden_nr_occ=nocc, molden_nr_virt=nvirt)
            
        print "localized orbitals saved to '%s'" % molden_file
            
def fermi_occupation(orbe, Nelec_paired, Nelec_unpaired=0, T=0.0):
    """
    Find the occupation of single-particle state a at finite temperature T
    according to the Fermi distribution:
      f_a = f(en_a) = 2 /(exp(en_a - mu)/(kB*T) + 1)
    The chemical potential is determined from the condition that
      sum_a f_a = Nelec

    Parameters:
    ===========
    orbe: 1d numpy array, orbital energies
    Nelec_paired: number of paired electrons, these electron will be placed in the same orbital
    Nelec_unpaired: number of unpaired electrons, these electrons will sit in singly occupied orbitals (only works at T=0)
    T: temperature in Kelvin

    Returns:
    ========
    mu: chemical potential
    f: 1d numpy array, list of occupations f[a] for orbital a
       (in the same order as the energies in orbe)
    """
    from math import ceil, floor
    if T == 0.0:
        return fermi_occupation_T0(orbe, Nelec_paired, Nelec_unpaired)
    Nelec = Nelec_paired + Nelec_unpaired
    def fermi(en,mu):
        return 2.0/(exp((en-mu)/(kBoltzmann*T))+1.0)
    def func(mu):
        """find the root of this function to enforce sum_a f_a = Nelec"""
        sumfa = 0.0
        for en_a in orbe:
            sumfa += fermi(en_a, mu)
        return sumfa - Nelec
    sort_indx = argsort(orbe)
#    print "Nelec = %s" % Nelec
#    print "index = %s" % max((int(floor(Nelec/2))-1),0)
    HOMO = orbe[sort_indx[max(int(floor(Nelec/2))-1,0)]] # highest doubly occupied orbital
    LUMOp1 = orbe[sort_indx[min(int(ceil(Nelec/2))+1, len(sort_indx)-1)]] # LUMO+1
    # look for fermi energy in the interval [HOMO, LUMO+1] 
    # find root of (sum_a fa - Nelec) by bisection
    a = HOMO
    b = LUMOp1
#    print "Search for Fermi energy in interval [%s, %s]" % (a,b)
#    print "Orbital energies: %s" % orbe
    fa,fb = func(a), func(b)
    assert fa*fb <= 0.0 # sign change within interval
    dx = b-a
    dy = max(fa,fb)
    tolerance = 1.0e-8
    while (dx**2+dy**2 > tolerance):
        c = (a+b)/2.0
#        print "c=%s, dx = %s, dy = %s" % (c,dx,dy)
        fa, fc = func(a), func(c)
        if fa*fc <= 0.0:
            b,dx,dy = c,c-a,fc
        else:
            a,dx,dy = c,b-c,fc
    dN = func(c)
    mu = c
    # plot func(mu)
    """
    from numpy import linspace
    from matplotlib.pyplot import plot, show
    mu_arr = linspace(-2.0,2.0, 1000)
    plot(mu_arr, func(mu_arr))
    show()

    print "Nelec = %s" % Nelec
    print "orbital energies = %s" % orbe
    print "|N - sum_a f_a| = %s" % dN
    print "Chemical potential mu = %s" % mu
    """
    #
    assert abs(dN) < 1.0e-4
    f = fermi(orbe, mu)
#    print "f = %s" % f
    return mu,f

def fermi_occupation_T0(orbe, Nelec_paired, Nelec_unpaired=0):
    """
    Find the occupation of single-particle states at T=0
    """
    sort_indx = argsort(orbe)
    f = zeros(len(orbe))

    for a in sort_indx:
        if Nelec_paired > 0.0:
            f[a] = min(2.0,Nelec_paired)
            Nelec_paired -= 2.0
        else:
            if Nelec_unpaired > 0.0:
                f[a] = min(1.0,Nelec_unpaired)
                Nelec_unpaired -= 1.0

    return 0.0, f

def annotated_hessian(atomlist, hessian, **opts):
    """
    format the hessian so that columns and rows are
    labels according to the atomic coordinates with
    respect to which the derivatives are taken.

    Parameters:
    ===========
    atomlist: list of atoms (Zi,(xi,yi,zi))
    hessian: 2D numpy array
    the same additional options can be passed as to utils.annotated_matrix(...).

    Returns:
    ========
    txt: formatted hessian
    """
    labels = []
    for i,(Zi,posi) in enumerate(atomlist):
        for c in ["x","y","z"]:
            labels.append("%s-%d %s" % (atom_names[Zi-1],i,c))
    return utils.annotated_matrix(hessian, labels, labels, **opts)


if __name__ == "__main__":
    from DFTB.XYZ import read_xyz, extract_keywords_xyz
    from DFTB.Analyse.Cube import CubeExporter
    from DFTB.ImplicitSolvent import SolventCavity
    from DFTB.optparse import OptionParserFuncWrapper
    
    import sys
    import os.path
    import string
    
    usage = "Usage: %s <xyz-file>\n" % sys.argv[0]
    usage += "   options can be passed on the command line or\n" 
    usage += "   in a configuration file called 'dftbaby.cfg'\n"
    usage += "   --help option will provide more information\n"
    
    parser = OptionParserFuncWrapper([DFTB2.__init__, DFTB2.runSCC, DFTB2.saveLocalizedOrbitals, MoldenExporter.export, CubeExporter.exportCubes, SolventCavity.__init__], usage)
    
    (options, args) = parser.parse_args()
    if len(args) < 1:
        print usage
        exit(-1)
        
    xyz_file = args[0]
    kwds = extract_keywords_xyz(xyz_file)
    # read first structure to get the atom types
    structures = read_xyz(xyz_file)
    atomlist = structures[0]

    dftb2 = DFTB2(atomlist, **options)
    (scf_options, args) = parser.parse_args(dftb2.runSCC)
    
    # solvent cavity
    (solvent_options, args) = parser.parse_args(SolventCavity.__init__)
    solvent_cavity = SolventCavity(**solvent_options)
    dftb2.setSolventCavity(solvent_cavity)
    
    for i,atomlist in enumerate(structures):
        dftb2.setGeometry(atomlist, kwds.get("charge", 0.0))
        dftb2.getEnergy(**scf_options)
        print "TOTAL ENERGY OF STRUCTURE %d              : %s %s\n" % (i, string.rjust("%.7f hartree" % dftb2.E_tot,20), string.rjust("%.7f eV" % (dftb2.E_tot*hartree_to_eV),20))
#        print "Dipole Moment (exact): %s" % dftb2.getDipoleMoment()
#        print "Dipole Moment (Mulliken approximation): %s" % dftb2.getMullikenDipoleMoment()

    # localize orbitals
    (options,args) = parser.parse_args(dftb2.saveLocalizedOrbitals)
    dftb2.saveLocalizedOrbitals(**options)

    # export molden file with orbitals
    (options, args) = parser.parse_args(MoldenExporter(None).export)
    molden = MoldenExporter(dftb2)
    molden.export(**options)
#    print fermi_occupation(array([-0.2, -10.0, -4.0, -1.0, -0.01]), 5, 30.0)
#    print fermi_occupation_T0(array([-0.2, -10.0, -4.0, -1.0, -0.01]), 5)
    
    # write density to grid
    cube = CubeExporter(dftb2, name=os.path.basename(xyz_file).replace(".xyz", ""))
    (options, args) = parser.parse_args(cube.exportCubes)
    cube.exportCubes(**options)

    # write partial dipoles
    if dftb2.mulliken_dipoles == 1:
        Mulliken.save_partial_dipoles("/tmp/mulliken_dipoles.dat", dftb2.atomlist, dftb2.ddip)
    # writing timings
    print T

    print "FINISHED"
    
