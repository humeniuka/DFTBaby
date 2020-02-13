#!/usr/bin/env python
"""

"""
import numpy
from numpy import array, dot, linspace, argsort, where, mean
from numpy.linalg import norm
from scipy import interpolate
from scipy import integrate
import os.path

from DFTB import XYZ
from DFTB.AtomicData import hartree_to_eV, bohr_to_angs, atomic_number, atom_names
from DFTB import utils
from pairwise_decomposition import pairwise_decomposition_analytical as pairwise_decomposition
#from pairwise_decomposition import pairwise_decomposition_3at_analytical as pairwise_decomposition
#from pairwise_decomposition import pairwise_decomposition_numerical as pairwise_decomposition

import matplotlib as mpl
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size

class ReppotFitter():
    """
    class for fitting the repulsive potential

    
    """
    def __init__(self, Z1, Z2, name=""):
        """atomic numbers of atom pair"""
        self.Z1 = Z1
        self.Z2 = Z2
        """name of the atom pair, e.g. h-c"""
        self.name_atom_pair = "%s-%s" % (atom_names[Z1-1], atom_names[Z2-1])
        """list of repulsive forces for different internuclear distances
         {Ri,V'_rep(Ri)}
        """  
        self.nuclear_separation = []
        self.repulsive_gradient = []
        self.weights = []
        # remember from which curve each data points originates
        self.curve_ID = [] 
        self.curve_counter = 0
        self.curve_names = [] # for labeling plots
    def add_cutoff_curve(self, bond_lengths, weight):
        """
        add a curve with 0.0s beyond Rcut to force the repulsive potential
        and its derivatives to vanish beyond the cutoff radius.
        """
        for bl in bond_lengths:
            self.nuclear_separation.append(bl)
            self.repulsive_gradient.append(0.0)
            self.weights.append(weight)
            self.curve_ID.append( self.curve_counter )
        self.curve_names.append( "cutoff" )
        self.curve_counter += 1
    def add_polyatomic_curve(self, geometries, forces_wo_frep, forces_with_frep, \
                                 curve_name="polyatomic", weight=1.0, max_error=0.01, active_atoms=None):
        """
        compute generalized forces that act along the bond lengths and
        calculate repulsive forces for arbitrary non-linear molecule with more
        than 2 atoms.

        Parameters:
        ===========
        geometries: list of geometries (in atomlist format)
        forces_wo_frep: list of forces (in atomlist format) containing the DFTB forces only
           from the electronic hamiltonian
        forces_with_frep: list of forces (in atomlist format) containing the forces from a 
           higher level method (including nuclear repulsion)
        curve_name: some identifier for plotting
        weight: assign a weight of this data set during fitting
        max_error: for large bond lengths and repulsive potentials that cannot be represented
           as a sum of pairwise potentials, the error incurred from a pairwise decomposition
           can be large. Those geometries where the error per atom exceeds max_error 
           are not included in the fit. 
        active_atoms: list of atomic indeces. Only atoms in this list take part in the fit.
           If set to [], all atoms are included.
        """
        for atomlist, f_wo_rep, f_with_rep in zip(geometries, forces_wo_frep, forces_with_frep):
            # repulsive potential is the difference between true forces and electronic DFTB forces
            # in cartesian coordinates
            frep = XYZ.atomlist2vector(f_with_rep) - XYZ.atomlist2vector(f_wo_rep)
            forcelist = XYZ.vector2atomlist(frep,atomlist)
            from numpy.linalg.linalg import LinAlgError
            if active_atoms == []:
                # by default all atoms are included into the fit
                active_atoms = [i for i in range(0, len(atomlist))]
            try:
                bond_lengths, pair_forces, error = pairwise_decomposition(atomlist, forcelist)
            except LinAlgError:
                print "LINALG ERROR => do not include this geometry"
                continue
            if error > max_error:
                print "Error (%s) from pairwise decomposition too large (threshold: %s) => do not include this geometry!" % (error, max_error)
                continue

            for (i,j) in bond_lengths.keys():
                if (i not in active_atoms) or (j not in active_atoms):
                    continue
                Zi = atomlist[i][0]
                Zj = atomlist[j][0]
                if (Zi == self.Z1 and Zj == self.Z2) or (Zi == self.Z2 and Zj == self.Z1):
                    self.nuclear_separation.append( bond_lengths[(i,j)] )
                    self.repulsive_gradient.append( - pair_forces[(i,j)] )
                    # each point is weighted by 1/sigma
                    self.weights.append(weight / max(error,1.0e-2))
                    self.curve_ID.append( self.curve_counter )
        self.curve_names.append(curve_name)
        self.curve_counter += 1
    def _plot_gradVrep(self):
        """
        make a scatter plot of {Ri,V'rep(Ri)}
        """
        print self.curve_ID
        
        from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, savefig, title
        
        title("DATA for gradient of repulsive potential for %s" % self.name_atom_pair)
        xlabel("internuclear distance R [$\\AA$]")
        ylabel("$V^{\'}_{rep}(R)$ [eV/$\\AA$]")
        for i in range(0, self.curve_counter):
            indx = where(array(self.curve_ID) == i)
            Rarr = array(self.nuclear_separation)[indx]*bohr_to_angs
            plot(    Rarr, \
                     array(self.repulsive_gradient)[indx]*hartree_to_eV/bohr_to_angs, \
                     "o", label="%s" % self.curve_names[i])
        plot(Rarr, 0.0*Rarr)
        legend(loc='lower right')
#        savefig("/tmp/repulsive_gradient.png")
        show()
    def _integrate_Vrep(self, Rcut, nr_knots):
        """
        Vrep(R) = - integral_R^(R_cut) V'_rep(r)dr
        """
        self.R = linspace(0.0, Rcut, 100)
        gradVrep = self._fit_gradient(Rcut, nr_knots)
        self.Vrep = -array([integrate.quad(gradVrep, Ri, Rcut)[0] for Ri in self.R])
        print "repulsive potential:"
        print self.Vrep
        for i,Ri in enumerate(self.R):
            print "%s %s" % (Ri,self.Vrep[i])
        #
        from matplotlib.pyplot import plot, show, xlabel, ylabel, title, legend
        #title("repulsive potential for %s" % self.name_atom_pair)
        xlabel("internuclear distance r [$\\AA$]", fontsize=17)
        ylabel("$V_{rep}(r)$ [$eV$]", fontsize=17)
        plot(self.R*bohr_to_angs, self.Vrep*hartree_to_eV, lw=2, label="%s" % self.name_atom_pair.upper())
        legend(fontsize=17)
        show()
        #
    def _fit_gradient(self, Rcut, nr_knots):
        """
        find a curve that smoothly approximates the data in {Ri,V'rep(Ri)}

        Parameters:
        ===========
        Rcut: cutoff radius
        nr_knots: number of equidistant knots. The fewer knots, the smoother the interpolation.

        Returns:
        ========
        a callable functions that gives the gradient V'_rep(x)
        """
        R = array(self.nuclear_separation)
        dVrep = array(self.repulsive_gradient)
        weights = array(self.weights)
        sort_indx = argsort(R)
        print "weights = "
        print "weights.max = %s" % weights.max()
        print "weights.min = %s" % weights.min()
        print weights
        # interior knots
        # should have much less knots than data points
        nr_knots = min(nr_knots, max(1,len(R)/2))
        # left and right knots are added automatically, compute spacing
        # to leftmost and rightmost interior knots
        dR = (R.max()-R.min())/float(nr_knots+1)
        knots = linspace(R.min()+dR, R.max()-dR, nr_knots)
        smoothing_spline = interpolate.LSQUnivariateSpline(R[sort_indx], dVrep[sort_indx], knots, w=weights[sort_indx])
        #
        from matplotlib.pyplot import plot, errorbar, show, xlabel, ylabel, legend, title
        #title("FIT to gradient of repulsive potential for %s" % self.name_atom_pair)
        xlabel(r"internuclear distance r [$\AA$]", fontsize=17)
        ylabel(r"$\frac{d V_{rep}}{dr}$ [$eV / \AA$]", fontsize=17)
        for i in range(0, self.curve_counter):
            sel_indx = where(array(self.curve_ID) == i)
            plot(    array(self.nuclear_separation)[sel_indx]*bohr_to_angs, \
                     array(self.repulsive_gradient)[sel_indx]*hartree_to_eV/bohr_to_angs, \
                     "o", label="%s" % self.curve_names[i])

        errorbar(R[sort_indx]*bohr_to_angs, dVrep[sort_indx]*hartree_to_eV/bohr_to_angs, 1.0/weights[sort_indx]*hartree_to_eV/bohr_to_angs, ls="", color="grey", lw=0.5)
        Rarr = linspace(0.0, R.max(), 1000)
        print smoothing_spline(Rarr)
        plot(Rarr*bohr_to_angs, smoothing_spline(Rarr)*hartree_to_eV/bohr_to_angs, ls="-.", label="fit (w/ %d knots)" % nr_knots, lw=2, color="black")
        legend(loc='lower right', fontsize=17)
        show()
        #
        return smoothing_spline
    def fit(self, Rcut=3.0, nr_knots=10, reppot_dir="/tmp/"):
        """
        Parameters:
        ===========
        Rcut: cut-off radius in Angstrom
        nr_knots: controls the number of knots. The fewer knots, the smoother the function, \
the more knots the greater the agreement between data and fit.
        reppot_dir: directory where the modules for repulsive potential are located. \
Arrays with interatomic distance and repulsive potential will be stored in \
a file called <reppot_dir>/<atom 1>_<atom 2>.py.
        """
        Rcut /= bohr_to_angs
        self.add_cutoff_curve(linspace(Rcut, Rcut+2.0, 10), mean(self.weights))
        self._plot_gradVrep()
        self._integrate_Vrep(Rcut, nr_knots)
        self.write_py_reppot(reppot_dir)
    def write_py_reppot(self, reppot_dir):
        """
        write table for repulsive potential to a python module that can be loaded.

        Parameters:
        ===========
        reppot_dir: directory where the modules for repulsive potential are located.
           arrays with interatomic distance and repulsive potential will be stored in
           a file called
               <reppot_dir>/<atom 1>_<atom 2>.py.
        """
        atname1 = atom_names[self.Z1-1]
        atname2 = atom_names[self.Z2-1]
        if reppot_dir == "reppot_tables":
            script_dir = os.path.dirname(os.path.realpath(__file__))
            reppot_dir = os.path.join(script_dir, "reppot_tables/")
        reppot_file = os.path.join(reppot_dir, "%s_%s.py" % (atname1, atname2))

        import pprint
        import sys
        numpy.set_printoptions(threshold=sys.maxint)
        pp = pprint.PrettyPrinter(depth=10)
        fh = open(reppot_file, "w")
        print>>fh, "# This file has been generated automatically by %s" % sys.argv[0]
        print>>fh, "# The repulsive potential has been fitted to the following data sets:"
        for curve_name in self.curve_names:
            print>>fh, "#   - %s" % curve_name
        print>>fh, "from numpy import array"
        print>>fh, "Z1 = %s" % self.Z1
        print>>fh, "Z2 = %s" % self.Z2
        print>>fh, "# grid for distance d between atomic centers in bohr"
        print>>fh, "d = \\\n%s" % pp.pformat(self.R)
        print>>fh, "# repulsive potential in hartree/bohr"
        print>>fh, "Vrep = \\\n%s" % pp.pformat(self.Vrep)
        
        fh.close()

    
if __name__ == "__main__":
    import sys
    usage = "Usage: %s <element1> <element2>\n" % sys.argv[0]
    usage += "element1 and element2 are the elements in the atom pair for which the repulsive potential should be fitted (e.g. h and c)."
    if len(sys.argv) < 3:
        print usage
        exit(-1)
    parser = utils.OptionParserFuncWrapper(ReppotFitter.fit,usage)
    (options, args) = parser.parse_args()
    Z1, Z2 = atomic_number(args[0]), atomic_number(args[1])
    Fitter = ReppotFitter(Z1,Z2)

    print "Files with geometries and forces are read from stdin."
    print "Each line should be of the form:"
    print "<identifier>  <xyz file with geometry> <xyz file with electronic dftb forces> <xyz file with forces from different method> \\"
    print "                                    weight=<weight, positive float> max_error=<max. error per atom> \\"
    print "                                            (active_atoms=<list of atoms IDs, e.g. 0,1,2>)"
    for line in sys.stdin.readlines():
        if line.strip()[0] == "#":
            # ignore comments
            continue
        molname, geom_file, dftb_force_file, force_file, keywords_str = line.strip().split(None,4)
        molname = molname.replace("__", " ")
        keywords = dict(map(lambda s: tuple(s.replace(" ","").replace("\t","").split("=")), keywords_str.split()))
        weight = float(keywords.get("weight", 1.0))
        max_error = float(keywords.get("max_error", 0.01))
        active_atoms_str = keywords.get("active_atoms", "")
        if active_atoms_str == "":
            # include all atoms
            active_atoms = []
        else:
            active_atoms = map(int, active_atoms_str.split(","))

        print "weight = %s" % weight
        print "max_error = %s" % max_error
        geometries = XYZ.read_xyz(geom_file)
        forces_wo_frep = XYZ.read_xyz(dftb_force_file, units="hartree/bohr")
        forces_with_frep = XYZ.read_xyz(force_file, units="hartree/bohr")
        print "add geometries from %s" % geom_file
        Fitter.add_polyatomic_curve(geometries, forces_wo_frep, forces_with_frep, \
                                 curve_name=molname, weight=weight, max_error=max_error, active_atoms=active_atoms)
#    Fitter.fit(3.0, 30.0)
    Fitter.fit(**options)
