#!/usr/bin/env python
"""
   least square fitting of repulsive potentials

The repulsive energy (interaction between ions, core-electrons, E_rep = E_ref - E_elec) is approximated
as a sum of two-body potentials V_AB and one-body repulsive terms V_A, which ensure the correct energy
upon dissociation when r_ij -> infinity and V_AB(r_ij) -> 0.

    E_rep({r_atoms}) = sum_i<j V_AB(r_ij)  + sum_i V_A

where A,B are the atom types of atoms i,j. To fit the repulsive potentials, V_AB is decomposed as a 
linear combination

    V_AB(r) = sum_(n=1)^nmax f_AB,n(r) x_AB,n

of inverse cutoff-polynomials 

                 (r-r_cutoff)^2
    f_AB,n(r) =  --------------.
                   r^(2+n)

The coefficients x_AB,n and the one-body terms V_A are chosen such that the
sum of the squares is minimized:
                                                                                  2
    error =    sum      | sum_i<j V_AB(r) + sum_i V_A  - (E_ref^(s) - E_elec^(s)) |
            fit step s

subject to the constraint that the pair-potential V_AB is strictly repulsive, V_AB(r)' <= 0.

This is a linear least squares problem with linear inequality constraints 
that can be solved with standard linear algebra routines.

References:

    [1] Z.Bodrog, B. Aradi, Th. Fraunheim,
        "Automated Repulsive Parametrization for the DFTB Method",
        J. Chem. Theory Comput. 2011, 7, 2654-2664

    [2] M.Lourenco, M. da Silva, A. Oliveira, M. Quintao, H. Duarte
        "FASP: a framework for automation of Slater-Koster file parameterization"
        Theor. Chem. Acc. 2016, 135-250


"""
from DFTB import XYZ
from DFTB import AtomicData

import numpy
import numpy as np
import numpy.linalg as la
import scipy.optimize

import os.path

def atom_pair(Za,Zb):
    if Za <= Zb:
        return (Za,Zb)
    else:
        return (Zb,Za)

# Basis functions for representing V_AB(r)

class ReppotBasisFunction(object):
    def __init__(self, Za,Zb):
        self.AB = atom_pair(Za,Zb)
    def f(self, r):
        pass
    def dfdr(self, r):
        """derivative of f(r) w/r/t r"""
        pass
    def getAtomPair(self):
        return self.AB

# basis function for one-body term V_A
# f(r) = 1, the coefficient of Constant(1) is the one-body term V_A
class Constant(ReppotBasisFunction):
    def __init__(self, Za):
        self.AB = (Za,Za)
    def f(self, r):
        return 1.0
    def dfdr(self, r):
        return 0.0
    def __str__(self):
        Za,Za = self.AB
        return "Constant(%s,1)" % AtomicData.atom_names[Za-1]

# basis functions for pair-potentials V_AB(r)

# f(r) = (r-r_cutoff)^n, n >= 2

class PolyCutoff(ReppotBasisFunction):
    def __init__(self, Za, Zb, cutoff, degree):
        self.AB = atom_pair(Za,Zb)
        self.cutoff = cutoff
        self.degree = degree
    def f(self, r):
        if r < self.cutoff:
            V = (r-self.cutoff)**self.degree
        else:
            V = 0.0
        return V
    def dfdr(self, r):
        if r < self.cutoff:
            dVdr = self.degree * (r-self.cutoff)**(self.degree-1)
        else:
            dVdr = 0.0
        return dVdr
    def __str__(self):
        Za,Zb = self.AB
        return "PolyCutoff(pair=%s-%s, n=%d, cutoff=%s)" % (
            AtomicData.atom_names[Za-1], AtomicData.atom_names[Zb-1],
            self.degree, self.cutoff)

#         (r-r_cutoff)^2
# f(r) = ----------------   , n >= 1
#             r^(2+n) 

class InversePolyCutoff(ReppotBasisFunction):
    def __init__(self, Za, Zb, cutoff, degree):
        self.AB = atom_pair(Za,Zb)
        self.cutoff = cutoff
        self.degree = degree
    def f(self, r):
        if r < self.cutoff:
            V = (r-self.cutoff)**2 / (r**(2 + self.degree))
        else:
            V = 0.0
        return V
    def dfdr(self, r):
        if r < self.cutoff:
            dVdr = - (r-self.cutoff) / (r**(3 + self.degree)) * (self.degree * (r-self.cutoff) - 2*self.cutoff)
        else:
            dVdr = 0.0
        return dVdr
    def __str__(self):
        Za,Zb = self.AB
        return "InversePolyCutoff(pair=%s-%s, n=%d, cutoff=%s)" % (
            AtomicData.atom_names[Za-1], AtomicData.atom_names[Zb-1],
            self.degree, self.cutoff)


    
class Basis(object):
    def __init__(self, unique_atomtypes, cutoffs, max_degree=8):        
        """
        build a basis for the repulsive potentials
        
        Parameters
        ----------
        unique_atomtypes: list of atomic numbers Z which should be covered
        by the basis functions
        cutoffs: dictionary with cutoffs in bohr for each atom pair (Za,Zb)
        
        Returns
        -------
        bfs: list of basis functions f_AB,nu
        """
        self.unique_atomtypes = unique_atomtypes
        self.cutoffs = cutoffs
        self.atom_pairs = []
        self.basis_pairs = {} # basis functions by keys (Za,Zb)
        self.basis_functions = []   # list of basis functions
        # The coefficient C(i) of basis function i has to lie in the range
        # lower_bounds[i] <= C(i) <= upper_bounds[i]
        self.lower_bounds = [] 
        self.upper_bounds = []
        # for each key (Za,Zb) a list with the distances
        # present in the fit paths
        self.rij_points = {}
        
        Ntype = len(unique_atomtypes)
        unique_atomtypes = np.sort(unique_atomtypes)

        # one-body terms V_A
        for Za in unique_atomtypes:
            bf = Constant(Za)
            # basis function for V_A = u * Constant(1)
            self.basis_functions.append(bf)
            # no bounds on coefficients of V_A
            self.lower_bounds.append(-np.inf)
            self.upper_bounds.append(+np.inf)
            
        # two-body potentials V_AB(r)
        for A in range(0, Ntype):
            Za = unique_atomtypes[A]
            for B in range(A, Ntype):
                Zb = unique_atomtypes[B]
                assert Za <= Zb

                ab_type = atom_pair(Za,Zb)
                self.atom_pairs.append( ab_type )
                # basis functions for V_AB(r) = sum_nu f_AB,nu(r)
                bfs = []
                for degree in range(2, max_degree):
                    #bf = PolyCutoff(Za,Zb, self.cutoffs[(Za,Zb)], degree)
                    bf = InversePolyCutoff(Za,Zb, self.cutoffs[(Za,Zb)], degree)
                    bfs.append(bf)
                    
                    # bounds for InversePolyCutoff
                    if bf.__class__ == InversePolyCutoff and degree >= max_degree - 3:
                        self.lower_bounds.append(0.0)
                        self.upper_bounds.append(+np.inf)
                    else:
                        self.lower_bounds.append(-np.inf)
                        self.upper_bounds.append(+np.inf)
                        
                self.basis_functions += bfs
                self.basis_pairs[ab_type] = bfs
                self.rij_points[ab_type] = []
                
    def structure_constants(self, atomlist):
        """
        compute the s-th row of the matrix of structure constants A

                         /  sum_(i<j,ij~AB) f_(AB,k)(r_ij)                 if     k >= 1
           A_s,(AB,k) = <
                         \   number of A-atoms in geometry                 if     k = 0
        """
        Nat = len(atomlist)
        Nbf = len(self.basis_functions)
        
        Ae = np.zeros( (1, Nbf) )         # 1 row for each energy
        Af = np.zeros( (3*Nat, Nbf) )     # 3*Nat rows for each force vector
        # one-body terms, for each atom of type A, we add 1*V_A.
        for Zi,posi in atomlist:
            for k,bf in enumerate(self.basis_functions):
                if bf.__class__ == Constant and bf.getAtomPair() == (Zi,Zi):
                    Ae[0,k] += bf.f(0.0)
                    # one-body terms have no spatial depencence, therefore
                    # Af does not change
                    
        # two-body terms
        # iterate over all distances r_ij
        for i in range(0, Nat):
            Zi,posi = atomlist[i]
            for j in range(i+1,Nat):
                Zj,posj = atomlist[j]
                rij_vec = np.array(posj) - np.array(posi)
                rij = la.norm(rij_vec)

                # unit vector
                unit_vec_ij = rij_vec / rij
                
                ij_type = atom_pair(Zi,Zj)

                r_cutoff = self.cutoffs[ij_type]
                if (rij > r_cutoff):
                    # distances outside the cutoff range are not considered
                    continue
                
                for k,bf in enumerate(self.basis_functions):
                    if bf.__class__ != Constant and bf.getAtomPair() == ij_type:
                        Ae[0,k] += bf.f(rij)
                        
                        # gradient of basis function
                        dfdr = bf.dfdr(rij)
                        # force on atom i
                        Af[3*i:3*(i+1), k] +=   dfdr * unit_vec_ij
                        # force on atom j
                        Af[3*j:3*(j+1), k] += - dfdr * unit_vec_ij
                # record available distances on the i-j potential curce
                self.rij_points[ij_type].append( rij )

        return Ae, Af

    def inequality_constraint(self, Npts):
        """
        build the matrix C for the inequality constraint
        
              C.x <= 0

        that ensures that the derivative of the repulsive potential 
        is non-positive at Npts equidistance points in the interval [0.001,r_cutoff].


        Parameters
        ----------
        Npts: number of points at which the inequality constraint should be enforced

        Returns
        -------
        C: Nab*Npts x Nbf matrix, where Nab is the number of distinct atom pairs
           and Nbf is the number of basis functions
        """
        Nab = len(self.atom_pairs)
        Nbf = len(self.basis_functions)
        C = np.zeros( (Nab*Npts, Nbf) )
        # Each atom pair AB has a different cutoff radius, so the division
        # into Npts points is different
        rs_dic = {}
        for (Z1,Z2) in self.atom_pairs:
            rs_dic[(Z1,Z2)] = np.linspace(0.4, self.cutoffs[(Z1,Z2)], Npts)

        mu = 0 # enumerate rows of C
        for (Z1,Z2) in self.atom_pairs:
            for i in range(0, Npts):
                for k,bf in enumerate(self.basis_functions):
                    if (Z1,Z2) == bf.getAtomPair():
                        C[mu,k] = bf.dfdr(rs_dic[(Z1,Z2)][i])
                mu += 1

        return C
        
    def evaluate_reppot(self, coefs, r, Za,Zb):
        Uab = 0.0

        AB_type = atom_pair(Za,Zb)

        for c_k, bf_k in zip(coefs, self.basis_functions):
            # leave out one-body terms, we are only interested in V_AB(r)

            # two-body terms
            if bf_k.__class__ != Constant and bf_k.getAtomPair() == AB_type:
                Uab += c_k * bf_k.f(r)

        return Uab

    def one_body_term(self, coefs, Za):
        Ua = 0.0
        
        for c_k, bf_k in zip(coefs, self.basis_functions):
            # filter out one-body terms, there should be only one, actually
            if bf_k.__class__ == Constant and bf_k.getAtomPair() == (Za,Za):
                Ua += c_k * bf_k.f(0.0)

        return Ua
    
    def evaluate_reppot_array(self, coefs, rs, Za,Zb):
        Uab = np.zeros(len(rs))
        for ir,r in enumerate(rs):
            Uab[ir] = self.evaluate_reppot(coefs, r, Za,Zb)
        return Uab
    
def get_unique_atomtypes(fit_dir):
    """
    given a list of xyz-files, determine the atomic numbers that are present
    in any file. For each atom combination a repulsive potential is fitted.
    """
    
    xyz_filenames = glob.glob("%s/FIT_PATHS/*.xyz" % fit_dir)
    
    unique_atomtypes = []
    for f in xyz_filenames:
        for atomlist in XYZ.read_xyz_it(f):
            # atomic numbers in first frame
            for Zi,posi in atomlist:
                if not Zi in unique_atomtypes:
                    unique_atomtypes.append( Zi )

    unique_atomtypes = np.sort(unique_atomtypes)
    
    print ""
    print "Atom Types:"
    print "-----------"
    for Zi in unique_atomtypes:
        print "  %s" % AtomicData.atom_names[Zi-1]
    print ""
    
    return unique_atomtypes

            
class Fitter:
    def __init__(self, fit_dir):
        """
        The directory 'fit_dir' should contain a specific folder structure:

             
              
             ./FIT_PATHS/                           #  contains geometries of fit paths 
                          ethane.xyz
                          ethene.xyz
                          ethyne.xyz
                          ...
             ./REFERENCE/                           #  contains total reference energies and forces
                          ethane.energies.dat       #  for each fit path
                          ethane.forces.xyz
                          ethene.energies.dat
                          ethene.forces.xyz
                          ...
             ./DFTB/                                #  contains electronic DFTB energies and forces
                          ethane.energies.dat       #  for each fit path
                          ethane.forces.xyz

                          
        """
        self.fit_dir = fit_dir
        self.unique_atomtypes = get_unique_atomtypes(self.fit_dir)

        # H C N O
        cutoffs = { (1,1) : 1.3 / AtomicData.bohr_to_angs,
                    (1,6) : 2.1 / AtomicData.bohr_to_angs,
                    (6,6) : 2.3 / AtomicData.bohr_to_angs,
                    (1,7) : 1.3 / AtomicData.bohr_to_angs,
                    (6,7) : 2.3 / AtomicData.bohr_to_angs,
                    (7,7) : 2.3 / AtomicData.bohr_to_angs,
                    (1,8) : 2.1 / AtomicData.bohr_to_angs,
                    (6,8) : 2.3 / AtomicData.bohr_to_angs,
                    (7,8) : 2.3 / AtomicData.bohr_to_angs,
                    (8,8) : 2.3 / AtomicData.bohr_to_angs,}
        # H O Br
        cutoffs[(1 ,35)] = 2.0 / AtomicData.bohr_to_angs
        cutoffs[(8 ,35)] = 2.5 / AtomicData.bohr_to_angs
        cutoffs[(35,35)] = 2.5 / AtomicData.bohr_to_angs
        # H C N O Si
        cutoffs[(1 ,14)] = 2.0 / AtomicData.bohr_to_angs
        cutoffs[(6 ,14)] = 2.5 / AtomicData.bohr_to_angs
        cutoffs[(7 ,14)] = 2.5 / AtomicData.bohr_to_angs
        cutoffs[(8 ,14)] = 2.5 / AtomicData.bohr_to_angs
        cutoffs[(14,14)] = 2.5 / AtomicData.bohr_to_angs
        # H C N O Ru
        cutoffs[(1 ,44)] = 1.8 / AtomicData.bohr_to_angs
        cutoffs[(6 ,44)] = 2.2 / AtomicData.bohr_to_angs
        cutoffs[(7 ,44)] = 2.2 / AtomicData.bohr_to_angs
        cutoffs[(8 ,44)] = 2.2 / AtomicData.bohr_to_angs
        cutoffs[(44,44)] = 4.0 / AtomicData.bohr_to_angs
        
        self.basis = Basis(self.unique_atomtypes, cutoffs)

        print "Basis functions:"
        print "----------------"
        for bf in self.basis.basis_functions:
            print bf
        print ""
        
    def fit(self, targets=["energies"]): #, "forces"]):
        """
        compute the matrix of structure constants `X` and the error vector `E`
        from all fit steps and solve the linear least squares problem that yields
        the optimum linear combination `A` of repulsive basis functions,

              min |X.A - E|^2 
               A

        The vector on the right hand side contains the differences between the 
        reference and the electronic DFTB energies:

          E[s] = (E_ref^(s) - E_elec^(s))

        The matrix element X[s,i] constains the structure factors for the step s
        and the basis function k, which belongs to the atom pair A,B. 
        The sum extends over all distances between atoms i<j, where the atom types are A,B. 

          X[s,k] = sum_(i<j == AB(k)) f_k(r^(s)_ij)

        A[k] contains the coefficients of the basis functions and we want to find the
        linear combination that minimizes the difference

           error = sum_s (sum_k X[s,k] A[k]  -  E[s])^2

        """
        A = []
        y = []
        # load steps of all fit paths
        for atomlist, en_rep, force_rep in self.load_fitpaths_it(self.fit_dir):

            Ae,Af = self.basis.structure_constants(atomlist)

            if "energies" in targets:
                A.append(Ae)
                y.append( [en_rep] )
            if "forces" in targets:
                A.append(Af)
                y.append( force_rep )

        self.A = np.vstack(A)
        self.y = np.hstack(y)

        # matrix C, such that C.x <= 0 ensures the repulsive potentials are monotonically
        # decreasing functions
        Npts = 50  # number of points r_mu at which V_AB(r_mu) <= 0
        self.C = self.basis.inequality_constraint(Npts)

        """
        # solve constrained linear least square problem, only some coefficients are
        # are required to be positive
        opt_res = scipy.optimize.lsq_linear(self.A, self.y, 
                        bounds=(self.basis.lower_bounds, self.basis.upper_bounds), method="bvls")
        self.x = opt_res.x
        """
        
        # solve constrained linear least squares problem, d/dr V_AB(r) is required to be positive
        # at 20 equidistant points between points r=0 and r=r_cutoff
        self.x = lsqlin(self.A, self.y, self.C)

        avg_error = la.norm(np.dot(self.A, self.x) - self.y)/len(self.y)
        max_error = abs(np.dot(self.A, self.x) - self.y).max()
        
        print "average error: %e" % avg_error
        print "maximal error: %e" % max_error

        # show one-body terms
        print "One-body terms (Hartree):"
        for Za in self.basis.unique_atomtypes:
            atname = AtomicData.atom_names[Za-1]
            Ua = self.basis.one_body_term(self.x, Za)
            print "  V_%s = %7.5f" % (atname, Ua) 
            
    def load_fitpaths_it(self, fit_dir):
        """
        load geometries of all fit paths from the folder `fit_dir`
        and return an iterator to the fit steps. Each fit step consists of the tuple

          (atomlist, en_rep, force_rep)
           
        The meaning of the elements of a fit step are as follows:

          atomlist : geometry of of fit step, list of tuples (Zi,[xi,yi,zi])
          en_rep   : energy difference (E_ref - E_elec) in Hartree
          force_rep: vector with force difference (F_ref - F_elec) in Hartree/bohr,
                     force_rep[3*i:3*(i+1)] is the force [Fx(i),Fy(i),Fz(i)] on atom i.
        
        """
        self.names = [os.path.basename(xyz)[:-4] for xyz in glob.glob("%s/FIT_PATHS/*.xyz" % fit_dir)]

        for name in self.names:
            print "Loading fit path '%s' " % name
            geom_xyz = os.path.join(fit_dir, "FIT_PATHS", "%s.xyz" % name)
        
            en_dftb_dat = os.path.join(fit_dir, "DFTB", "%s.energies.dat" % name)        
            en_ref_dat = os.path.join(fit_dir, "REFERENCE", "%s.energies.dat" % name)
            
            forces_dftb_xyz = os.path.join(fit_dir, "DFTB", "%s.forces.xyz" % name)
            forces_ref_xyz = os.path.join(fit_dir, "REFERENCE", "%s.forces.xyz" % name)

            try:
                geometries = XYZ.read_xyz(geom_xyz)
                
                energies_dftb = np.loadtxt(en_dftb_dat)
                energies_ref = np.loadtxt(en_ref_dat)
                
                forces_dftb = XYZ.read_xyz(forces_dftb_xyz)
                forces_ref = XYZ.read_xyz(forces_ref_xyz)
                
                for atomlist, en_dftb, en_ref, f_dftb, f_ref in zip(
                        geometries, energies_dftb, energies_ref, forces_dftb, forces_ref):
                    # compute E_rep = E_ref - E_el
                    en_rep = en_ref - en_dftb
                    # compute E_rep = F_ref - F_el
                    force_rep = XYZ.atomlist2vector(f_ref) - XYZ.atomlist2vector(f_dftb)
                    
                    yield atomlist, en_rep, force_rep
                
            except IOError as e:
                print "Loading of fit path '%s' failed!" % name
        
    def show_fit(self):
        """
        plot the repulsive potentials for each atom combination that is 
        present in the fit paths.
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        mpl.rcParams['xtick.labelsize'] = 17
        mpl.rcParams['ytick.labelsize'] = 17
        
        fs = 18.0

        #plt.title("repulsive potential $V_{rep}(r)$", fontsize=fs)
        plt.xlabel("$r / \AA$", fontsize=fs)
        plt.ylabel("$V_{AB}$ / Hartree", fontsize=fs)
        plt.xlim((0.2, 5.0))
        plt.ylim((-0.03, 0.5))
        
        rs = np.linspace(0.1, 10.0, 3000)
        for (Za,Zb) in self.basis.atom_pairs:
            if len(self.basis.rij_points[(Za,Zb)]) == 0:
                # None the fit paths contained the atom combination Za-Zb, so
                # the repulsive potential for this atom pair could not be determined
                continue
            label = "%s-%s" % (AtomicData.atom_names[Za-1].capitalize(), AtomicData.atom_names[Zb-1].capitalize())
            print "%s" % label
            # evaluate Uab on equidistant grid
            Uab = self.basis.evaluate_reppot_array(self.x, rs, Za,Zb)
            l = plt.plot(rs*AtomicData.bohr_to_angs,
                         Uab, label=label)
            # evaluate Uab at rij's
            ab_type = atom_pair(Za,Zb)
            Uab = self.basis.evaluate_reppot_array(self.x, self.basis.rij_points[ab_type], Za,Zb)
            rs_pts = np.array(self.basis.rij_points[ab_type])
            plt.plot(rs_pts*AtomicData.bohr_to_angs,
                     Uab, "o", color=l[0].get_color(), alpha=0.3)

        plt.legend(ncol=3, fontsize=20)
        plt.show()

    def save_reppots(self, reppot_dir="REPPOT_TABLES"):
        """
        write table for repulsive potential to a python module that can be loaded.

        Parameters:
        ===========
        reppot_dir: directory where the modules for repulsive potential are located.
           arrays with interatomic distance and repulsive potential will be stored in
           a file called
               <reppot_dir>/<atom 1>_<atom 2>.py.
        """
        Npts = 500
        for (Z1,Z2) in self.basis.atom_pairs:
            atname1 = AtomicData.atom_names[Z1-1]
            atname2 = AtomicData.atom_names[Z2-1]
            if len(self.basis.rij_points[(Z1,Z2)]) == 0:
                # None the fit paths contained the atom combination Z1-Z2, so
                # the repulsive potential for this atom pair could not be determined
                print "Atom combination %s-%s not present in any fit path" % (atname1, atname2)
                continue
            dmin = min(self.basis.rij_points[(Z1,Z2)])
            d = np.linspace(dmin, self.basis.cutoffs[(Z1,Z2)], Npts)

            reppot_file = os.path.join(reppot_dir, "%s_%s.py" % (atname1, atname2))

            import pprint
            import sys
            numpy.set_printoptions(threshold=sys.maxint)
            pp = pprint.PrettyPrinter(depth=10)
            fh = open(reppot_file, "w")
            print>>fh, "# This file has been generated automatically by %s" % sys.argv[0]
            print>>fh, "# The repulsive potential has been fitted to the following data sets:"
            for name in self.names:
                print>>fh, "#   - %s" % name
            print>>fh, "from numpy import array"
            print>>fh, "Z1 = %s" % Z1
            print>>fh, "Z2 = %s" % Z2
            print>>fh, "# grid for distance d between atomic centers in bohr"
            print>>fh, "d = \\\n%s" % pp.pformat(d)
            print>>fh, "# repulsive potential in hartree/bohr"
            Vrep = self.basis.evaluate_reppot_array(self.x, d, Z1,Z2)
            print>>fh, "Vrep = \\\n%s" % pp.pformat(Vrep)
        
            fh.close()

            print "saved repulsive potential %s_%s.py to %s" % (atname1,atname2,reppot_dir)

def lsqlin(A,y, C):
    """
    solve the constrained linear least square problem

      min 1/2 || A.x - y ||^2   subject to   C.x <= 0
       x

    using the SLSQP algorithm (sequential least squares programming) as implemented
    in scipy.optimize.
    """
    # define objective function
    def f(x):
        """ 1/2 || A.x - y ||^2 """
        fx = 0.5 * la.norm(np.dot(A,x) - y)**2
        # The objective function is badly scaled. To avoid the error "Positive directional derivative for linesearch"
        # in the optimizer one has to rescale f(x) and f'(x) so that the values fall into a smaller range.
        fx /= 10.0
        return fx
    # jacobian of objective function
    def fprime(x):
        """ df/dx = A^T.(A.x-y) """
        dfdx = np.dot(A.transpose(), np.dot(A,x) - y )
        dfdx /= 10.0
        return dfdx
    # inequality constraints
    def f_ieqcons(x):
        return -np.dot(C,x)
    def fprime_ieqcons(x):
        return -C
    # initial guess
    x0 = np.zeros( A.shape[1] )
    
    xopt,fx,its,imode, smode = scipy.optimize.fmin_slsqp(
        f, x0, fprime=fprime,
        f_ieqcons=f_ieqcons, fprime_ieqcons=fprime_ieqcons,
        iprint=2,
        iter=1000,
        full_output=True)
    print "exit mode of optimizer: %s" % imode
    assert imode == 0
    # check that constraint C.x <= 0 is fulfilled
    print "C.x"
    print np.dot(C, xopt)
    #assert np.all(np.dot(C, xopt) <= 0.0)
    
    
    return xopt
    
    
            
if __name__ == "__main__":
    import glob
    import sys
    import os.path
    
    if len(sys.argv) < 2:
        print "Usage: %s  <fit directory>" % os.path.basename(sys.argv[0])
        print ""
        exit(-1)
    
    fit_dir = sys.argv[1]

    fitter = Fitter(fit_dir)
    fitter.fit()
    fitter.show_fit()
    fitter.save_reppots(reppot_dir="REPPOT_TABLES")
    
