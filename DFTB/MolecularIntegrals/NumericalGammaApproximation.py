#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The gamma approximation is a form of density fitting for approximating 2-center electron
integrals. Usually the charge fluctuation functions are taken to be Slater of Gaussian functions,
whose width is determined by the Hubbard parameter for the particular atom. 

The fluctuation functions can also be calculated directly from the atomic valence orbitals. 
Since the gamma-function should only depend on the distance, the charge density in the valence shell 
has to be made spherically symmetric by averaging over all m-values. 

 * l and m are the angular and magnetic quantum numbers of an atomic orbital

 * wfn^A_{n,l,m}(r) is a valence orbital of a confined pseudo atom A

The spherically symmetric charge fluctuation function for atom A and shell l is then given by

                 1                                      2
  F_{A,l}(r) = ----- sum_{m=-l,...,l} |wfn^A_{n,l,m}(r)|
               2*l+1

The gamma-function for two atoms A and B separated by a distance R=|R_A-R_B| is the following
Coulomb+xc integral
                      /     /                   1
  gamma (R)        =  | dr1 | dr2  F  (r1) { -------  + f  [rho0  + rho0 ] } F (r2)
       (A,lA),(B,lB)  /     /       A,lA     |r1-r2|     xc     A       B     B,lB

where rho0_A and rho0_B are the unperturbed atomic densities of atoms A and B (including core electrons). 
The integral is obtained numerically using Becke's scheme and tabulated for a range of distances 
for each atom combination A-B. At runtime gamma_AB(R) is interpolated using splines, for large R the 
asymptotic limit gamma_AB(R) --> 1/R is used. The numerical gamma-matrix obviates the need for
a Hubbard parameter.
But it means that the gamma-matrix depends not only on the atom type but also on the angular momenta
of the interacting orbitals.

References
----------
[1] Dominguez et.al., "Extensions of the Time-Dependent Density Functional Based Tight-Binding Approach",
    J. Chem. Theory Comput. 2013, 9, 4901−4914

""" 

from DFTB.BasisSets import AtomicBasisSet, AuxiliaryBasisSet
from DFTB.AtomicData import atomic_number
from DFTB.MolecularIntegrals.ERIs import electron_repulsion_integral_rho, AtomicDensitySuperposition
from DFTB.SlaterKoster import XCFunctionals

import numpy as np

def spherically_averaged_density(shell, l):
    """
    define a function for evaluating F_{A,l}(r) = 1/(2*l+1) sum_{m=-l,...,0,...,l} |wfn_{l,m}(r)|^2

    Parameters
    ----------
    shell      : list of instances of AtomicBasisFunction belonging to the same shell (same l)
    l          : angular momentum of shell

    Returns
    -------
    F          : callable, F(x,y,z) evaluates the spherically averaged density of the valence shell
    """
    def F(x,y,z):
        f = 0*x
        for bf in shell:
            assert bf.l == l
            # compute F_{A,l} = 1/(2*l+1) sum_{m=-l,...,0,...,l} |wfn_{l,m}(r)|^2
            amp = bf.amp(x,y,z)
            f += abs(amp)**2 / (2*l+1.0)
        return f
    
    return F


def charge_fluctuation_functions(atom_name='h', confined=True):
    """
    create charge fluctuation functions for valence shells

    Parameters
    ----------
    atom_name     :   string with atom name

    Optional
    --------
    confined      :   controls where confined atoms (True) or free atoms (False)
                      are used

    Returns
    -------
    ls            :   list of angular momenta of the valence shells
    Fs            :   list of charge fluctuation functions for each shell 
                      callable, Fs[i](x,y,z) evaluates the spherically averaged 
                      charged fluctuation function for the shell with angular
                      momentum ls[i]
    """
    Zat = atomic_number(atom_name)
    atomlist = [(Zat,(0,0,0))]
    # load atomic valence orbitals
    basis = AtomicBasisSet(atomlist, confined=confined)
    # group basis functions by angular momentum
    ls = []
    shells = []  # contains list of basis functions for each shell
    for bf in basis.bfs:
        if len(ls) == 0 or bf.l > ls[-1]:
            # new shell begins
            ls.append(bf.l)
            # create new shell with current basis function as first element
            shells.append( [bf] )
        else:
            # append basis function to current shell
            shells[-1].append(bf)
    # define charge fluctuation function for each shell
    Fs = []   # list of charge fluctuation functions    
    for l,shell in zip(ls,shells):
        assert len(shell) == 2*l+1
        Fs.append( spherically_averaged_density(shell, l) )
        
    return ls, Fs    


l2spec = {0: 's', 1: 'p', 2: 'd'}

def numerical_gamma_integrals(atom_nameA, atom_nameB, distances, confined=True):
    """
    compute the integrals

      gamma_{A,lA,B,lB} = (F_{A,lA}|1/r12 + f_xc[rho0A+rho0B]|F_{B,lB})

    numerically on a multicenter grid

    Parameters
    ----------
    atom_nameA, atom_nameB  :  names of interacting atoms, e.g. 'h' or 'c'
    distances               :  numpy array with interatomic separations (in bohr) for which
                               the gamma integrals should be calculated

    Optional
    --------
    confined      :   controls where confined atoms (True) or free atoms (False)
                      are used

    Returns
    -------
    gamma_dic               :  dictionary with values of gamma integrals on the distance grid,
                               gamma_dic[(lA,lB)] is a numpy array holding the integrals between
                               shell lA on atom A and shell lB on atom B
    """
    # atomic numbers
    Za = atomic_number(atom_nameA)
    Zb = atomic_number(atom_nameB)
    # charge fluctuation functions for each shell
    lsA, FsA = charge_fluctuation_functions(atom_nameA, confined=confined)
    lsB, FsB = charge_fluctuation_functions(atom_nameB, confined=confined)

    xc_functional = XCFunctionals.libXCFunctional("lda_x", "lda_c_pw")

    gamma_dic = {}
    for la,Fa in zip(lsA,FsA):
        for lb,Fb in zip(lsB,FsB):
            print "  integral between %s-shell and %s-shell" % (l2spec[la], l2spec[lb])
            gamma_ab = np.zeros(len(distances))
            for i,r_ab in enumerate(distances):

                print "   %3.d of %d   interatomic distance : %4.7f bohr" % (i+1,len(distances), r_ab)
                # atoms A and B are placed symmetrically on the z-axis,
                # separated by a distance of r_ab
                atomlist = [(Za,(0,0,-0.5*r_ab)),
                            (Zb,(0,0,+0.5*r_ab))]
                # define displaced charge fluctuation functions
                def rhoAB(x,y,z):
                    return Fa(x,y,z-0.5*r_ab)
                def rhoCD(x,y,z):
                    return Fb(x,y,z+0.5*r_ab)
                # 
                rho0 = AtomicDensitySuperposition(atomlist, confined=confined)
        
                # evaluate integral
                gamma_ab[i] = electron_repulsion_integral_rho(atomlist, rhoAB, rhoCD, rho0, xc_functional)

            # save gamma integral for the interaction of a shell with angular momentum la on atom A
            # and a shell with angular momentum lb on atom B
            gamma_dic[(la,lb)] = gamma_ab
            
    return gamma_dic

def tabulate_gamma_integrals(atom_names, filename, confined=True):
    """
    Compute the gamma integrals

      gamma_{A,lA,B,lB} = (F_{A,lA}|1/r12 + f_xc[rho0A+rho0B]|F_{B,lB})

    for a range of distances for each atom combination and save them to a python file
    that can be imported

    Parameters
    ----------
    atom_names  :   list of atom names, e.g. ['h','c',...]
    filename    :   path to a python file, where the integrals will be stored in a
                    human-readable form

    Optional
    --------
    confined      :   controls where confined atoms (True) or free atoms (False)
                      are used

    Returns
    -------
    nothing, data is written to a file
    """
    import pprint
    import sys
    
    # interatomic distances for which gamma-function is computed numerically,
    # values in between need to be interpolated
    distances = np.linspace(0.001, 6.0, 60)

    gamma_integrals_dic = {}
    # enumerate all unique atom combinations
    n = len(atom_names)
    for a in range(0, n):
        for b in range(a,n):
            print "computing gamma integrals for atom combination %s-%s" % (atom_names[a].upper(), atom_names[b].upper())
            gamma_dic = numerical_gamma_integrals(atom_names[a], atom_names[b], distances, confined=confined)

            Za = atomic_number(atom_names[a])
            Zb = atomic_number(atom_names[b])

            if Za > Zb:
                # swap atom, so that Za <= Zb
                tmp = Zb
                Zb = Za
                Za = tmp
            gamma_integrals_dic[(Za, Zb)] = gamma_dic

    fh = open(filename, "w")
    np.set_printoptions(threshold=sys.maxint)
    pp = pprint.PrettyPrinter(depth=10)
    print>>fh, "# This file has been generated automatically by %s." % sys.argv[0]
    print>>fh, "from numpy import array"
    print>>fh, ""
    print>>fh, "# atoms for which gamma-integrals are available"
    print>>fh, "atom_names = %s" % atom_names
    print>>fh, "# distances in bohr for which gamma-integrals are tabulated"
    print>>fh, "r = %s" % pp.pformat(distances)
    print>>fh, "# gamma-integrals gamma_{A,lA,B,lB} = (F_{A,lA}|1/r12 + f_xc[rho0A+rho0B]|F_{B,lB}) for each"
    print>>fh, "# atom combination A-B. The integrals are stored in a nested dictionary, the keys"
    print>>fh, "# on the first level are the atomic numbers (Za,Zb) with Za <= Zb, the keys on the"
    print>>fh, "# second level are the angular momenta of the valence shell on atom A and B, i.e. (la,lb)"
    print>>fh, "# For example the integral between the s-shell on carbon and the p-shell on nitrogen"
    print>>fh, "# would be stored in gamma_integrals[(6,7)][(0,1)]"
    print>>fh, "#                                      |      |   "
    print>>fh, "#                                   (Za,Zb) (lA,lB)"
    print>>fh, "gamma_integrals = \\\n%s" % pp.pformat(gamma_integrals_dic)

    print "gamma integrals for atoms %s written to '%s'" % (" ".join(atom_names), filename)
    fh.close()

#############################################################
#
# plotting and testing
#
#############################################################

def plot_charge_fluctuation_functions(atom_names, confined=True):
    """
    plot numerical charge fluctuation functions computed from valence orbital density
    """
    import matplotlib.pyplot as plt
    from DFTB.Parameters import get_hubbard_parameters
    
    plt.title("charge fluctuation functions")
    plt.xlabel("r / bohr", fontsize=17)
    plt.ylabel("numerical F(r)")
    
    r = np.linspace(0.01, 4.0, 200)
    x = r
    y = 0*r
    z = 0*r
    
    for atom_name in atom_names:
        lsA, FsA_numerical = charge_fluctuation_functions(atom_name, confined=confined)
        for la,Fa_numerical in zip(lsA, FsA_numerical):
            line, = plt.plot(r, Fa_numerical(x,y,z), lw=2, label=r"$F_{%s,%s}(r)$ (numerical)" % (atom_name.upper(), l2spec[la]))

    # compare with Gaussian fluctuation function
    for atom_name in atom_names:
        Zat = atomic_number(atom_name)
        hubbard_U = get_hubbard_parameters(None)
        # load Gaussian fluctuation functions, width is defined by the Hubbard parameter
        aux_basis = AuxiliaryBasisSet([(Zat, (0,0,0))], hubbard_U)
        # only the s-function is used
        Fa_gaussian = aux_basis.bfs[0]
        plt.plot(r, Fa_gaussian.amp(x,y,z), lw=2,
                 label=r"$F_{%s}(r)$ (Gaussian)" % atom_name.upper(),
                 ls="-.")
        
    plt.legend(ncol=2)
    plt.show()

def plot_numerical_gamma_integrals(atom_names, confined=True):
    import matplotlib.pyplot as plt

    plt.xlabel(r"$r = \vert R_A - R_B \vert$ / bohr")
    plt.ylabel(r"$\gamma_{A,l_A;B,l_B}(r)$")
    
    # interatomic distances for which gamma-function is computed numerically,
    # values in between need to be interpolated
    distances = np.linspace(0.01, 6.0, 20)
    
    # enumerate all unique atom combinations
    n = len(atom_names)
    for a in range(0, n):
        for b in range(a,n):
            print "computing gamma integrals for atom combination %s-%s" % (atom_names[a].upper(), atom_names[b].upper())
            gamma_dic = numerical_gamma_integrals(atom_names[a], atom_names[b], distances, confined=confined)

            for (la,lb),gamma_ab in gamma_dic.iteritems():
                plt.plot(distances, gamma_ab, lw=2, label=r"$\gamma_{%s%s,%s%s}(r)$" % (atom_names[a].upper(), l2spec[la], atom_names[b].upper(), l2spec[lb]))

    # asymptotic limit  1/r
    r = np.linspace(5.0, 15.0, 100)
    plt.plot(r, 1.0/r, label=r"$\frac{1}{r}$", lw=2, ls="-.")
    
    plt.legend()
    plt.show()

                    
if __name__ == "__main__":

    # compute on-site integrals for these atoms
    atom_names = ['h', 'c', 'n', 'o']
    #tabulate_gamma_integrals(atom_names, "/tmp/gamma_integrals.py", confined=True)
    #tabulate_gamma_integrals(atom_names, "/tmp/gamma_integrals.py", confined=False)

    # for confined pseudo atoms
    plot_charge_fluctuation_functions(atom_names, confined=True)
    #plot_numerical_gamma_integrals(atom_names, confined=True)
    # for free pseudo atoms
    plot_charge_fluctuation_functions(atom_names, confined=False)
    #plot_numerical_gamma_integrals(atom_names, confined=False)
    
