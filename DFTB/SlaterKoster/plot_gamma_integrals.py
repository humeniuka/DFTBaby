#!/usr/bin/env python
"""
compare gamma-integrals derived from numerical atomic valence orbitals with
those derived from Gaussian charge fluctuation functions
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

from DFTB.AtomicData import atomic_number
from DFTB.MolecularIntegrals.NumericalGammaApproximation import l2spec

from DFTB.Parameters import get_hubbard_parameters

import matplotlib as mpl
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size

if __name__ == "__main__":
    import sys
    import os.path
    
    if len(sys.argv) < 3:
        print "Usage: %s  <atom A> <atom B>" % os.path.basename(sys.argv[0])
        print "  plot gamma-integrals for atom combination A-B."
        print "  <atom A> and <atom B> should be the names of the atoms, i.e. 'h' or 'c' etc."
        exit(-1)

    atom_nameA = sys.argv[1].upper()
    atom_nameB = sys.argv[2].upper()
    Za = atomic_number(atom_nameA)
    Zb = atomic_number(atom_nameB)

    # plot gamma functions for
    # ... confined pseudo atoms 
    #from DFTB.SlaterKoster.confined_pseudo_atoms import gamma_integrals
    # ... free pseudo atoms
    from DFTB.SlaterKoster.free_pseudo_atoms import gamma_integrals

    try:
        gamma_dic = gamma_integrals.gamma_integrals[(Za,Zb)]
    except KeyError as e:
        print "ERROR: No numerical gamma-functions found for atom pair %s-%s." % (atom_nameA, atom_nameB)
        print "You have to add these atoms in the script 'generate_gamma_integrals.py' and recompute."
        raise e

    plt.xlabel(r"$r=\vert R_A - R_B \vert$ / bohr", fontsize=17)
    plt.ylabel(r"$\gamma$ integrals / Hartree", fontsize=17)
    
    # distances (in bohr) for which gamma_ab's are tabulated
    rg = gamma_integrals.r

    for (la,lb),gamma_ab in gamma_dic.iteritems():
        plt.plot(rg, gamma_ab, lw=2, label=r"$\gamma_{%s%s}(r)$ $(%s%s)$" % (atom_nameA, atom_nameB, l2spec[la], l2spec[lb]))

    # gamma-integrals for Gaussian fluctuation functions
    r = np.linspace(0.0001, 10.0, 200)
    
    hubbard_U = get_hubbard_parameters(None)
    Ua = hubbard_U[Za]
    Ub = hubbard_U[Zb]

    sigma_a = 1.0/(np.sqrt(np.pi) * Ua)
    sigma_b = 1.0/(np.sqrt(np.pi) * Ub)

    Cab = 1.0/np.sqrt(2*(sigma_a**2 + sigma_b**2))

    gamma_ab = erf(Cab*r)/r

    plt.plot(r, gamma_ab, lw=2, ls="-.", label=r"$\gamma_{%s%s}(r)$ (Gaussian)" % (atom_nameA, atom_nameB))

    # asymptotic limit
    r = np.linspace(2.0, 10.0, 100)
    plt.plot(r, 1/r, ls="--", label=r"$\frac{1}{r}$")
    
    plt.legend(fontsize=16)
    plt.tight_layout()
    plt.show()

    
                 
