#!/usr/bin/env python
"""
To speed up SCF calculations we need a good guess for the initial
density. The radial density can be parametrized as 
n(r) = Nelec * a^3/(8*pi) * exp(-a*r)

Here we fit for each pseudo atom the exponent a 
to the converged density from a previous calculation.
For atoms with p-orbitals the fit gives a bad initial density.
"""

from numpy import *
from matplotlib.pyplot import *
from scipy.optimize import fmin

def fit_density(atom):
    """
    Find the exponent a that adjusts the parametrized
    density n(r) = Nelec * (a/pi)^(3/2) * exp(-a*r^2)
    best to the converged pseudo atom density of atom.

    Parameters:
    ===========
    atom: atom module, as written by generate_pseudoatoms.py
    """
    def n0(a,r):
        return atom.Nelec * pow(a, 3.0)/(8.0*pi) * exp(-a*atom.r)
    def deviation(a):
        return sum(abs(n0(a,atom.r) - atom.radial_density)**2)

    aopt = fmin(deviation, 1.0)[0]
 
    title("Z = %s, Nelec = %s" % (atom.Z, atom.Nelec))
    plot(atom.r, atom.radial_density, label="LDA density")
    plot(atom.r, n0(aopt, atom.r), label="fitted density")
    legend()
    show()

    return aopt

from pseudo_atoms import h,he,li,    c

density_exponents = {(1,1): fit_density(h),
                     (2,2): fit_density(he),
                     (3,3): fit_density(li),
                     (6,6): fit_density(c)}

print "density_exponents = %s" % density_exponents

    
