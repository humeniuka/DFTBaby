#!/usr/bin/env python
"""
Compute the gamma integrals

      gamma_{A,lA;B,lB} = (F_{A,lA}|1/r12 + f_xc[rho0A+rho0B]|F_{B,lB})

numerically on a multicenter grid for atom types A-B. The charge fluctuation
functions F_{A,lA}(r) are derived from the numerical atomic orbitals. 
F_{A,lA} is taken to be the spherically averaged density of the 
valence shell (with angular momentum lA and 2*lA+1 degenerate orbitals)

                     1                                 2
      F_{A,l}(r) = ----- sum_{m=-l,..,l} |phi_{l,m}(r)|
                   2*l+1
"""
from DFTB.MolecularIntegrals.NumericalGammaApproximation import tabulate_gamma_integrals
import os.path

if __name__ == "__main__":

    # compute gamma integrals for these atoms
    atom_names = ['h', 'c', 'n', 'o']

    script_dir = os.path.dirname(os.path.realpath(__file__))

    # gamma integrals for confined pseudo atoms
    #confined_orbdir = os.path.join(script_dir, "confined_pseudo_atoms/")
    #tabulate_gamma_integrals(atom_names, confined_orbdir+"gamma_integrals.py", confined=True)

    # and for free pseudo atoms
    free_orbdir = os.path.join(script_dir, "free_pseudo_atoms/")
    tabulate_gamma_integrals(atom_names, free_orbdir+"gamma_integrals.py", confined=False)
