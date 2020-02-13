#!/usr/bin/env python
"""
compute the unique on-site (i.e. 1-center) electron integrals for each atom
using the valence s- and p-orbitals of the confined or free pseudoatoms
"""
from DFTB.MolecularIntegrals.OnSiteIntegrals import tabulate_onsite_integrals
import os.path

if __name__ == "__main__":

    # compute on-site integrals for these atoms
    atom_names = ['h', 'c', 'n', 'o']

    script_dir = os.path.dirname(os.path.realpath(__file__))

    # for confined pseudo atoms
    confined_orbdir = os.path.join(script_dir, "confined_pseudo_atoms/") 
    tabulate_onsite_integrals(atom_names, confined_orbdir+"onsite_electron_integrals.py", confined=True)
    # and for free pseudo atoms
    free_orbdir = os.path.join(script_dir, "free_pseudo_atoms/") 
    tabulate_onsite_integrals(atom_names, free_orbdir+"onsite_electron_integrals.py", confined=False)
    
