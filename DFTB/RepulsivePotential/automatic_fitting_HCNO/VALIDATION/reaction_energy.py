#!/usr/bin/env python
import sys
import os.path
"""
For a reaction

   m*A + n*B <---> k*C + ...     m,n,k integers

the reaction energy E(reaction) = E(products) - E(educts) should be determined.
The reaction equation is rewritten into

   m A  + n B  - k C ... <---> 0

so that the coefficients of the educts are positive and the coefficients of the products
are negative.

Example:

   %s   1 ethene  2 h2  -2 methane

loads the files ethene.xyz, h2.xyz and methane.xyz in the DFTB and REFERENCE folders
and computes the energy for the reaction

       ethene + 2*H2 <---> 2*methane

from the total energies given in the comment line of the xyz-files.
""" % os.path.basename(sys.argv[0])

from DFTB import XYZ, AtomicData

if __name__ == "__main__":
    import sys
    args = sys.argv[1:]
    
    educts = []
    coefs_educts = []
    energies_educts = []
    products = []
    coefs_products = []
    energies_products = []

    assert len(args) % 2 == 0, "For each molecule the stoichiometry and the name should be specified, e.g. '2 H2'"

    nmol = len(args)/2

    for i in range(0, nmol):
        coef = int(args[2*i])
        molecule = str(args[2*i+1])
        if coef < 0:
            products.append(molecule)
            coefs_products.append(coef)
        else:
            educts.append(molecule)
            coefs_educts.append(coef)

    # chemical equation
    educts_str   = " + ".join(["%d %s " % (coef, educt) for (coef,educt) in zip(coefs_educts, educts)])
    products_str = " + ".join(["%d %s " % (-coef, prod) for (coef,prod) in zip(coefs_products, products)])
    reaction_str =  educts_str + " <---> " + products_str
    
    # compute reference energy
    en_ref = 0.0
    for coef, mol in zip(coefs_educts + coefs_products, educts + products):
        en = float(XYZ.extract_keywords_xyz(os.path.join("REFERENCE", mol+".xyz"))["energy"])
        en_ref -= coef * en

    # compute DFTB energy
    en_dftb = 0.0
    for coef, mol in zip(coefs_educts + coefs_products, educts + products):
        en = float(XYZ.extract_keywords_xyz(os.path.join("DFTB", mol+".xyz"))["energy"])
        en_dftb -= coef * en
    
    # convert from Hartree to kcal/mol
    en_ref  *= AtomicData.hartree_to_kcalmol
    en_dftb *= AtomicData.hartree_to_kcalmol

    print "%40.40s         %+15.1f       %+15.1f" % (reaction_str, en_ref, en_dftb)
    

