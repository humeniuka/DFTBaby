#!/usr/bin/env python
"""

Common atom types in the DREIDING force field:

   ID    name             description
   ===================================================
   0     H_        hydrogen
   1     H__HB     hydrogen capable of forming hydrogen bonds

   5     C_3       tetrahedral carbon (sp^3)
   6     C_R       sp^2 carbon involved in resonance (e.g. benzene)
   7     C_2       trigonal sp^2 carbon
   8     C_1       linear sp^1 carbon
   9     N_3       tetrahedral nitrogen
  10     N_R       nitrogen in resonance situtation (e.g. peptide bond)
  11     N_2       trigonal nitrogen
  12     N_1       linear nitrogen
  13     O_3       tetrahedral oxygen
  14     O_R       oxygen in resonance situtation
  15     O_2       trigonal 
  16     O_1       linear

  36     Zn        Zn^2+

There are problems with distinguishing trigonal (_2) from resonance (_R).

If a structure consists of many identical units (molecular aggregate) one can
first assign the atom types for a single monomer unit with the script

    assign_atom_types.py

and then copy the atom types to all other units in the large structure with

    copy_atom_types.py

"""
from DFTB import XYZ, AtomicData
from DFTB.Modeling import BondOrders
import numpy as np

atom_types2num = {"H_": 0, "H__HB": 1, "H_b": 2,
                  "B_3": 3, "B_2": 4,
                  "C_3": 5, "C_R": 6, "C_2": 7, "C_1": 8,
                  "N_3": 9, "N_R": 10, "N_2": 11, "N_1": 12,
                  "O_3": 13, "O_R": 14, "O_2": 15, "O_1": 16, "F_": 17,
                  "Al3": 18, "Si3": 19, "P_3": 20, "S_3": 21, "Cl": 22,
                  "Ga3": 23, "Ge3": 24, "As3": 25, "Se3": 26, "Br": 27,
                  "In3": 28, "Sn3": 29, "Sb3": 30, "Te3": 31, "I_": 32,
                  "Na": 33, "Ca": 34, "Fe": 35, "Zn": 36}

if __name__ == "__main__":
    import sys
    import os
    if len(sys.argv) < 3:
        print "Usage: %s .xyz-file .ff-file" % (os.path.basename(sys.argv[0]))
        print "  assigns atom types automatically for geometry in .xyz-file."
        print ""
        print "  This script simply adds a 5th column with the atom type"
        print "  and writes the result to the .ff-file ."
        print " "
        print "  WARNING: The type assignment and the force field implementation itself probably have lots of bugs!"
        exit(-1)
    # laod xyz-file
    xyz_file = sys.argv[1]
    atomlist = XYZ.read_xyz(xyz_file)[0]
    # find total charge
    kwds = XYZ.extract_keywords_xyz(xyz_file)
    charge = kwds.get("charge", 0.0)
    
    nat=len(atomlist)
    # determine bond orders
    print "connectivity matrix"
    ConMat = XYZ.connectivity_matrix(atomlist, search_neighbours=200)
    print "bond order assignment"
    bondsTuples, bond_orders, lone_pairs, formal_charges = BondOrders.assign_bond_orders(atomlist, ConMat, charge=charge)
    # list of bond orders 
    bond_orders_atomwise = [ [] for i in range(0, nat) ]
    # names of atoms connected to each 
    bonded_atomnames = [ [] for i in range(0, nat) ]
    for i,(atI,atJ) in enumerate(bondsTuples):
        BO = bond_orders[i]
        bond_orders_atomwise[atI].append( BO )
        bond_orders_atomwise[atJ].append( BO )
        # name of atoms involved in this bond
        Zi = atomlist[atI][0]
        Zj = atomlist[atJ][0]
        atnameI = AtomicData.atom_names[Zi-1].capitalize()
        atnameJ = AtomicData.atom_names[Zj-1].capitalize()
        bonded_atomnames[atI].append( atnameJ )
        bonded_atomnames[atJ].append( atnameI )
    print "atom type assignment"
    ff_file = sys.argv[2]
    fh = open(ff_file, "w")
    print>>fh, nat
    print>>fh, ""
    for i,(Zi,posi) in enumerate(atomlist):
        posi = np.array(posi) * AtomicData.bohr_to_angs
        atname = AtomicData.atom_names[Zi-1].capitalize()
        # list of bond orders for all bonds connected to atom i
        BOs = bond_orders_atomwise[i]
        # number of bonds
        nr_bonds = len(BOs)
        #
        # Type assignment in part taken from 'pysimm' package
        #
        atomtype = None
        # partial charges
        charge = 0.0

        if atname == "H":
            if (len(bonded_atomnames[i]) == 1) and bonded_atomnames[i][0] in ['N', 'O', 'F']:
                # hydrogen bonded to 1 electronegative atom  (N, O, F)
                atomtype = 'H__HB'
            else:
                atomtype = 'H_'
        elif atname == "B":
            if 2.0 in BOs:
                atomtype = "B_2"
            else:
                atomtype = "B_3"
        elif atname == "C":
            if 1.5 in BOs:
                atomtype = 'C_R'
            elif nr_bonds == 4:
                atomtype = 'C_3'
            elif nr_bonds == 3:
                atomtype = 'C_2'
            elif nr_bonds == 2:
                atomtype = 'C_1'
        elif atname == "N":
            if nr_bonds == 2 and (2.0 in BOs):
                atomtype = 'N_R'
            elif 2.0 in BOs:
                atomtype = 'N_2'
            elif 3.0 in BOs:
                atomtype = 'N_1'
            elif 1.0 in BOs:
                if nr_bonds == 3:
                    atomtype = 'N_3'
                elif nr_bonds == 2:
                    atomtype = 'N_2'
        elif atname == "O":
            if 1.5 in BOs:
                atomtype = 'O_R'
            elif 2.0 in BOs:
                if nr_bonds == 2:
                    atomtype = 'O_2'
                elif nr_bonds == 1:
                    atomtype = 'O_1'
            elif 1.0 in BOs:
                atomtype = 'O_3'
        elif atname in ["Cl", "Br", "Na", "Ca", "Fe", "Zn"]:
            atomtype = atname
        elif atname in ["F", "I"]:
            atomtype = atname + "_"
        elif atname in ["Al", "Si", "Ga", "Ge", "As", "Se", "In", "Sn", "Sb", "Te"]:
            atomtype = atname + "3"
        elif atname in ["P", "S"]:
            atomtype = atname + "_3"
            
        if atomtype == None:
                raise RuntimeError("Cannot assign atom type to '%s(%d)'" % (atname, i+1))

        type_num = atom_types2num[atomtype]

        print>>fh, " %s     %10.7f %10.7f %10.7f    %3.1d   %3.2f  %s" % (atname, posi[0], posi[1], posi[2], type_num, charge, atomtype)
    """
    # Lattice vectors
    lattice_constant = 100000.0  #  very large box
    print>>fh, "Tv	   %10.7f	  0.00000        0.00000" % lattice_constant
    print>>fh, "Tv	   0.00000	  %10.7f	 0.00000" % lattice_constant
    print>>fh, "Tv	   0.00000	  0.00000	 %10.7f " % lattice_constant
    fh.close()
    """
    print "Atom types written to 7th column of '%s'" % ff_file
    
