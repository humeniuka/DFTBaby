#!/usr/bin/env python
"""
Generate Slater-Koster files from precomputed pseudo atoms. You have to regenerate
the tables if you change the pseudo atoms or confinement potentials.
"""
from SKIntegrals import AtomPair
import os.path
import numpy
from DFTB.AtomicData import atom_names

script_dir = os.path.dirname(os.path.realpath(__file__))
slako_dir = os.path.join(script_dir, "slako_tables/")

print "-load polar two-center grid"
try:
    from slako_tables.double_polar_grid import d, grid
except ImportError as e:
    print e
    raise Exception("Maybe you first have to generate a double polar grid for integration. Try: \n  python DFTB/SlaterKoster/generate_ptcgrid.py")
        
def slako_table(atom1,atom2):
    """write and plot Slater-Koster table for atom pair"""
    print "-compute Slater-Koster integrals for atom pair %s-%s" \
        % (atom_names[atom1.Z-1], atom_names[atom2.Z-1])
    dimer = AtomPair(atom1,atom2)

#    slako_dir = "/tmp/"
    dimer.getSKIntegrals(d, grid)
    dimer.write_py_slako(slako_dir)
    # disable plotting when working through ssh
#    dimer.plotSKIntegrals(slako_dir)
    

from confined_pseudo_atoms import \
       h,  he,                                            \
       li, be,                         b, c, n, o, f, ne, \
       na, mg,                         al,si,p, s, cl,ar, \
               sc, ti,  fe,  cu,  zn,              br,    \
                        ru,  ag

"""
# generate slako-tables for all possible combinations of pseudo atoms
from confined_pseudo_atoms import pseudo_atom_list
for at1 in pseudo_atom_list:
    slako_table(at1,at1)
    for at2 in pseudo_atom_list:
        if at1.Z < at2.Z:
            slako_table(at1,at2)
"""

if __name__ == "__main__":
    """
    # CHNO atoms for bio-organic molecules
    # H-H
    slako_table(h,h)

    # C-X
    slako_table(h,c)
    slako_table(c,c)
    slako_table(c,n)
    slako_table(c,o)

    # N-X
    slako_table(h,n)
    slako_table(n,n)
    slako_table(n,o)
    
    # O-X
    slako_table(h,o)
    slako_table(o,o)
    """
    
    """
    # sulphur and silver are contained in the ligand protected Ag_44(p-MBA)_30^(4-) cluster
    # What about the alkali counter ions that carry 4+ charges?
    # sulphur
    # S-X
    slako_table(h,s)
    slako_table(c,s)
    slako_table(n,s)
    slako_table(o,s)
    slako_table(s,s)

    # silver 
    # Ag-X
    slako_table(h,ag)
    slako_table(c,ag)
    slako_table(n,ag)
    slako_table(o,ag)
    slako_table(s,ag)
    slako_table(ag,ag)
    """
    """
    # Transition metals that are important in biomolecules
    # iron
    # Fe-X
    slako_table(h,fe)
    slako_table(c,fe)
    slako_table(n,fe)
    slako_table(o,fe)
    slako_table(s,fe)
    slako_table(ag,fe)
    slako_table(fe,fe)
    """
    """
    # sodium
    # Na-X
#    slako_table(h,na)
#    slako_table(c,na)
#    slako_table(n,na)
#    slako_table(o,na)
#    slako_table(na,na)
    """    
    """
    # silicon
    # Si-X
    slako_table(h,si)
    slako_table(c,si)
    slako_table(n,si)
    slako_table(o,si)
    slako_table(si,si)
    """
    """
    # zinc
    # Zn-X
    slako_table(h,zn)
    slako_table(c,zn)
    slako_table(n,zn)
    slako_table(o,zn)
#    slako_table(s,zn)
#    slako_table(ag,zn)
#    slako_table(fe,zn)
    slako_table(zn,zn)
    """
    """
    # Halogenes F, Cl and Br
    # F-X
    slako_table(h,f)
    slako_table(c,f)
    slako_table(n,f)
    slako_table(o,f)
    slako_table(f,f)
    #slako_table(f,s)
    #slako_table(f,zn)
    #slako_table(f,ag)

    # Cl-X
    slako_table(h,cl)
    slako_table(c,cl)
    slako_table(n,cl)
    slako_table(o,cl)
    slako_table(f,cl)
    #slako_table(s,cl)
    slako_table(cl,cl)
    #slako_table(cl,zn)
    #slako_table(cl,ag)
    """
    """
    # Br-X
    slako_table(h,br)
    slako_table(c,br)
    slako_table(n,br)
    slako_table(o,br)
    #slako_table(f,br)
    #slako_table(s,br)
    #slako_table(cl,br)
    #slako_table(zn,br)
    slako_table(br,br)
    #slako_table(br,ag)
    """
    """
    # Boron
    slako_table(h,b)
    slako_table(b,b)
    slako_table(c,b)
    slako_table(n,b)
    slako_table(o,b)
    # Phosphorous
    slako_table(h,p)
    slako_table(b,p)
    slako_table(c,p)
    slako_table(n,p)
    slako_table(o,p)
    slako_table(p,p)
    """
    """
    # Ruthenium 2+
    slako_table(h,ru)
    slako_table(c,ru)
    slako_table(n,ru)
    slako_table(o,ru)
    slako_table(ru,ru)
    """
    """
    # Silver
    slako_table(h,ag)
    slako_table(b,ag)
    slako_table(c,ag)
    slako_table(n,ag)
    slako_table(o,ag)
    slako_table(p,ag)
    slako_table(ag,ag)
    """
    """
    # Titanium - H,C,N,O
    slako_table(h,ti)
    slako_table(c,ti)
    slako_table(n,ti)
    slako_table(o,ti)
    slako_table(ti,ti)
    """
    """
    # Lithium
    slako_table(h,li)
    # Lithium dimer
    slako_table(li,li)
    """

    # Scandium - H,C,N,O,F,S
    slako_table(h,sc)
    slako_table(c,sc)
    slako_table(n,sc)
    slako_table(o,sc)
    slako_table(f,sc)
    slako_table(s,sc)
    slako_table(sc,sc)

    
