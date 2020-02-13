#!/usr/bin/env python
"""
assign bond orders to a molecular structure only based on the atom types and the connectivity
"""
import numpy as np
from scipy import optimize
import os

from DFTB import XYZ, AtomicData


def is_metal_ion(Z):
    if Z in [12, 28, 30]: # magnesium, nickel, zinc
        return True
    else:
        return False

def solve_linprog(c, A_ub, b_ub, A_eq, b_eq, bounds):
    """
    solve a integer linear program either using scipy's Simplex algorithm
    or the external library lpsolve
    """
    try:
        # use lpsolve
        from lp_maker import lp_maker, lpsolve
        # convert to format expected by lpsolve

        c = c.astype(int).tolist()  # c.x is the linear objective function
        a = np.vstack((A_eq, A_ub)).astype(int).tolist() # matrix which contains both equality and inequality constraints, A_eq.x == b_eq, A_ub.x <= b_ub
        b = np.hstack((b_eq, b_ub)).astype(int).tolist() # right hand side of equality and inequality constraints
        e = np.hstack(([0 for bi in b_eq], [-1 for bi in b_ub])).tolist()
        nvars = len(b_eq) + len(b_ub) # number of integer variables xi
        xint = [i for i in range(1, nvars+1)]
        vlb = [bound[0] for bound in bounds]
        vub = [bound[1] for bound in bounds]

        # build lpsolve object
        lp = lp_maker(c, a, b, e, vlb, vub, xint)
        lpsolve('solve', lp)
        lpsolve('get_objective', lp)
        X = lpsolve('get_variables', lp)[0]
        lpsolve('delete_lp', lp)
        return X
    except ImportError as e:
        # use scipy's linprog
        res = optimize.linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, options={"maxiter": 10000000, "disp": True})
        X = res.x
        print "WARNING: simplex algorithm may not respect integer constraints! Some bond orders may not be integers."
        return X
        
    
def assign_bond_orders(atomlist, ConMat, charge=0, latte_file=""):
    """
    Determine bond orders only based on element types and connectivity.
    The problem of distributing the electrons over the bonds in an optimal fashion
    is formulated as a integer linear program following the article

    Froeyen,~M. Herdewijn,P.
    'Correct Bond Order Assignment in a Molecular Framework Using Integer Linear
    Programming with Application to Molecules Where Only Non-Hydrogen Atom Coordinates
    Are Available'
    J. Chem. Inf. Model., 2005, 45, 1267-1274.

    The optimal solution is found using scipy's optimize.linprog function.
    The number of feasible solutions, which equals the number of lattice points inside
    the polytope defined by the constraints on electron counts, can be counted by
    the Latte program.

    Parameters:
    ===========
    atomlist: list of Nat tuples (Zi,[xi,yi,zi]) with atomic number and positions
    ConMat: connectivity matrix, if atom I is bonded to atom J, then ConMat[I,J] = 1

    Results:
    ========
    bonds: list of tuples (I,J) for each bond
    bond_orders: list of integers with bond orders (1=single, 2=double, 3=triple)
      for each bond. Since the Simplex algorithm cannot handle integer constraints,
      some bond orders might not be integers, e.g. in C60 some bond orders are found
      to be 1.5 .
    lone_pairs: list with the number of free electron pairs on each atom
    formal_charges: list of formal charges for each atom
    """
    # Xij is the number of bonds between atoms i and j
    # Xii is the number of lone pairs on atom i
    
    Nat = len(atomlist)
    # octet electrons
    OCT = 0.0
    # number of valence electrons
    V = 0.0 - charge
    ks = []
    valence_electrons = []
    for i,(Zi,posi) in enumerate(atomlist):
        atname = AtomicData.atom_names[Zi-1]
        # octet electrons
        if (Zi == 1):
            #
            ki = 2
        elif (Zi == 30):
            # Zinc 2+
            ki = 10
        else:
            ki = 8
        ks.append(ki)
        OCT += ki
        # valence electrons
        vi = AtomicData.valence_electrons[atname]
        valence_electrons.append(vi)
        V += vi
    # number of bond electrons
    B = OCT-V
    # number of electrons not participating in bonds
    F = V-B
    # constraints
    bonds = []
    # create list of all bonds (I,J)
    i = 0
    for I in range(0, Nat):
        for J in range(I+1, Nat):
            if ConMat[I,J] == 1:
                #
                Zi = atomlist[I][0]
                Zj = atomlist[J][0]
                if is_metal_ion(Zi) or is_metal_ion(Zj):
                    # neglect ionic bonds
                    print "neglect ionic bond between %s%d and %s%d" % (AtomicData.atom_names[Zi-1], I+1, AtomicData.atom_names[Zj-1], J+1)
                    continue
                #
                bonds.append( (I,J) )

                print "BOND %d  = %s%d-%s%d" % (i, AtomicData.atom_names[Zi-1], I+1, AtomicData.atom_names[Zj-1], J+1)
                i += 1
    # number of bonds + number of free electron pairs = number of variables Xk
    Nb = len(bonds)
    N = Nb + Nat
    # EQUALITY CONSTRAINTS
    constraintsD = []
    constraintsB = []
    # 1.  sum_(i,j) X_ij = V
    #     2 * (bonds + lone pairs) = total number of electrons
    d = 2 * np.ones(N)
    b = V
    constraintsD.append(d)
    constraintsB.append(b)
    # 2. total number of bond electrons
    d = np.zeros(N)
    for i in range(0, Nb):
        d[i] = 2    # 2 per bond
    b = B
    constraintsD.append(d)
    constraintsB.append(b)
    # 3. total number of free electrons
    d = np.zeros(N)
    for i in range(Nb,N):
        d[i] = 2    # 2 per lone pair
    b = F
    constraintsD.append(d)
    constraintsB.append(b)
    # 4. octet rule for each atom
    for I in range(0, Nat):
        d = np.zeros(N)
        d[Nb+I] = 2 # 2 electrons per lone pair
        for i,bond in enumerate(bonds):
            if bond[0] == I or bond[1] == I:
                d[i] = 2   # 2 electrons per bond 
        b = ks[I]
        
        constraintsD.append(d)
        constraintsB.append(b)
    # 5. hydrogen and carbon should have no formal chagres, zinc should have formal charge +2
    Nc = 0 # number of hydrogens and carbons
    formal_charge_constraints = {1: 0, 6: 0, 30: +2} # hydrogen and carbon: 0, zinc: +2
    for I in range(0, Nat):
        Zi = atomlist[I][0]
        if Zi in formal_charge_constraints.keys():
            Nc += 1
            # hydrogen or carbon
            d = np.zeros(N)
            b = valence_electrons[I] - formal_charge_constraints[Zi]  # no formal charge FC=valence electrons - electrons attributed to atom I 
            # electrons in bonds attached to atom I
            for i,bond in enumerate(bonds):
                if bond[0] == I or bond[1] == I:
                    d[i] = 1  # bond electrons are split equally between partners
            # electrons in free electron pairs attached to atom I
            d[Nb+I] = 2 # all 2 electrons from a lone pair belong to the atom

            constraintsD.append(d)
            constraintsB.append(b)
    # UPPER BOUND CONSTRAINTS
    ub_constraintsD = []
    ub_constraintsB = []
    # 1. Each bond should have at least a bond order of 1
    #      bo(i) >= 1    <=>  -bo(i) <= -1
    #
    for i,bond in enumerate(bonds):
        d = np.zeros(N)
        d[i] = -1
        b = -1
        ub_constraintsD.append(d)
        ub_constraintsB.append(b)

    # 2. nitrogen/oxygen should have a formal charge of at most -1/-2
    Nfc = 0 # number of atoms which can have formal charges
    fc_lower_limit = {7: -1, 8: -2}  # nitrogen can have at most 1 addition electronal, oxygen at most 2
    for I in range(0, Nat):
        Zi = atomlist[I][0]
        if Zi in fc_lower_limit.keys():
            Nfc += 1
            #
            d = np.zeros(N)
            b = valence_electrons[I]-fc_lower_limit[Zi]  # formal charge FC=-(valence electrons - electrons attributed to atom I) of at most 1 or 2 electrons
            # electrons in bonds attached to atom I
            for i,bond in enumerate(bonds):
                if bond[0] == I or bond[1] == I:
                    d[i] = 1  # bond electrons are split equally between partners
            # electrons in free electron pairs attached to atom I
            d[Nb+I] = 2 # all 2 electrons from a lone pair belong to the atom

            ub_constraintsD.append(d)
            ub_constraintsB.append(b)
    # 3. nitrogen should have at most +1 and oxygen should never have a positive formal charge
    fc_upper_limit = {7: +1, 8: 0}
    for I in range(0, Nat):
        Zi = atomlist[I][0]
        if Zi in fc_upper_limit.keys():
            d = np.zeros(N)
            b = valence_electrons[I]-fc_upper_limit[Zi]  # formal charge FC=-(valence electrons - electrons attributed to atom I) of at least 0
            # electrons in bonds attached to atom I
            for i,bond in enumerate(bonds):
                if bond[0] == I or bond[1] == I:
                    d[i] = 1  # bond electrons are split equally between partners
            # electrons in free electron pairs attached to atom I
            d[Nb+I] = 2 # all 2 electrons from a lone pair belong to the atom

            ub_constraintsD.append(-d)
            ub_constraintsB.append(-b)
    
    
    print "number of atoms: %d" % Nat
    print "number of bonds: %d" % Nb
    print "number of hydrogens+carbons: %s" % Nc
    print "number of equality constraints: %d" % len(constraintsB)
    print "number of upper-bound constraints: %d" % len(ub_constraintsB)
    print "number of valence electrons V=%d" % V
    print "number of bond electrons B=%d" % B
    print "number of variables: %d" % N
    # additional constraints on free electron pairs
    bounds = []

    for i,bond in enumerate(bonds):
        bounds.append( (1,3) ) # at least single at most triple bond
    for (Zi,posi) in atomlist:
        if Zi == 7:
            bounds.append((0,2))  # at most two lone pairs on nitrogen
        elif Zi == 8:
            # oxygen
            bounds.append((0,3)) # at most 3 lone pairs on oxygen
        elif Zi == 30:
            bounds.append((0,5)) # 10 electrons on zinc
        else:
            # no free electron pairs on hydrogen, carbon, etc.
            bounds.append((0,0))

    # integer linear program
    c = np.ones(N)
    # equality constraints
    Nc = len(constraintsB)
    A_eq = np.array(constraintsD)  # number of bonds = number of bonding electron pairs = 2*number of electrons
    b_eq = np.array(constraintsB)
    # upper-bound constraints
    Nc_ub = len(ub_constraintsB)
    A_ub = np.array(ub_constraintsD)
    b_ub = np.array(ub_constraintsB)

    if latte_file != "":
        # create input for Latte to determine number of resonance structures compatible with the constraints
        latte_input(latte_file, A_eq, b_eq, A_ub, b_ub)
        print "Latte input written to %s" % latte_file
        #nr_resonances = latte_count(latte_file)
        #print "number of resonance structures: %d" % nr_resonances
    """
    # valid Lewis structure for nitromethane
    X = np.zeros(N)
    X[0:5] = 1
    X[5] = 2
    X[6+0] = 0
    X[6+1] = 0
    X[6+2] = 2
    X[6+3] = 3
    X[6+4:6+7]= 0
    print "test solution X"
    print X
    print "equality constraints A_eq.X:"
    print (np.dot(A_eq, X) - b_eq)
    print "upper-bound constraints A_ub.X:"
    print (b_ub - np.dot(A_ub, X))

    print "A_eq"
    print A_eq
    print "b_eq"
    print b_eq
    """
    """
    # input for LattE
    print "Input for LattE"
    print "------------------------------------"
    m,d = A_eq.shape
    print "%d %d" % (m,d+1)
    for i in range(0, m):
        print "%2d  " % b_eq[i],
        for j in range(0, d):
            print "%2d " % (-A_eq[i,j]),
        print ""
    print "linearity %d " % m,
    for i in range(0, m):
        print "%d " % (i+1),
    print ""
    print "nonnegative %d " % d,
    for j in range(0, d):
        print "%d " % (j+1),
    print ""
    print "-------------------------------------"
    """
    """
    # 
    X = np.zeros(N)
    for i in range(0, Nb):
        X[i] = 2  # single bonds with two electrons
        
    print "A_eq"
    print A_eq
    print "b_eq"
    print b_eq
    """
    X = solve_linprog(c, A_ub, b_ub, A_eq, b_eq, bounds)
    print "c.X = %s" % np.dot(c, X)
    print "equality constraints:"
    dif = np.dot(A_eq, X) - b_eq
    for n in range(0, Nc):
        print "constraint %d  deviation = %s" % (n, dif[n])
    print "upper bound constraints:"
    dif = np.dot(A_ub, X) - b_ub
    for n in range(0, Nc_ub):
        print "ub constraint %d  deviation = %s <? 0" % (n, dif[n])
    
    bond_orders = X[:Nb]
    print "Bond Orders"
    print "==========="
    for i,bond in enumerate(bonds):
        bo = bond_orders[i]
        I,J = bond
        #
        Zi = atomlist[I][0]
        Zj = atomlist[J][0]

        print "BOND %s%d-%s%d   BO=%f" % (AtomicData.atom_names[Zi-1], I+1, AtomicData.atom_names[Zj-1], J+1, bo)
    lone_pairs = X[Nb:Nb+Nat]
    print "Free Electron Pairs"
    print "==================="
    for I in range(0, Nat):
        fI = lone_pairs[I]
        if (fI > 0):
            Zi = atomlist[I][0]
            atname = AtomicData.atom_names[Zi-1]
            print "Atom %s%d  number of lone pairs = %d" % (atname, I+1,fI)
    formal_charges = []
    print "Formal Charges"
    print "=============="
    for I in range(0, Nat):
        fcI = 0.0
        for i,bond in enumerate(bonds):
            # bonding electrons are shared
            if bond[0] == I or bond[1] == I:
                fcI += bond_orders[i]
        fcI += 2.0*lone_pairs[I]
        fcI -= valence_electrons[I]
        if abs(fcI) > 0.0:
            Zi = atomlist[I][0]
            atname = AtomicData.atom_names[Zi-1]
            # change sign of fcI so that having less electrons leads to a positive number
            print "Atom %s%d  formal charge = %f" % (atname, I+1, -fcI)
        formal_charges.append( -fcI )
            
    return bonds, bond_orders, lone_pairs, formal_charges

def latte_input(latte_file, A_eq, b_eq, A_ub, b_ub):
    """
    write input file in half-space format for Latte
    """
    fh = open(latte_file, "w")
    m,d = A_eq.shape
    n,d2 = A_ub.shape
    assert d == d2
    # equality constraints
    print>>fh, "%d %d" % (m+n,d+1)
    for i in range(0, m):
        print>>fh, "%2d  " % b_eq[i],
        for j in range(0, d):
            print>>fh, "%2d " % (-A_eq[i,j]),
        print>>fh, ""
    # upper bound constraints
    for i in range(0, n):
        print>>fh, "%2d " % b_ub[i],
        for j in range(0, d):
            print>>fh, "%2d " % (-A_ub[i,j]),
        print>>fh, ""
    # Which constraints should be treated as equality constraints?
    print>>fh, "linearity %d " % m,
    for i in range(0, m):
        print>>fh, "%d " % (i+1),
    print>>fh, ""
    print>>fh, "nonnegative %d " % d,
    for j in range(0, d):
        print>>fh, "%d " % (j+1),
    print>>fh, ""
    fh.close()

def latte_count(latte_file):
    os.system("count %s &> latte.out" % latte_file)
    fh = open("latte.out", "r")
    lines = fh.readlines()
    print "LINES = %s" % lines
    fh.close()
    for l in lines:
        print "LINE=%s" % line.strip()
        if ("The polytope has" in l) and ("vertices." in l):
            nr_vertices = int(l.split()[3])
        elif "The number of lattice points is:" in l:
            print "##################################################"
            exit(-1)
            if len(l.split()) == 9:
                nr_inside = int(l.split()[7])
            else:
                nr_inside = 0
            break
    else:
        raise Exception("Latte calculation failed! See output")
    print "DONE"
    return nr_inside+nr_vertices
    
if __name__ == "__main__":
    import sys
    
    xyz_file = sys.argv[1]
    atomlist = XYZ.read_xyz(xyz_file)[0]
    ConMat = XYZ.connectivity_matrix(atomlist)
    assign_bond_orders(atomlist, ConMat, latte_file=xyz_file.replace(".xyz", ".hrep.latte"))
    
