"""
assigns partial charges and partial dipoles to
individual atoms based on Mulliken's partitioning
"""

# WARNING: this code is slow, the fortran implementation extensions.mulliken should be used instead

import numpy as np

def monopoles(atomlist, P, P0, S, valorbs):
    """
    Parameters:
    ============
    atomlist: list of tuples (Zi,(xi,yi,zi)) for each atom of type Zi
    P: density matrix, shape (Norb,Norb)
    P0: reference density matrix, shape (Norb,Norb)
    S: overlap matrix, shape (Norb,Norb) 
    valorbs: dictionary with list of valence orbitals 
             with quantum numbers (ni,li,mi)

    Returns:
    ========
    q: Mulliken charges for each atom, shape (Nat)
    dq: partial Mulliken charges dq=q-q0 for each atom, shape (Nat)

    The monopoles are computed as

       q_A = sum_(mu in A) sum_nu (P-P0)_(mu,nu) S_(mu,nu)
    """
    Nat = len(atomlist)
    q = np.zeros(Nat)
    dq = np.zeros(Nat)
    #
    dP = P-P0
    # iterate over atoms
    mu = 0
    for i,(Zi,posi) in enumerate(atomlist):
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbs[Zi]:
            # iterate over atoms
            nu = 0
            for j,(Zj,posj) in enumerate(atomlist):
                # iterate over orbitals on center j
                for (nj,lj,mj) in valorbs[Zj]:
                    q[i] += P[mu,nu]*S[mu,nu]
                    dq[i] += dP[mu,nu]*S[mu,nu]
                    nu += 1
            mu += 1
    return q, dq

def dipoles(atomlist, P, P0, S, D, valorbs, Rcc):
    """
    Parameters:
    ============
    atomlist: list of tuples (Zi,(xi,yi,zi)) for each atom of type Zi
    P: density matrix, shape (Norb,Norb)
    P0: reference density matrix, shape (Norb,Norb)
    S: overlap matrix, shape (Norb,Norb)
    D: dipole matrix elements in the atomic basis, shape (Norb,Norb,3) 
    valorbs: list of valence orbitals with quantum numbers (ni,li,mi)
    Rcc: center of nuclear charge, Rcc = (sum_i Zi*Ri)/(sum_i Zi)

    Returns:
    ========
    dip: Mulliken dipoles for each atom, shape (Nat,3)
    ddip: partial Mulliken dipoles ddip=dip-dip0 for each atom, shape (Nat,3)

    The monopoles are computed as

       dip_A = sum_(mu in A) sum_nu (P-P0)_(mu,nu) D_(mu,nu)
    """
    Nat = len(atomlist)
    dip = np.zeros((Nat,3))
    ddip = np.zeros((Nat,3))
    #
    dP = P-P0
    # iterate over atoms
    mu = 0
    for i,(Zi,posi) in enumerate(atomlist):
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbs[Zi]:
            # iterate over atoms
            nu = 0
            for j,(Zj,posj) in enumerate(atomlist):
                # iterate over orbitals on center j
                for (nj,lj,mj) in valorbs[Zj]:
                    dip[i,:] += P[mu,nu]*D[mu,nu,:]
                    ddip[i,:] += dP[mu,nu]*D[mu,nu,:]
                    # subtract dipole moment of point charges
                    dip[i,:] -= P[mu,nu]*S[mu,nu] * Rcc
                    ddip[i,:] -= dP[mu,nu]*S[mu,nu] * Rcc

                    nu += 1
            mu += 1
            
    return dip, ddip

def multipoles(atomlist, P, P0, S, D, valorbs, Rcc):
    """
    Parameters:
    ============
    atomlist: list of tuples (Zi,(xi,yi,zi)) for each atom of type Zi
    P: density matrix, shape (Norb,Norb)
    P0: reference density matrix, shape (Norb,Norb)
    S: overlap matrix, shape (Norb,Norb) 
    D: dipole matrix elements in the atomic basis, shape (Norb,Norb,3) 
    valorbs: list of valence orbitals with quantum numbers (ni,li,mi)
    Rcc: center of nuclear charge, Rcc = (sum_i Zi*Ri)/(sum_i Zi)

    Returns:
    ========
    qm: Mulliken multipoles (mono- & dipoles) for each atom, shape (4*Nat)
    dqm: partial Mulliken multipoles for each atom, shape (4*Nat)
    """
    Nat = len(atomlist)
    qm = np.zeros(4*Nat)
    # qm[4*i]  is the monopole of the i-th atom
    # qm[4*i+1:4*(i+1)] is the dipole vector on the i-th atom
    dqm = np.zeros(4*Nat)
    #
    dP = P-P0
    # iterate over atoms
    mu = 0
    for i,(Zi,posi) in enumerate(atomlist):
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbs[Zi]:
            # iterate over atoms
            nu = 0
            for j,(Zj,posj) in enumerate(atomlist):
                # iterate over orbitals on center j
                for (nj,lj,mj) in valorbs[Zj]:
                    # full charges and dipoles
                    qm[4*i] += P[mu,nu]*S[mu,nu]
                    qm[4*i+1:4*(i+1)] += P[mu,nu]*D[mu,nu,:]
                    # dipole moments have to be calculated relative to the center of charge
                    qm[4*i+1:4*(i+1)] -= P[mu,nu]*S[mu,nu] * Rcc
                    # partial 
                    dqm[4*i] += dP[mu,nu]*S[mu,nu]
                    dqm[4*i+1:4*(i+1)] += dP[mu,nu]*D[mu,nu,:]
                    # partial dipole moments have to be calculated
                    # relative to the center of charge
                    dqm[4*i+1:4*(i+1)] -= dP[mu,nu]*S[mu,nu] * Rcc

                    nu += 1
            mu += 1

    return qm, dqm

from DFTB import AtomicData

def save_partial_dipoles(filename, atomlist, ddip):
    """
    write the molecular geometry (in bohr) followed by the 
    dipole vectors to a file. This file can be visualized
    with blender
    """
    # write partial dipoles
    fh = open(filename, "w")
    nat = len(atomlist)
    print>>fh, "%s" % nat
    for (Zi,(x,y,z)) in atomlist:
        print>>fh, "%s   %s %s %s" % (AtomicData.atom_names[Zi-1], x, y, z)
    for i in range(0, nat):
        print>>fh, "%s %s %s" % tuple(ddip[i,:])
    fh.close()
    print "wrote partial dipoles to file %s" % filename

