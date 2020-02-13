"""
Pipek-Mezey orbital localization

A unitary transformation of the MO coefficients is found by maximizing a localization 
metric, which is based on Mulliken charges, so that molecular orbitals are localized on fragments. 

References
----------
[1] J. Pipek, P. Mezey, 
    "A fast intrinsic localization procedure applicable for ab initio and semiempirical
     linear combination of atomic orbital wave functions",
     J. Chem. Phys. 90, 4916 (1989)
[2] I-M.Hoyvik, B. Jansik, P. Jorgensen
    "Pipek-Mezey Localization of Occupied and Virtual Orbitals"
     J. Comp. Chem. 34,  1456-1462 (2013)
"""
import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
from scipy import optimize

from DFTB.Analyse import MolecularGraph as MG

def orbital_rotation(x, n):
    """
    create unitary n x n matrix U from parameter vector x, which contains
    unique elements of a skew-symmetric n x n matrix u, such that

      U = exp(u)

    """
    # The Lie algebra corresponding to the special orthogonal group SO(n)
    # consists of the skew-symmetric n x n matrices. An element of the
    # Lie algebra is parametrized by n*(n-1)/2 real numbers, which
    # are stored in the parameter vector x.
    assert len(x) == (n*(n-1))/2

    # construct skew-symmetric n x n matrix. The parameter vector x contains
    # the unique matrix elements above the diagonal.
    u = np.zeros((n,n))
    ij = 0
    for i in range(0, n):
        for j in range(i+1,n):
            # antisymmetric matrix
            u[i,j] =  x[ij]
            u[j,i] = -x[ij]
            ij += 1

    # construct element of Lie group U = exp(u)
    U = sla.expm(u)
    
    return U
    

def localize_pipek_mezey(atomlist, orbs, orbe, focc, S, valorbs):
    """
    localize orbitals onto diconnected fragments using the Pipek-Mezey method.

    Parameters
    ----------
    atomlist   : list of tuples (Zat, [xi,yi,zi]) with molecular geometry
    orbs       : MO coefficients, orbs[:,i] are the coefficients of orbital i
    orbe       : MO energies, orbe[i] is the Kohn-Sham energy of the i-th orbital
    focc       : occupation numbers, if focc[i] > 0, the i-th orbital is occupied
    S          : overlap matrix in AO basis
    valorbs    : dictionary with list of quantum numbers of valence
                 orbitals for each atom

    Returns
    -------
    orbs_loc   : MO coefficients of localized orbitals, the sets of
                 occupied and virtual orbitals are localized independently
    orbe_loc   : expectation values of the energy for the localized orbitals,
                 orbe_loc[i] = <i|H|i>. The Kohn-Sham Hamiltonian is not
                 diagonal in the basis of localized orbitals, but the total
                 energy of the ground state Slater determinant is not changed
                 by the unitary transformation that localizes the occupied
                 orbitals.
    U          : numpy array of shape (nmo,nmo) with unitary transformation
                 that localizes both occupied and virtual orbitals
    frags      : Each localized orbital is assigned to a fragment, frags[i]
                 gives the index (0-based) of the fragment to which orbital
                 i belongs
    """
    # number of atoms
    nat = len(atomlist)
    # identify disconnected fragments
    fragment_graphs = MG.atomlist2graph(atomlist, hydrogen_bonds=False)
    # number of fragments
    nfrag = len(fragment_graphs)
    # list of indices into atomlist of atoms belonging to each fragment
    fragment_indices = [MG.graph2indeces(g) for g in fragment_graphs]
    # mapping from atom indices to fragment index
    atom2frag = np.zeros(nat, dtype=int)
    for ifrag, atom_indices in enumerate(fragment_indices):
        atom2frag[atom_indices] = ifrag

    # Show which atom indices belong to which fragment
    print ""
    print "  Fragments for Orbital Localization"
    print "  =================================="
    print "  For each disconnected fragment the indices (1-based) "
    print "  of the atoms belonging to that fragment are listed."
    print ""
    for ifrag, atom_indices in enumerate(fragment_indices):
        print "  Fragment %2.0d :   " % (ifrag+1),
        for i,iatom in enumerate(atom_indices):
            print "%d " % (iatom+1),
            if (i+1) % 10 == 0:
                print ""
                print "                   ",
        print ""

        
    # occupied and virtual orbitals are localized separately
    # MO coefficients
    orbs_occ = orbs[:,focc > 0]
    orbs_virt = orbs[:, focc == 0.0]
    # MO eigenenergies
    orbe_occ = orbe[focc > 0]
    orbe_virt = orbe[focc == 0.0]
    
    # localize occupied orbitals
    Uocc, orbs_occ_loc, Q_occ_loc = _localize_orbital_set(atomlist, orbs_occ, S, valorbs, atom2frag, nfrag)
    nao,nocc = orbs_occ_loc.shape
    # localize virtual orbitals
    Uvirt, orbs_virt_loc, Q_virt_loc = _localize_orbital_set(atomlist, orbs_virt, S, valorbs, atom2frag, nfrag)
    nao,nvirt = orbs_virt_loc.shape
    
    orbs_loc = np.zeros(orbs.shape)
    orbs_loc[:,focc > 0] = orbs_occ_loc
    orbs_loc[:,focc == 0.0] = orbs_virt_loc

    # expectation values of energy for localized orbitals
    #  ~   ~                 *                *
    # <i|H|i> = sum <k|H|j> U   U   = sum e  U   U
    #           k,j          ki  ji    j   j  ji  ji
    orbe_occ_loc = np.zeros(orbe_occ.shape)
    orbe_virt_loc = np.zeros(orbe_virt.shape)
    # Each localized orbital can be assign to one fragment, the one for
    # which Q_ii^F is largest
    frags_occ = np.zeros(nocc, dtype=int)
    frags_virt = np.zeros(nvirt, dtype=int)
    for i in range(0, nocc):
        orbe_occ_loc[i] = np.sum(Uocc[:,i].conjugate() * orbe_occ * Uocc[:,i])
        frags_occ[i] = np.argmax(Q_occ_loc[i,i,:])
    for a in range(0, nvirt):
        orbe_virt_loc[a] = np.sum(Uvirt[:,a].conjugate() * orbe_virt * Uvirt[:,a])
        frags_virt[a] = np.argmax(Q_virt_loc[a,a,:])
    # combine occupied and virtual orbitals again
    orbe_loc = np.zeros(orbe.shape)
    orbe_loc[focc > 0] = orbe_occ_loc
    orbe_loc[focc == 0.0] = orbe_virt_loc
    frags = np.zeros(nocc+nvirt, dtype=int)
    frags[focc > 0] = frags_occ
    frags[focc == 0.0] = frags_virt

    # combine unitary transformations for occupied and virtual orbitals
    #      ( Uocc  0     )
    #  U = (             )
    #      (  0    Uvirt )
    nmo = nocc+nvirt
    Uloc = np.zeros((nmo, nmo))
    Uloc[np.ix_(focc > 0   ,focc > 0)]    = Uocc
    Uloc[np.ix_(focc == 0.0,focc == 0.0)] = Uvirt

    print ""
    print " Pipek-Mezey localization of orbitals"
    print " ------------------------------------"
    print "     MO(F)        Fragment Charges Q_{i,i}^F          "
    print "            i   ",
    for ifrag in range(0, nfrag):
        print "      F=%2.1d    " % (ifrag+1),
    print ""
    for imo in range(0, nocc):
        ifrag = np.argmax(Q_occ_loc[imo,imo,:])
        print " occ(%d)   %3.1d" % (ifrag+1,imo+1),
        for ifrag in range(0, nfrag):
            print "        %+4.3f" % Q_occ_loc[imo,imo,ifrag],
        print ""
    for imo in range(0, nvirt):
        ifrag = np.argmax(Q_virt_loc[imo,imo,:])
        print " virt(%d)  %3.1d" % (ifrag+1,nocc+imo+1),
        for ifrag in range(0, nfrag):
            print "        %+4.3f" % Q_virt_loc[imo,imo,ifrag],
        print ""
    print ""

    return orbs_loc, orbe_loc, Uloc, frags
    
def _localize_orbital_set(atomlist, C, S, valorbs, atom2frag, nfrag):
    """
    localize a set of orbitals (occupied or virtual) 
    given as the nao x nmax matrix C to fragments 

    Returns
    -------
    Uopt     :  unitary transformation that localizes the orbitals
    Copt     :  Uopt.C, localized orbitals, matrix of shape (nao,nmo)
    Qopt     :  Mulliken MO charges for localized orbitals,
                numpy array of shape (nmo,nmo,nfrag)
    """
    nao,nmo = C.shape
    # Mulliken fragment charges
    #  F         AO        AO
    # Q0    = sum       sum   2 C    S    C
    #  i,j      m on F     n     m,i  m,n  n,j
    Q0 = np.zeros((nmo,nmo,nfrag))

    SC = np.dot(S,C)

    # iterate over atoms
    mu = 0
    for iat,(Zi,posi) in enumerate(atomlist):
        # atom i belongs to fragment f
        f = atom2frag[iat]
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbs[Zi]:
            # add contribution from orbital mu on fragment f
            Q0[:,:,f] += 2*np.outer(C[mu,:], SC[mu,:])
            
            # increase orbital counter
            mu += 1

    # sum of Mulliken charges should equal number of electrons
    #print np.sum(Q0[:,:,0])
    #print np.sum(Q0[:,:,1])

    # objective function that should be maximized w/r/t the orbital
    # rotations U
    #
    #             MO         F     2
    #   L(U) = sum    sum  (Q(U)  )  
    #             i=1    F   ii
    # with
    #    F                     F
    #   Q (U) = sum     U    Q0    U
    #    ii        j,k   j,i   j,k  k,i
    #
    # The orthogonal transformation U is parametrized by the corresponding
    # element of the Lie algebra, the antisymmetric matrix u, which
    # has n*(n-1)/2 parameters, (n == number of MOs).

    # number of parameters of orbital rotation
    npar = (nmo*(nmo-1))/2

    # Mulliken MO charges on fragments Q^F_{i,j}
    def mullikenQ(x):
        # orthogonal n x n transformation
        U = orbital_rotation(x, nmo)
        Ut = U.transpose()
        # transformed Q = U^T . Q0 . U
        Q = np.zeros((nmo,nmo,nfrag))
        for f in range(0, nfrag):
            Q[:,:,f] = np.dot(Ut, np.dot(Q0[:,:,f], U))
        return Q
    
    # loc_metric(x) = L(U(x)) as a function of uni
    def funcL(x):
        Q = mullikenQ(x)
        # L(U)
        L = 0.0
        for f in range(0, nfrag):
            for i in range(0, nmo):
                L += Q[i,i,f]**2
        # Maximizing L is equivalent to minimizing -L
        return -L

    def gradL(x):
        # gradient of objective function L
        Q = mullikenQ(x)
        # gradient of L(U) w/r/t parameter vector x
        gradL = np.zeros(npar)
        ij = 0
        for i in range(0, nmo):
            for j in range(i+1,nmo):
                gradL[ij] = 2 * np.sum( (Q0[j,j,:] - Q0[i,i,:]) * (Q[i,j,:] + Q[j,i,:]) )
                ij += 1
        # Maximizing L is equivalent to minimizing -L        
        return -gradL

    def hessL(x):
        # Hessian of objective function
        # not implemented yet, needed for 'trust-krylov'
        pass
        
    # The parameter vector x0 = 0 corresponds to U=Id. This is the initial
    # guess for the optimization. 
    x0 = np.zeros(npar)
    # MO Mulliken charges before localization
    Q0 = mullikenQ(x0)
    
    #optimize.minimize(funcL, x0, jac=gradL, hess=hessL, method="trust-krylov")
    res = optimize.minimize(funcL, x0, jac=gradL, method="L-BFGS-B", options={'disp': False, 'gtol': 1.0e-7})

    xopt = res.x
    
    # MO Mulliken charges
    Qopt = mullikenQ(xopt)
    
    # unitary transformation that localized orbitals
    Uopt = orbital_rotation(xopt, nmo)
    # transform MO coefficients of orbital set
    Copt = np.dot(C, Uopt)

    return Uopt, Copt, Qopt

class OrbitalLocalization:
    def __init__(self, tddftb):
        self.tddftb = tddftb
    def projectLocalizedExcitation(self, ifrag,iorb, afrag,aorb, threshold=0.01):
        """
        A localized excited state is represented by a single excitation from an 
        occupied orbital localized on one fragment to a virtual orbital localized 
        on the same or another fragment. If this CIS state is projected onto the 
        basis of adiabatic eigenstates, a superposition of states results. 

        Parameters
        ----------
        ifrag      :  index of fragment for occupied orbital (1-based)
        iorb       :  occupied orbital counted from the HOMO downward, so iorb=0 corresponds to
                      the HOMO on fragment `ifrag`, iorb=1 corresponds to HOMO-1, etc.
                      Only the occupied orbitals on fragment `ifrag` are considered.
        afrag      :  index of fragment for virtual orbital (1-based)
        aorb       :  virtual orbital counted from the LUMO upward, so aorb=0 corresponds to
                      the LUMO, aorb=1 to LUMO+1 etc. 
                      Only the virtual orbitals on fragment `afrag` are considered.

        Returns
        -------
        proj       :  coefficients of the localized excitation in the basis of adiabatic
                      eigenstates, proj[0] corresponds to the ground state and is always 0.
        """
        dftb = self.tddftb.dftb2
        # localize orbitals and assign them to fragments
        orbs_loc, orbe_loc, U, frags = localize_pipek_mezey(dftb.atomlist, dftb.orbs, dftb.orbe,
                                                            dftb.f, dftb.S, dftb.valorbs)        
        
        # overlap between localized and canonical MOs equals the unitary transformation matrix
        #    ~
        # <i|j> = <i| sum |k> U    = sum <i|k> U    = U
        #              k       k,j    k         k,j    i,j
        Smo = U
        #
        
        # indeces of occupied spatial orbitals
        Nelec = self.tddftb.dftb2.Nelec_val
        assert Nelec % 2 == 0
        # ground state occupation
        occs_0 = np.array([i for i in range(0, int(Nelec/2))])
        # indices active occupied and virtual orbitals 
        occ, virt = self.tddftb.getActiveOrbitals()
        nocc = len(occ)
        nvirt = len(virt)

        # overlap between localized excitation (ifrag,iorb) --> (afrag, aorb)
        # of the i-th highest occupied orbital on fragment ifrag
        # to the a-th lowest unoccupied orbital on fragment afrag
        occs_ia = np.copy(occs_0)
        # indices of orbitals that are occupied and belong to fragment
        assert iorb >= 0
        assert aorb >= 0
        occ_i = occ[frags[occ] == ifrag-1][-iorb-1]
        virt_a = virt[frags[virt] == afrag-1][aorb]
        #i = np.where((frags == ifrag-1) & (dftb.f > 0.0))[0][-iorb-1]
        #a = np.where((frags == afrag-1) & (dftb.f == 0.0))[0][aorb]

        print "excitation from HOMO-%d on fragment %d  to LUMO+%d on fragment %d" % (iorb,ifrag,aorb,afrag)
        print "occupied orbital corresponds to MO %d" % (occ_i+1)
        print "virtual orbital corresponds to MO %d" % (virt_a+1)
        
        # coefficients of adiabatic eigenstates
        C = self.tddftb.Cij
        # number of adiabatic eigenstates (depends on size of active space)
        Nst = C.shape[0]

        # projection of localized excited state onto adiabatic eigenstates
        # The projection on ground state, proj[0], is always zero.

        proj = np.zeros(Nst+1, dtype=complex)
        # Adiabatic eigenstates are superposition of single-excitations j->b
        # with singlet spin pairing.
        # The projection of the prepared state |Psi_ia> onto the adiabatic basis
        # is given by
        #                (I)*      1      1
        #  c  = sum     C     < Psi  | Psi  >
        #   I      j,b   j,b       jb     ia

        # The formulae for overlaps between singlet configuration state functions
        # can be found in chapter 9.1. 'Scalar non-adiabatic couplings' in
        #    A. Humeniuk, PhD thesis (2018),
        #    "Methods for Simulating Light-Induced Dynamics in Large Molecular Systems"
        
        #   1      1
        #  < Psi  | Psi  > =  <1,...,b,...,N/2|1,...,a,...,N/2> <1,...,j,...,N/2|1,...,i,...,N/2>
        #       jb     ia
        #                   + <1,...,j,...,N/2|1,...,a,...,N/2> <1,...,b,...,N/2|1,...,i,...,N/2>
        #
        #                  =  det_ba * det_ji + det_ja * det_bi
        #
        # where i,j and a,b are occupied and virtual spatial orbitals, respectively and <...|...> are
        # overlaps of Slater determinants. 
        
        # <1,...,j,...|1,...,i,...>
        S_ji = Smo[occs_0,:][:,occs_0]
        det_ji = la.det(S_ji)
        
        # occupation numbers of single excitation
        occs_ia[occ_i] = virt_a
        # <1,...,j,...|1,...,a,...>
        S_ja = Smo[occs_0,:][:,occs_ia]
        det_ja = la.det(S_ja)
        
        for j in range(0, nocc):
            for b in range(0, nvirt):
                if abs(C[:,j,b]).max() < threshold:
                    # skip the j->b excitation
                    continue
                # occupied orbitals in the configuration state function
                # |Psi_jb>
                occs_jb = np.copy(occs_0)
                occs_jb[occ[j]] = virt[b]
                # select part of overlap matrix for orbitals
                # in |Psi_ia> and |Psi_jb>
                # <1,...,b,...|1,...,a,...>
                S_ba = Smo[occs_jb,:][:,occs_ia]
                det_ba = la.det(S_ba)
                # <1,...,b,...|1,...,i,...>
                S_bi = Smo[occs_jb,:][:,occs_0]
                det_bi = la.det(S_bi)

                # loop over excited states
                for I in range(1, Nst+1):
                    # coefficient c_I
                    # see eqn. (9.39) in thesis
                    proj[I] += C[I-1,j,b].conjugate() * (det_ba * det_ji + det_ja * det_bi)

        print ""
        print " Projection onto adiabatic eigenstates "
        print " ------------------------------------- "
        print "   adiab.state     Coefficients      "
        for I in range(0,Nst+1):
            print "      %3.1d           %+5.4f" % (I, proj[I].real)
        print "   Norm          %e" % sla.norm(proj)
        
        return proj
