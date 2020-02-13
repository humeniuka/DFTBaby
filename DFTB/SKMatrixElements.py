"""
computes matrix elements of the reference Hamiltonian H0,
the overlap matrix and the dipoles D between the valence orbitals
using Slater Koster rules
"""
import numpy as np

# There a two codes for calculating matrix elements using Slater-Koster rules:
#   - a more or less well tested python implementation ('PY')
#   - and a much fast OpenMP-parallelized Fortran implementation ('F90')
# This flag allows to switch between the two
SLAKO_IMPLEMENTATION = "F90" #"PY"  # "F90"

def H0andS(atomlist, valorbs, SKT, orbital_energies, Mproximity):
    """
    # check that F90 and Python implementation give the same results
    # parallelized fortran code
    import time
    ta = time.time()
    S_f90, H0_f90 = _H0andS_f90(atomlist, atomlist, valorbs, SKT, orbital_energies, Mproximity)
    tb = time.time()
    time_f90 = tb-ta
    # serial code
    ta = time.time()
    S_py, H0_py = _H0andS(atomlist, atomlist, valorbs, SKT, orbital_energies, Mproximity)
    tb = time.time()
    time_py = tb-ta
    # compare
    from DFTB import utils
    orb_labels = ["%d" % i for i in range(0, S_f90.shape[0])]
#    print "S F90"
#    print utils.annotated_matrix(S_f90, orb_labels, orb_labels)
#    print "S py"
#    print utils.annotated_matrix(S_py, orb_labels, orb_labels)
    errS = (abs(S_py-S_f90)).max()
    assert errS < 1.0e-3, "python and Fortran implementations disagree for S, error = %s" % errS
#    print "H F90"
#    print utils.annotated_matrix(H0_f90, orb_labels, orb_labels)
#    print "H py"
#    print utils.annotated_matrix(H0_py, orb_labels, orb_labels)
    errH0 = (abs(H0_py-H0_f90)).max()
    assert errH0 < 1.0e-3, "python and Fortran implementations disagree for H0, error = %s" % errH0

    print "Timing:"
    print "Fortran: %s seconds" % time_f90
    print "Python : %s seconds" % time_py
    exit(-1)
    """
    if SLAKO_IMPLEMENTATION == "F90":
        # fortran implementation
        S, H0 = _H0andS_f90(atomlist, atomlist, valorbs, SKT, orbital_energies, Mproximity)
    else:
        # python implementation
        S, H0 = _H0andS(atomlist, atomlist, valorbs, SKT, orbital_energies, Mproximity)

    return S, H0
    
def count_orbitals(atomlist, valorbs):
    Norb = 0
    for i,(Zi,posi) in enumerate(atomlist):
        Norb += len(valorbs[Zi])
    return Norb

def _H0andS(atomlist1, atomlist2, valorbs, SKT, orbital_energies, Mproximity):
    """
    construct the matrix elements of H0 and the overlaps S for a fixed
    nuclear geometry using Slater-Koster transformations. 
    This has to be performed only once, since during the DFTB iterations 
    only the mo coefficients and the partial charges change.

    Parameters:
    ===========
    atomlist1, atomlist2: list of tuples (Zi,[xi,yi,zi]) of atom types and positions
    valorbs: list of valence orbitals with quantum numbers (ni,li,mi)
    SKT: Slater Koster table
    orbitals_energies: dictionary with orbital energies for each atom type
    Mproximity: M[i,j] == 1, if the atoms i and j are close enough 
      so that the overlap between orbitals on i and j should be 
      computed
    
    Returns:
    ========
    S: overlap matrix
    H0: reference KS hamiltonian
    """
    # count valence orbitals
    Norb1 = count_orbitals(atomlist1, valorbs)
    Norb2 = count_orbitals(atomlist2, valorbs)

    H0 = np.zeros((Norb1,Norb2))
    S = np.zeros((Norb1,Norb2))
    # iterate over atoms
    mu = 0
    for i,(Zi,posi) in enumerate(atomlist1):
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbs[Zi]:
            # iterate over atoms
            nu = 0
            for j,(Zj,posj) in enumerate(atomlist2):
                # iterate over orbitals on center j
                for (nj,lj,mj) in valorbs[Zj]:
                    if Mproximity[i,j] == 1:
                        if mu < nu:
                            if Zi <= Zj:
                                # the first atom given to getHamiltonian() or getOverlap()
                                # has to be always the one with lower atomic number
                                if i == j:
                                    # different orbitals on the same atom should be orthogonal
                                    assert mu != nu
                                    s = 0.0
                                    h0 = 0.0
                                else:
                                    s  = SKT[(Zi,Zj)].getOverlap(li,mi,posi, lj,mj,posj)
                                    h0 = SKT[(Zi,Zj)].getHamiltonian0(li,mi,posi, lj,mj,posj)
                            else:
                                # swap atoms if Zj > Zi
                                s  = SKT[(Zj,Zi)].getOverlap(lj,mj,posj, li,mi,posi)
                                h0 = SKT[(Zj,Zi)].getHamiltonian0(lj,mj,posj, li,mi,posi)
                            H0[mu,nu] = h0
                            S[mu,nu] = s
                        elif mu == nu:
                            assert Zi == Zj
                            H0[mu,nu] = orbital_energies[Zi][(ni,li)] # use the true single particle orbitals energies
                            S[mu,nu] = 1.0 # orbitals are normalized to 1
                        else:
                            # S and H0 are hermitian/symmetric
                            H0[mu,nu] = H0[nu,mu]
                            S[mu,nu] = S[nu,mu]
                    nu += 1
            mu += 1
    return S, H0



# FASTER SLATER-KOSTER RULES WITH FORTRAN
from DFTB import XYZ
from DFTB.extensions import slako
from DFTB.SlaterKoster.SKIntegrals import combine_slako_tables_f90

def atomlist2orbitals(atomlist, valorbs, atom_type_dic):
    """
    Parameters:
    ===========
    atomlist: 

    atom_types_dic: translates atomic numbers into atom types
    """
    # quantum numbers l and m for all orbitals
    ls = []
    ms = []
    # atom_indeces[mu] is the index of the atom which the mu-th orbital sits on
    atom_indeces = []
    # atom_types[mu] is the atom type of the atom which the mu-th orbital sits on
    atom_types = []
    # iterate over atoms
    for i,(Zi,posi) in enumerate(atomlist):
        attype = atom_type_dic[Zi]
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbs[Zi]:
            ls.append(li)
            ms.append(mi)
            atom_indeces.append(i)
            atom_types.append(attype)

    return atom_indeces,atom_types,ls,ms

def atomlist2orbital_energies(atomlist, valorbs, orbital_energies_dic):
    # diagonal energies of the orbitals
    orbital_energies = []
    # iterate over atoms
    for i,(Zi,posi) in enumerate(atomlist):
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbs[Zi]:
            orbital_energies.append( orbital_energies_dic[Zi][(ni,li)] )
    return orbital_energies

def _H0andS_f90(atomlist1, atomlist2, valorbs, SKT, orbital_energies_dic, Mproximity):
    """
    This function calls an external Fortran function that computes the matrix elements.
    Since hashes and tuples cannot be used easily in Fortran, the arguments are brought
    into a form that can be fed into the Fortran function.
    """
    atom_type_dic, spline_deg, tab_filled_SH0, tab_filled_D, \
        (S_knots, S_coefs, H_knots, H_coefs, D_knots, D_coefs) = \
                            combine_slako_tables_f90(SKT)
    
    # number of atoms
    Nat1 = len(atomlist1)
    Nat2 = len(atomlist2)
    # orbitals
    atom_indeces1,atom_types1,ls1,ms1 = atomlist2orbitals(atomlist1, valorbs, atom_type_dic)
    atom_indeces2,atom_types2,ls2,ms2 = atomlist2orbitals(atomlist2, valorbs, atom_type_dic)
    orben = atomlist2orbital_energies(atomlist1, valorbs, orbital_energies_dic)
    assert atom_indeces1 == atom_indeces2, "atomlist1 and atomlist2 should contain the same atoms in the same order!"
    # count valence orbitals
    Norb1 = len(ls1)
    Norb2 = len(ls2)
    # distances and direction
    pos1 = XYZ.atomlist2vector(atomlist1)
    pos1 = np.reshape(pos1,(Nat1,3)).transpose()  # pos(:,i) is the 3d position of atom i
    pos2 = XYZ.atomlist2vector(atomlist2)
    pos2 = np.reshape(pos2,(Nat2,3)).transpose() 
    r,x,y,z = slako.slako.directional_cosines(pos1,pos2)
    # call Fortran function
    S,H0 = slako.slako.h0ands(atom_indeces1,atom_types1,ls1,ms1,
                                           atom_indeces2,atom_types2,ls2,ms2,
                                           r,x,y,z,
                                           Mproximity,
                                           orben,
                                           spline_deg, tab_filled_SH0,
                                           S_knots, S_coefs, H_knots, H_coefs)
    return S,H0
    
def DipoleMatrix(atomlist, valorbs, SKT, Mproximity, S):
    """
    only the Fortran code can use the overlap matrix S, in the python code
    the matrix elements for S are recalculated
    """
    """
    # compare python and fortran code
    # Dipole matrices will disagree slightly, error ~ 1.0e-3, since the interpolation
    # algorithm implemented in Fortran and in numpy are not exactly the same.
    import time
    ta = time.time()
    Dipole_f90 = _DipoleMatrix_f90(atomlist, atomlist, valorbs, SKT, Mproximity, S)
    tb = time.time()
    time_f90 = tb-ta
    # serial code iny python
    ta = time.time()
    Dipole_py = _DipoleMatrix(atomlist, atomlist, valorbs, SKT, Mproximity)
    tb = time.time()
    time_py = tb-ta
    # compare
    from DFTB import utils
    orb_labels = ["%d" % i for i in range(0, Dipole_f90.shape[0])]
#    xyz = ["X","Y","Z"]
#    for i in range(0, 3):
#        print "D F90 %s" % xyz[i]
#        print utils.annotated_matrix(Dipole_f90[:,:,i], orb_labels, orb_labels)
#        print "D py  %s" % xyz[i]
#        print utils.annotated_matrix(Dipole_py[:,:,i], orb_labels, orb_labels)
    # show differences
    absD = abs(Dipole_py-Dipole_f90)
    print "Differing matrix elements"
    print Dipole_py[absD > 1.0e-3]
    print Dipole_f90[absD > 1.0e-3]
    errD = (abs(Dipole_py-Dipole_f90)).max()
    assert errD < 1.0e-3, "python and Fortran implementations disagree for Dipoles, error = %s" % errD

    print "Timing:"
    print "Fortran: %s seconds" % time_f90
    print "Python : %s seconds" % time_py
    exit(-1)
    """
    if SLAKO_IMPLEMENTATION == "F90":
        # fortran code, parallelized with OpenMP
        Dipole = _DipoleMatrix_f90(atomlist, atomlist, valorbs, SKT, Mproximity, S)
    else:
        # serial code iny python
        Dipole = _DipoleMatrix(atomlist, atomlist, valorbs, SKT, Mproximity)
    return Dipole

def _DipoleMatrix(atomlist1, atomlist2, valorbs, SKT, Mproximity):
    """
    compute matrix elements of dipole operator between valence orbitals
    using Slater-Koster Rules

    Parameters:
    ===========
    atomlist1, atomlist2: list of tuples (Zi,[xi,yi,zi]) of atom types and positions
    valorbs: list of valence orbitals with quantum numbers (ni,li,mi)
    SKT: Slater Koster table
    Mproximity: M[i,j] == 1, if the atoms i and j are close enough 
      so that the overlap between orbitals on i and j should be 
      computed
    
    Returns:
    ========
    Dipoles: matrix with shape (Norb,Norb,3)
       Dipoles[i,j,0] for instance would be <i|x|j>

    """
    # count valence orbitals
    Norb1 = count_orbitals(atomlist1, valorbs)
    Norb2 = count_orbitals(atomlist2, valorbs)
    
    Dipole = np.zeros((Norb1,Norb2,3))

    # iterate over atoms
    mu = 0
    for i,(Zi,posi) in enumerate(atomlist1):
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbs[Zi]:
            # iterate over atoms
            nu = 0
            for j,(Zj,posj) in enumerate(atomlist2):
                # iterate over orbitals on center j
                for (nj,lj,mj) in valorbs[Zj]:
                    if Mproximity[i,j] == 1:
                        if mu <= nu:
                            if Zi <= Zj:
                                # the first atom given to getDipole()
                                # has to be always the one with lower atomic number
                                Dx,Dy,Dz  = SKT[(Zi,Zj)].getDipole(li,mi,posi, lj,mj,posj)
                            else:
                                # swap atoms if Zj > Zi
                                Dx,Dy,Dz  = SKT[(Zj,Zi)].getDipole(lj,mj,posj, li,mi,posi)
                            Dipole[mu,nu,0] = Dx
                            Dipole[mu,nu,1] = Dy
                            Dipole[mu,nu,2] = Dz
                        else:
                            # Dipole matrix is symmetric 
                            Dipole[mu,nu,0] = Dipole[nu,mu,0]
                            Dipole[mu,nu,1] = Dipole[nu,mu,1]
                            Dipole[mu,nu,2] = Dipole[nu,mu,2]
                    nu += 1
            mu += 1
    return Dipole

# FASTER DIPOLE MATRIX ELEMENTS WITH FORTRAN
def _DipoleMatrix_f90(atomlist1, atomlist2, valorbs, SKT, Mproximity, S):
    """
    This function calls an external Fortran function that computes the matrix elements.
    Since hashes and tuples cannot be used easily in Fortran, the arguments are brought
    into a form that can be fed into the Fortran function.
    """
    atom_type_dic, spline_deg, tab_filled_SH0, tab_filled_D, \
        (S_knots, S_coefs, H_knots, H_coefs, D_knots, D_coefs) = \
                            combine_slako_tables_f90(SKT)
    
    # number of atoms
    Nat1 = len(atomlist1)
    Nat2 = len(atomlist2)
    # orbitals
    atom_indeces1,atom_types1,ls1,ms1 = atomlist2orbitals(atomlist1, valorbs, atom_type_dic)
    atom_indeces2,atom_types2,ls2,ms2 = atomlist2orbitals(atomlist2, valorbs, atom_type_dic)
    assert atom_indeces1 == atom_indeces2, "atomlist1 and atomlist2 should contain the same atoms in the same order!"
    # count valence orbitals
    Norb1 = len(ls1)
    Norb2 = len(ls2)
    # distances and direction
    pos1 = XYZ.atomlist2vector(atomlist1)
    pos1 = np.reshape(pos1,(Nat1,3)).transpose()  # pos(:,i) is the 3d position of atom i
    pos2 = XYZ.atomlist2vector(atomlist2)
    pos2 = np.reshape(pos2,(Nat2,3)).transpose() 
    r,x,y,z = slako.slako.directional_cosines(pos1,pos2)

    Dipole = slako.slako.dipolematrix(atom_indeces1,atom_types1,ls1,ms1,
                              atom_indeces2,atom_types2,ls2,ms2,
                              r,x,y,z,pos1,
                              Mproximity,
                              S,
                              spline_deg, tab_filled_D,
                              D_knots, D_coefs)
    # In fortran the faster indeces come first, roll axes D(xyz,mu,nu) -> D(mu,nu,xyz)
    Dipole = np.rollaxis(Dipole,0,3)
    return Dipole


def OnSiteMatrix(atomlist, valorbs):
    """
    fill in matrix with 1-center electron integrals for on-site correction

       G[m,n] = (mn|mn)   for  m != n and m,n on the same atom
    
    """
    from DFTB.SlaterKoster.free_pseudo_atoms import onsite_electron_integrals
    from DFTB.AtomicData import atom_names
    
    Norb = count_orbitals(atomlist, valorbs)
    G_onsite = np.zeros((Norb,Norb))

    # iterate over atoms
    mu = 0
    for i,(Zi,posi) in enumerate(atomlist):
        # on-site integrals for atom i
        ssss, spsp, sspp, pppp, pqpq = onsite_electron_integrals.onsite_integrals[atom_names[Zi-1]]
        # only (sp|sp) and (pp'|pp') are needed
        
        # iterate over valence orbitals on atom i
        nu0 = mu
        for (ni,li,mi) in valorbs[Zi]:
            nu = nu0
            for (nj,lj,mj) in valorbs[Zi]:
                if li > 1 or lj > 1:
                    # on-site correction is only for s- and p-type orbitals
                    nu += 1
                    continue
                if mu == nu:
                    # on-site correction is only for different orbitals
                    nu += 1
                    continue
                #
                if li != lj:
                    # (sp|sp)
                    G_onsite[mu,nu] = spsp
                elif li == lj:
                    # (pp'|pp')
                    G_onsite[mu,nu] = pqpq
                    
                nu += 1
            mu += 1
                
    return G_onsite
