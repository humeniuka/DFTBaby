""""
approximate non-adiabatic coupling matrix elements by the overlap
of the electronic wavefunctions at different positions

V_IJ = <Psi_I(t)|Psi_K(t+dt)>

see Mitric et.al., J. Chem. Phys. A, 113, 2009
"""
import numpy as np
from numpy import linalg as la

from DFTB.SKMatrixElements import SLAKO_IMPLEMENTATION, atomlist2orbitals

def count_valence_orbitals(atomlist, valorbs):
    # count valence orbitals
    Norb = 0
    for i,(Zi,posi) in enumerate(atomlist):
        Norb += len(valorbs[Zi])
    return Norb

def overlap_AB(atomlistA, atomlistB, valorbsA, valorbsB, SKT):
    """
    # compare python and fortran implementations
    # parallelized fortran code
    import time
    ta = time.time()
    S_f90 = _overlap_AB_f90(atomlistA, atomlistB, valorbsA, valorbsB, SKT)
    tb = time.time()
    time_f90 = tb-ta
    # serial code
    ta = time.time()
    S_py = _overlap_AB(atomlistA, atomlistB, valorbsA, valorbsB, SKT)
    tb = time.time()
    time_py = tb-ta
    # compare
    from DFTB import utils
    orb_labels = ["%d" % i for i in range(0, S_f90.shape[0])]
    print "S F90"
    print utils.annotated_matrix(S_f90, orb_labels, orb_labels)
    print "S py"
    print utils.annotated_matrix(S_py, orb_labels, orb_labels)
    errS = (abs(S_py-S_f90)).max()
    assert errS < 1.0e-3, "python and Fortran implementations disagree for S, error = %s" % errS

    print "Timing:"
    print "Fortran: %s seconds" % time_f90
    print "Python : %s seconds" % time_py
    exit(-1)
    """
    if SLAKO_IMPLEMENTATION == "F90":
        # fortran implementation
        S = _overlap_AB_f90(atomlistA, atomlistB, valorbsA, valorbsB, SKT)
    else:
        # python implementation
        S = _overlap_AB(atomlistA, atomlistB, valorbsA, valorbsB, SKT)
    return S

def _overlap_AB(atomlistA, atomlistB, valorbsA, valorbsB, SKT):
    """
    compute overlap matrix elements between two sets of atoms using
    Slater-Koster rules

    Parameters:
    ===========
    valorbs: 
    atomlistA, atomlistB: list of (Zi,(xi,yi,zi)) for each atom
    valorbsA, valorbsB: list of valence orbitals with quantum numbers (ni,li,mi)
    SKT: Slater Koster table

    Returns:
    ========
    S_AB: overlap between the atomic orbitals of structure A and B
    """
    # count valence orbitals
    NorbA = count_valence_orbitals(atomlistA, valorbsA)
    NorbB = count_valence_orbitals(atomlistB, valorbsB)
    S = np.zeros((NorbA,NorbB))
    # iterate over atoms
    mu = 0
    for i,(Zi,posi) in enumerate(atomlistA):
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbsA[Zi]:
            # iterate over atoms
            nu = 0
            for j,(Zj,posj) in enumerate(atomlistB):
                # iterate over orbitals on center j
                for (nj,lj,mj) in valorbsB[Zj]:
                    if Zi <= Zj:
                        s  = SKT[(Zi,Zj)].getOverlap(li,mi,posi, lj,mj,posj)
                    else:
                        # swap atoms if Zj > Zi
                        s  = SKT[(Zj,Zi)].getOverlap(lj,mj,posj, li,mi,posi)
                    S[mu,nu] = s
                    nu += 1
            mu += 1
    return S

# FASTER FORTRAN CODE
from DFTB import XYZ
from DFTB.extensions import slako
from DFTB.SlaterKoster.SKIntegrals import combine_slako_tables_f90

def _overlap_AB_f90(atomlist1, atomlist2, valorbs1, valorbs2, SKT):
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
    atom_indeces1,atom_types1,ls1,ms1 = atomlist2orbitals(atomlist1, valorbs1, atom_type_dic)
    atom_indeces2,atom_types2,ls2,ms2 = atomlist2orbitals(atomlist2, valorbs2, atom_type_dic)
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
    S = slako.slako.overlap12(atom_indeces1,atom_types1,ls1,ms1,
                              atom_indeces2,atom_types2,ls2,ms2,
                              r,x,y,z,
                              spline_deg, tab_filled_SH0,
                              S_knots, S_coefs)
    return S
#

class ScalarCoupling:
    def __init__(self, tddftb):
        self.tddftb = tddftb
    def ci_overlap(self, atomlist1, orbs1, C1, atomlist2, orbs2, C2, Nst, \
                       threshold=1.0e-3):
        """
        Parameters:
        ===========
        atomlist1, atomlist2: two different geometries of the same molecule
        orbs1, orbs2: KS orbitals for geometries 1 and 2
        C1, C2: TD-DFT coefficients for geometry 1 and 2
           C[I,i,a] is the coefficients of the single excitation i->a
           in the TD-DFTB "wavefunction" of the (I+1)-th excited state
        Nst: number of states to consider (including ground state)
        threshold: excitations with |amplitudes| below this threshold
           are not considered

        atomlist1 and atomlist2 should contain the same atoms
        in the same order as in dftb.atomlist
        """
        if self.tddftb.dftb2.verbose > 0:
            print "Compute CI overlap between TD-DFT 'wavefunctions'"
            print "Excitations i->a with coefficients |C_ia| < %e will be neglected" % threshold
        # QM/MM
        qmmm = self.tddftb.dftb2.qmmm
        if qmmm != None:
            # filter the atomlist so that we only get the QM part
            atomlist1 = qmmm.partitionGeometry(atomlist1)
            atomlist2 = qmmm.partitionGeometry(atomlist2)
        #
        valorbs = self.tddftb.dftb2.valorbs
        SKT = self.tddftb.dftb2.SKT
        # overlap between atomlist orbitals
        Sao = overlap_AB(atomlist1, atomlist2, valorbs, valorbs, SKT)
        # overlap between molecular orbitals
        Smo = np.dot(orbs1.transpose(), np.dot(Sao, orbs2))
        # overlap between CI-wavefunctions 1 and 2
        Sci = np.zeros((Nst, Nst))
        # indeces of occupied spatial orbitals
        Nelec = self.tddftb.dftb2.Nelec_val
        assert Nelec % 2 == 0
        # ground state occupation
        occs_0 = np.array([i for i in range(0, int(Nelec/2))])
        # indices active occupied and virtual orbitals 
        occ, virt = self.tddftb.getActiveOrbitals()
        nocc = len(occ)
        nvirt = len(virt)

        S_ij = Smo[occs_0,:][:,occs_0]
        det_ij = la.det(S_ij)
        
        # overlaps between excited states  Sci[1:,1:]
        for i in range(0, nocc):
            for a in range(0, nvirt):
                #
                if abs(C1[:,i,a]).max() < threshold:
                    # skip the i->a excitation
                    continue
                # occupied orbitals in the configuration state function
                # |Psi_ia>
                occs_ia = np.copy(occs_0)
                occs_ia[occ[i]] = virt[a]
                # <1,...,a,...|1,...,j,...>
                S_aj = Smo[occs_ia,:][:,occs_0]
                det_aj = la.det(S_aj)

                for j in range(0, nocc):
                    for b in range(0, nvirt):
                        if abs(C2[:,j,b]).max() < threshold:
                            # skip the j->b excitation
                            continue
                        # occupied orbitals in the configuration state function
                        # |Psi_jb>
                        occs_jb = np.copy(occs_0)
                        occs_jb[occ[j]] = virt[b]
                        # select part of overlap matrix for orbitals
                        # in |Psi_ia> and |Psi_jb>
                        # <1,...,a,...|1,...,b,...>
                        S_ab = Smo[occs_ia,:][:,occs_jb]
                        det_ab = la.det(S_ab)
                        # <1,...,i,...|1,...,b,...>
                        S_ib = Smo[occs_0,:][:,occs_jb]
                        det_ib = la.det(S_ib)
                        # loop over excited states
                        for I in range(1, Nst):
                            for J in range(1, Nst):
                                cc = C1[I-1,i,a].conjugate() * C2[J-1,j,b]
                                # see eqn. (9.39) in
                                #  A. Humeniuk, PhD thesis (2018),
                                #  "Methods for Simulating Light-Induced Dynamics in Large Molecular Systems"
                                Sci[I,J] += cc * (det_ab * det_ij + det_aj * det_ib)
        # overlaps between ground state <Psi0|PsiJ'> and excited states
        for i in range(0, nocc):
            for a in range(0, nvirt):
                #
                if abs(C2[:,i,a]).max() < threshold: 
                    # skip the i->a excitation
                    continue
                # occupied orbitals in the configuration state function
                # |Psi_ia>
                occs_ia = np.copy(occs_0)
                occs_ia[occ[i]] = virt[a]
                S_ia = Smo[occs_0,:][:,occs_ia]
                det_ia = la.det(S_ia)
                for J in range(1,Nst):
                    c0 = C2[J-1,i,a]
                    Sci[0,J] += c0 * np.sqrt(2) * (det_ia * det_ij)
        # overlaps between ground state <PsiI|Psi0'> and excited states
        for i in range(0, nocc):
            for a in range(0, nvirt):
                #
                if abs(C1[:,i,a]).max() < threshold:
                    # skip the i->a excitation
                    continue
                # occupied orbitals in the configuration state function
                # |Psi_ia>
                occs_ia = np.copy(occs_0)
                occs_ia[occ[i]] = virt[a]
                S_ia = Smo[occs_ia,:][:,occs_0]
                det_ia = la.det(S_ia)
                for I in range(1,Nst):
                    c0 = C1[I-1,i,a]
                    Sci[I,0] += c0 * np.sqrt(2) * (det_ia * det_ij)
                    
        # overlap between ground states <Psi0|Psi0'>
        Sci[0,0] = det_ij

        return Sci
