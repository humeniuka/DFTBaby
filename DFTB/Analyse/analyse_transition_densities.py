"""
Analyse the transition density of a linear polymer in terms of the transition
densities of the monomer. 
It is important that the hydrogen atoms come at the very end in the geometries of the monomer
and the polymers. 
"""

import numpy as np
import numpy.linalg as la
from matplotlib import pyplot as plt
import glob

from DFTB import XYZ, utils

def count_hydrogens(atomlist):
    hydrogen_indeces = []
    # count number of hydrogens
    nhyd = 0
    for i,(Zi, posi) in enumerate(atomlist):
        if Zi == 1:
            nhyd += 1
    #
    for (Zi, posi) in atomlist[i:]:
        assert Zi == 1, "All hydrogen atoms should be placed at the end, got Z=%d" % Zi
    return nhyd

if __name__ == "__main__":
    import sys
    from optparse import OptionParser
    
    usage = "Usage: python %s <monomer geometry .xyz> <polymer geometry .xyz> <pattern for monomer transition densities> <pattern for polymer transition densities>\n" % sys.argv[0]
    usage += "  The excited states of the polymer are analysed in terms of the excited states of the monomer\n"

    parser = OptionParser(usage)

    (opts, args) = parser.parse_args()
    if len(args) < 4:
        print usage
        exit(-1)

    atomlist_mono = XYZ.read_xyz(args[0])[0]
    atomlist_poly = XYZ.read_xyz(args[1])[0]
    monomer_pattern = args[2]
    polymer_pattern = args[3]
    
    print "load transition densities for monomer states"
    Ptrans_mono = []  # list of transition densities matrices in AO basis, one for each state
    for tdense_file in glob.glob(monomer_pattern+"*.mat"):
        print tdense_file
        PtransI = np.loadtxt(tdense_file)
        # remove hydrogen orbitals
        nao,nao = PtransI.shape
        nhyd = count_hydrogens(atomlist_mono)
        print "remove %d hydrogens" % nhyd
        PtransI = PtransI[:(nao-nhyd),:][:,:(nao-nhyd)]
        Ptrans_mono.append(PtransI)
    print "load transition densities for polymer states"
    Ptrans_poly = []  # list of transition densities matrices in AO basis, one for each state
    for tdense_file in glob.glob(polymer_pattern+"*.mat"):
        print tdense_file
        PtransI = np.loadtxt(tdense_file)
        # remove hydrogen orbitals
        nao,nao = PtransI.shape
        nhyd = count_hydrogens(atomlist_poly)
        print "remove %d hydrogens" % nhyd
        PtransI = PtransI[:(nao-nhyd),:][:,:(nao-nhyd)]
        Ptrans_poly.append(PtransI)
        
    # check that transition density matrices form orthonormal basis
    for Ptrans in [Ptrans_mono, Ptrans_poly]:
        Nst = len(Ptrans)
        S = np.zeros((Nst,Nst))
        for i in range(0, Nst):
            for j in range(0, Nst):
                S[i,j] = np.dot(Ptrans[i].flatten(), Ptrans[j].flatten())
        # relative error
        err =  np.sum(abs(S - np.eye(Nst))) / np.sum(abs(S))
        #assert err < 1.0e-10, "err = %s" % err
#        assert err < 1.0e-1, "err = %s" % err  # Since hydrogen orbitals have been removed, the P's are not exactly orthogonal anymore


    Nst_mono = len(Ptrans_mono)
    Nst_poly = len(Ptrans_poly)

    # number of AOs in the monomer
    nAOm,nAOm = Ptrans_mono[0].shape
    print "nAOm = %s" % nAOm
    # number of AOs in polymer
    nAOp,nAOp = Ptrans_poly[0].shape
    print "nAOp = %s" % nAOp
    assert nAOp % nAOm == 0
    # number of monomeric units in polymer
    nunits = nAOp/nAOm
    print "Linear polymer consists of %d monomeric units" % nunits

    for I in range(0, Nst_poly):
        PI = Ptrans_poly[I]
        
#        plt.matshow(PI)
        
        print "Analyse state %d of polymer" % I
        MI = np.zeros((Nst_mono, nunits))
        for n in range(0, nunits):
            # extract diagonal block of the n-th subunit
            Pn = PI[n*nAOm:(n+1)*nAOm, n*nAOm:(n+1)*nAOm]
            vvec = Pn.flatten()

            for Im in range(0, Nst_mono):
                bvec = Ptrans_mono[Im].flatten()
                MI[Im,n] = np.dot(bvec, vvec)
            
        unit_labels = ["U-%d" % nu for nu in range(0, nunits)]
        state_labels = ["St-%d" % st for st in range(0, Nst_mono)]
        print utils.annotated_matrix(MI, state_labels, unit_labels)
        plt.matshow(MI, interpolation="nearest", cmap=plt.cm.jet)
        #plt.show()
    plt.show()


    """
    plt.matshow(Ptrans_mono[0], interpolation="nearest", cmap=plt.cm.jet)
    for Ptrans in Ptrans_poly:
        plt.matshow(Ptrans, interpolation="nearest", cmap=plt.cm.jet)
    plt.show()
    """
