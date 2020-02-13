"""
create Wigner distributions from Gaussian 09 or Turbomole Hessians
"""
from scipy import optimize
from scipy.misc import factorial
import numpy as np
import numpy.linalg as la
import scipy.linalg as sla

import sympy

from DFTB import AtomicData, XYZ
from DFTB.Formats.Gaussian2py import GaussianWavepacket
from DFTB.Formats.Gaussian2py import Gaussian
from DFTB.Formats.Gaussian2py import Checkpoint

def center_of_mass(pos, masses):
    """
    Parameters:
    ===========
    x: 
    """
    # xm = sum_i ri * mi
    xm = pos*masses
    M = np.sum(masses[0::3])
    assert abs(np.sum(masses[1::3]) - np.sum(masses[2::3])) < 1.0e-8
    com = np.array([np.sum(xm[0::3]), np.sum(xm[1::3]), np.sum(xm[2::3])])/M
    return com
    

def tensor_of_inertia(pos, masses):
    """
    find tensor of inertia which relates angular velocity omega and angular momentum
    of a rigid body: L = I*w
    """
    I = np.zeros((3,3))
    for i in xrange(0, len(pos)/3):
        xi,yi,zi = pos[3*i:3*i+3]
        ri2 = xi*xi+yi*yi+zi*zi
        mi = masses[3*i]
        I[0,0] += mi*(ri2-xi*xi)
        I[0,1] -= mi*xi*yi
        I[0,2] -= mi*xi*zi
        I[1,1] += mi*(ri2-yi*yi)
        I[1,2] -= mi*(yi*zi)
        I[2,2] += mi*(ri2-zi*zi)
    I[1,0] = I[0,1]
    I[2,0] = I[0,2]
    I[2,1] = I[1,2]
    return I

def shift_to_com(pos,masses):
    """
    shift center of mass to the origin
    """
    com = center_of_mass(pos, masses)
    pos_shifted = np.zeros(pos.shape)
    for i in xrange(0, len(pos)/3):
        pos_shifted[3*i:3*i+3] = pos[3*i:3*i+3] - com
    return pos_shifted

def expected_zero_modes(x0, masses):
    """
    How many eigen values of the Hessian should be zero
    because of the structure of the molecule (atom, linear, polyatomic)?
    """
    # assume that coordinate vector with 3*N components belongs to a molecule
    # with N atom
    assert len(x0) % 3 == 0
    Nat = len(x0) / 3
    # shift the origin to the center of mass
    x0_shift = shift_to_com(x0, masses)
    # diagonalize the tensor of inertia to obtain the principal moments and 
    # the normalized eigenvectors of I
    Inert = tensor_of_inertia(x0_shift, masses)
    principle_moments, X = la.eigh(Inert)
    # check that the number of rotational and translational modes
    # is what is expected:
    #   single atom  => 3 translational modes only = 3
    #   dimer or linear molecule => 3 translational and 2 rotational modes = 5
    #   otherwise => 3 translational and 3 rotational modes = 6
    
    # If the molecule is linear two of the principle moments of inertia
    # will be zero
    Ntrans = 3    # number of translational
    Nrot = 3      # and rotational modes which should be very close to zero
    is_linear = False
    pmom_abs = np.sort(abs(principle_moments))
    # In a linear triatomic molecule we have Icc = Ibb > Iaa = 0
    if abs(pmom_abs[2]-pmom_abs[1]) < 1.0e-6 and abs(pmom_abs[0]) < 1.0e-6:
        is_linear = True
        print "Molecule is linear"
    if Nat == 1:
        Nrot = 0
    elif is_linear == True or Nat == 2:
        Nrot = 2
    else:
        Nrot = 3
        
    return (Ntrans+Nrot)


def vibrational_analysis(hess, masses, zero_threshold=1.0e-9, is_molecule=False):
    """
    compute vibrational frequencies and modes from Hessian

    Parameters:
    ===========
    hess: Hessian in cartesian coordinates
    masses: Ndim-vector with masses for each cartesian coordiante
    zero_threshold: modes with frequencies^2 below that threshold
      are treated as zero
    is_molecule: If this flag is set, the program checks whether the number
      of zero modes is as expected for that molecule

    Returns:
    ========
    freq: vibrational frequencies in atomic units
    modes: displacement vectors in mass-weighted coordinates
    """
    # convert Hessian to mass-weighted coordinates
    hess_mwc = hess / np.sqrt(np.outer(masses, masses))
    print "HESSIAN for mass weighted cartesian coordinates"
    print "==============================================="
    print hess_mwc
    # mass weighted coordinates are now qi = sqrt(mi) dxi
    # compute eigen values of hess_mwc
    omega2,modes = la.eigh(hess_mwc)

    # modes that are zero within numerical accuracy
    zero_modes = np.where((omega2 < zero_threshold))[0]
    vib_modes = np.where((omega2 >= zero_threshold))[0]
    Nzero = len(zero_modes)

    if is_molecule == True:
        Nzero_expected = expected_zero_modes(x0, masses)
        if Nzero != Nzero_expected:
            raise Exception("Expected %d modes with 0 frequency (translation + rotation) but got %d modes" % (Nzero_expected, Nzero))

    freqs = np.sqrt(omega2*(1.0+0.0j))
    print "Frequencies"
    print "==========="
    print "- Zero modes (should be close to zero)"
    for fr in freqs[zero_modes]:
        print "   %s    " % (str(fr).ljust(15))
    print "- Vibrations"
    for fr in freqs[vib_modes].real:
        print "   %5.7f Hartree      %5.7f cm-1" % (fr, fr*AtomicData.hartree_to_wavenumbers)

    return freqs[vib_modes], modes[:,vib_modes]

def wigner_distribution(x0, hess, masses, zero_threshold=1.0e-9, is_molecule=False):
    """
    compute wigner distribution for quadratic potential. 

    The hessian can have zero vibrational frequencies (translation and rotation) which
    whould result in a very broad wave packet along those zero modes. Therefore only
    the non-zero modes are transformed to the Wigner representation while the zero modes
    are constrained by delta functions:

      W(Q1,...,QN;P1,...,PN) ~ Prod(vib.modes i) exp(-Omega_i*Qi^2 - 1/Omega_i * Pi^2)
                              *Prod(zero modes j) delta(Qj)*delta(Pj)


    Returns:
    ========
    Aw, Bw
    """
    # convert Hessian to mass-weighted coordinates
    hess_mwc = hess / np.sqrt(np.outer(masses, masses))
    # mass weighted coordinates are now qi = sqrt(mi) dxi
    # compute eigen values of hess_mwc
    omega2,modes = la.eigh(hess_mwc)

    # modes that are zero within numerical accuracy
    zero_modes = np.where((omega2 < zero_threshold))[0]
    vib_modes = np.where((omega2 >= zero_threshold))[0]
    Nzero = len(zero_modes)

    if is_molecule == True:
        Nzero_expected = expected_zero_modes(x0, masses)
        if Nzero != Nzero_expected:
            raise Exception("Expected %d modes with 0 frequency (translation + rotation) but got %d modes" % (Nzero_expected, Nzero))

    Ndim = len(masses)
    Aq = np.zeros(hess.shape)#, dtype=complex)
    Ap = np.zeros(hess.shape)#, dtype=complex)

    # vibrational modes
    Msq = np.sqrt(np.outer(masses, masses))
    MsqInv = 1.0/Msq
    Oi = np.sqrt(omega2[vib_modes])
    Li = modes[:,vib_modes]
    print "Wigner distribution for vibrational modes"
    Oii = np.diag(Oi)
    OiiInv = np.diag(1.0/Oi)
    Aq = 2 * np.dot(Li, np.dot(Oii, Li.transpose())) * Msq
    Ap = 2 * np.dot(Li, np.dot(OiiInv, Li.transpose())) * MsqInv
    # constrain zero modes by delta-functions  
    #   delta(Qi)*delta(Pi) ~ lim_(Oconstr->infty) exp(-Oconstr*(Qi^2 + Pi^2))
    print "delta distribution for zero modes"
    Oconstr = 1.0e6 * np.eye(len(zero_modes))  # very large number, e.g. 1.0e10
    Li = modes[:,zero_modes]
    Aq += np.dot(Li, np.dot(Oconstr, Li.transpose())) * Msq
    Ap += np.dot(Li, np.dot(Oconstr, Li.transpose())) * MsqInv * Msq.max()/MsqInv.max()

    Bq = np.dot(x0, Aq)
    Bp = np.zeros(x0.shape)

    Zo = np.zeros((Ndim, Ndim))
    Aw = np.bmat([[Aq, Zo],
                  [Zo, Ap]])
    Aw = np.asarray(Aw)
    Bw = np.hstack([Bq, Bp]).transpose()
    
    return Aw, Bw


############ Wigner distributions from ab-initio calculations ###########
def wigner_from_G09_hessian(g09_file, Nsample=100, zero_threshold=1.0e-9):
    """
    create Wigner ensemble based on hessian matrix from Gaussian 09 calculation
    """
    suffix = g09_file.split(".")[-1]
    if suffix in ["out", "log"]:
        print "Reading Gaussian 09 log file %s" % g09_file
        atomlist = Gaussian.read_geometry(g09_file)
        forces = Gaussian.read_forces(g09_file)
        hess = Gaussian.read_force_constants(g09_file)
    elif suffix in ["fchk"]:
        print "Reading formatted Gaussian 09 checkpoint file %s" % g09_file
        Data = Checkpoint.parseCheckpointFile(g09_file)
        # cartesian coordinates
        pos = Data["_Current_cartesian_coordinates"]
        atnos = Data["_Atomic_numbers"]
        # forces
        frc = - Data["_Cartesian_Gradient"]
        atomlist = []
        forces = []
        for i,Zi in enumerate(atnos):
            atomlist.append( (Zi, tuple(pos[3*i:3*(i+1)])) )
            forces.append( (Zi, tuple(frc[3*i:3*(i+1)])) )
        # Hessian
        hess = Data["_Cartesian_Force_Constants"]

    masses = np.array(AtomicData.atomlist2masses(atomlist))
    
    x0 = XYZ.atomlist2vector(atomlist)
    x0 = shift_to_com(x0, masses)

    grad = -XYZ.atomlist2vector(forces)
    
    grad_nrm = la.norm(grad)
    print "  gradient norm = %s" % grad_nrm
#    assert grad_nrm < 1.0e-3, "Gradient norm too large for minimum!"
    vib_freq, vib_modes = vibrational_analysis(hess, masses, zero_threshold=zero_threshold)
    Aw,Bw = wigner_distribution(x0, hess, masses, zero_threshold=zero_threshold)
    gw = GaussianWavepacket.Gaussian(Aw,Bw)
    qs, ps = gw.sample(Nsample)
    mx = np.outer(masses, np.ones(Nsample)) * qs
    avg_com = np.mean(np.sum(mx[::3,:], axis=0)), np.mean(np.sum(mx[1::3,:], axis=0)), np.mean(np.sum(mx[2::3,:], axis=0))
    print avg_com

    geometries = [XYZ.vector2atomlist(qs[:,i], atomlist) for i in range(0, Nsample)]

    return geometries

if __name__ == "__main__":
    from DFTB import XYZ

    import sys
    import optparse
    usage = "python %s <G09 .log of .fchk file> <wigner ensemble .xyz>\n" % sys.argv[0]
    usage += "  reads the Hessian from the log or formatted checkpoing file and\n"
    usage += "  writes a Wigner ensemble to the specified xyz-file.\n"

    parser = optparse.OptionParser(usage)
    parser.add_option("--Nsample", dest="Nsample", type=int, default=30, help="Number of samples from the Wigner distribution [default: %default]")

    (opts, args) = parser.parse_args()
    if len(args) < 2:
        print usage
        exit(-1)
    g09_file = args[0]
    wigner_xyz = args[1]
    geometries = wigner_from_G09_hessian(g09_file, Nsample=opts.Nsample)

    XYZ.write_xyz(wigner_xyz, geometries)
    print "Wigner ensemble written to file %s" % wigner_xyz
