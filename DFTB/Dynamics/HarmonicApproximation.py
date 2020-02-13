"""
Approximate a potential energy surface V(x) around its minimum x0
by a second order polynomial. 
  V(x) ~ V(x0) + 1/2 sum_ij Hij (x-x0)_i (x-x0)_j
"""
from scipy import optimize
from scipy.misc import factorial
import numpy as np
import numpy.linalg as la
import scipy.linalg as sla

from os.path import join

from DFTB import AtomicData, XYZ
from DFTB.Modeling import MolecularCoords as MolCo
from DFTB.Dynamics import GaussianWavepacket


def numerical_hessian_G(grad,x0,h=1.0e-5):
    """
    compute hessian by numerical differentiation of gradients

    Parameters:
    ===========
    grad: grad(x) computes the gradient at position x
    """
    n = len(x0)
    hess = np.zeros((n,n))
    g0 = grad(x0)
    print "numerical differentiation of forces requires %d gradient calculations" % (2*n)
    for i in range(0, n):
        print "%s of %d" % (2*(i+1), 2*n)
        ei = np.zeros(n)
        ei[i] = 1.0
        x_phi = x0 + h*ei
#        hess[i,:] = (grad(x_phi) - g0)/h       
        x_mhi = x0 - h*ei
        hess[i,:] = (grad(x_phi) - grad(x_mhi))/(2*h)

    # hessian should be symmetric
    hess = 0.5 * (hess + hess.transpose())
    
    return hess


##############################################

def expected_zero_modes(x0, masses):
    """
    How many eigen values of the Hessian should be zero
    because of the structure of the molecule (atom, linear, polyatomic)?
    """
    # assume that coordinate vector with 3*N components belongs to a molecule
    # with N atoms
    assert len(x0) % 3 == 0
    Nat = len(x0) / 3
    # shift the origin to the center of mass
    x0_shift = MolCo.shift_to_com(x0, masses)
    # diagonalize the tensor of inertia to obtain the principal moments and 
    # the normalized eigenvectors of I
    Inert = MolCo.inertial_tensor(masses, x0_shift)
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


def vibrational_analysis(x0, hess, masses, zero_threshold=1.0e-9, is_molecule=False):
    """
    compute vibrational frequencies and modes from Hessian

    Parameters:
    ===========
    x0: 
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
            print "WARNING: Expected %d modes with 0 frequency (translation + rotation) but got %d modes!" % (Nzero_expected, Nzero)

    freqs = np.sqrt(omega2*(1.0+0.0j))
    print "Frequencies"
    print "==========="
    print "- Zero modes (should be close to zero)"
    for fr in freqs[zero_modes]:
        print "   %s    " % (str(fr).ljust(15))
    print "- Vibrations"
    for fr in freqs[vib_modes].real:
        print "   %5.7f Hartree      %5.7f cm-1" % (fr, fr*AtomicData.hartree_to_wavenumbers)
    print ""
    en_zp = np.sum(freqs[vib_modes].real)/2.0
    print "zero-point energy:   %5.7f Hartree      %5.7f cm-1" % (en_zp, en_zp*AtomicData.hartree_to_wavenumbers)
    return freqs[vib_modes], modes[:,vib_modes]
#    return freqs, modes

def ho_wavefunction(x0, hess, masses):
    """
    computes the Gaussian that represents the ground state of 
    the harmonic oscillator potential defined by the hessian at x0.

    Returns:
    ========
    Apsi,Bpsi: matrix A and vector B such that
     psi(x1,...,xD) ~ exp(-1/2 x.A.x + B.x)
    """
    # sqrt(mi*mj)
    Msq = np.sqrt(np.outer(masses, masses))
    # convert Hessian to mass-weighted coordinates
    hess_mwc = hess / Msq
    # compute matrix square root of mass-weighted hessian and
    # multiply again by masses
    # psi(x1,...,xD) ~ exp(-1/2 (x-x0)^T.A.(x-x0)) 
    Apsi = sla.sqrtm(hess_mwc) * Msq
    Bpsi = np.dot(x0, Apsi)
    return Apsi,Bpsi

def wigner_distribution(x0, hess, masses, zero_threshold=1.0e-9, is_molecule=True):
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
            print "WARNING: Expected %d modes with 0 frequency (translation + rotation) but got %d modes!" % (Nzero_expected, Nzero)

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

def initial_conditions_wigner(x0, hess, masses, Nsample=50, zero_threshold=1.0e-8, power_wigner=1.0):
    """
    sample positions and momenta from Wigner distribution
    """
    vib_freq, vib_modes = vibrational_analysis(x0, hess, masses, zero_threshold=zero_threshold)
    Aw,Bw = wigner_distribution(x0, hess, masses, zero_threshold=zero_threshold)
    gw = GaussianWavepacket.Gaussian(Aw,Bw)
    qs, ps = gw.sample_wigner(Nsample, pw=power_wigner)
    return qs, ps

def dynamics_in_format(atomlist, q, p, fname):
    atomlist_q = XYZ.vector2atomlist(q, atomlist)
    masses = AtomicData.atomlist2masses(atomlist)
    Nat = len(atomlist_q)
    txt = "%d\n" % Nat
    # positions in bohr
    for Zi, (x,y,z) in atomlist_q:
        atname = AtomicData.atom_names[Zi-1]
        txt += "%4s    %+10.15f    %+10.15f    %+10.15f\n" % (atname, x, y, z)
    # velocities in bohr/s
    vel = p/masses
    for i in range(0, Nat):
        vx,vy,vz = vel[3*i:3*(i+1)]
        txt += "   %+10.15f    %+10.15f   %+10.15f\n" % (vx,vy,vz)
    return txt

def qmf3_in_format(atomlist, q, p, fname):
    """
    Save initial conditions in the format expected by Matthias' QMF3 program
    """
    atomlist_q = XYZ.vector2atomlist(q, atomlist)
    masses = AtomicData.atomlist2masses(atomlist)
    Nat = len(atomlist_q)
    txt = "$coordinates [\n"
    # positions in bohr
    for Zi, (x,y,z) in atomlist_q:
        atname = AtomicData.atom_names[Zi-1]
        txt += "%4s    %+10.15f    %+10.15f    %+10.15f\n" % (atname, x, y, z)
    txt += "]\n"
    # velocities in bohr/s
    vel = p/masses
    txt += "$velocities [\n"
    for i in range(0, Nat):
        vx,vy,vz = vel[3*i:3*(i+1)]
        txt += "   %+10.15f    %+10.15f   %+10.15f\n" % (vx,vy,vz)
    txt += "]\n"
    txt += "$end\n"
    return txt

    
def save_initial_conditions(atomlist, qs, ps, inidir, name, frmt="FISH"):
    """
    save initial conditions in the format used by Jens' program

    The names will be inidir/name_####.in

    Parameters:
    ===========
    inidir: directory where the initial conditions should  be saved to
    name: initial conditions will be written to files with the name name_####.in
    frmt: 'FISH' format as expected by Jens' FISH program, 'QMF3' format as expected by Matthias' QMF3 program
    """
    assert frmt in ["FISH", "QMF3"]
    Ncoord, Nsample = qs.shape
    for i in range(0, Nsample):
        fname = join(inidir, "%s_%.4d.in" % (name, i))
        if frmt == "FISH":
            txt = dynamics_in_format(atomlist, qs[:,i], ps[:,i], fname)
        elif frmt == "QMF3":
            txt = qmf3_in_format(atomlist, qs[:,i], ps[:,i], fname)
        fh = open(fname, "w")
        fh.write(txt)
        fh.close()
    print "Wrote initial conditions %s_%.4d.in to %s_%.4d.in to directory %s" % (name, 0, name, Nsample-1, inidir)

if __name__ == "__main__":
    pass

