"""
Stability of closed-shell ground state wavefunction

The closed-shell restricted wavefunction is tested for restricted/unrestricted instabilities.
The stability condition is that the matrix

  3         3
   A      +  B      =  (e  - e ) delta   delta    -  (ij|ab)  - (ib|aj)
     ia,jb    ia,jb      a    i       ij      ab                       lr

is positive definite (i.e. has no negative eigenvalue)


References
----------
[1] R. Reeger, J. Pople
    "Self-consistent molecular orbital methods. XVIII. Constraints and stability in Hartree-Fock theory",
    J. Chem. Phys. 66, 3045 (1977)
"""
import numpy as np
import numpy.linalg as la
import scipy.sparse.linalg as sla

def check_stability(gamma, gamma_lr,
                    qtrans_oo, qtrans_vv, qtrans_ov,
                    omega, nocc, nvirt):
    print "Stability of closed-shell DFTB ground state"

    # Coulomb integrals
    # (ij|ab)
    ijab = np.tensordot(qtrans_oo, np.tensordot(gamma, qtrans_vv, axes=(1,0)), axes=(0,0))
    ijab = np.swapaxes(ijab, 1,2)
    # exchange integrals
    # (ia|jb)_lr
    iajb_lr = np.tensordot(qtrans_ov, np.tensordot(gamma_lr, qtrans_ov, axes=(1,0)), axes=(0,0))

    # double indices
    #  s = (i->a)     t = (j->b)
    shape = (nocc*nvirt, nocc*nvirt)
    ijab    = np.reshape(ijab,       shape)
    iajb_lr = np.reshape(iajb_lr,    shape)

    omega = np.reshape(omega, nocc*nvirt)
    
    # restricted/unrestricted stability matrix
    R = np.diag(omega) - ijab - iajb_lr

    # The restricted ground state is stable if the R matrix is positive semidefinite
    w, v = sla.eigsh(R, which='SA', k=min(6, nocc*nvirt))
    print w
    print "  lowest eigenvalue of restricted/unrestricted stability matrix R: %e" % w[0]
    if (w[0] < 0.0):
        print "   !!! restricted/unrestricted instability detected !!!"


    
