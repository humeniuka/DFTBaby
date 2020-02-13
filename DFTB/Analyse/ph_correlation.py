#!/usr/bin/env python
"""

  electron-hole correlation function

    g_eh(r) = int rho_e(r') rho_h(r'+r) dr'

  The function is integrated over the orientation so that it only
  depends on the |r|

"""
import numpy as np
import numpy.linalg as la

from DFTB import Parameters, XYZ

def G_matrix(atomlist, Nr):
    Nat = len(atomlist)
    # sigmas
    # find unique atom types
    atomtypes = list(set([Zi for (Zi,posi) in atomlist]))
    atomtypes.sort()

    print "Computing sigmas..."
    sigma = np.zeros(max(atomtypes)+1)
    for Zi in atomtypes:
        Ui = Parameters.hubbard_U_byZ[Zi]
        sigma[Zi] = 1.0/(np.sqrt(np.pi) * Ui)

    Nat = len(atomlist)
    print "Computing distance matrix..."
    Dmat = np.zeros((Nat,Nat))
    for A,(ZA,posA) in enumerate(atomlist):
        for B in range(A+1, Nat):
            ZB, posB = atomlist[B]
            R_AB = la.norm(np.array(posA) - np.array(posB))
            Dmat[A,B] = R_AB
            Dmat[B,A] = R_AB
    rs = np.linspace(0.0, Dmat.max(),Nr)
    print "Computing G-matrix ..."
    # G-matrix
    G = np.zeros((Nat,Nat,Nr))
    for A,(ZA,posA) in enumerate(atomlist):
        for B in range(A, Nat):
            ZB, posB = atomlist[B]
            R_AB = Dmat[A,B]
            s2_AB = sigma[ZA]**2 + sigma[ZB]**2
            f_AB = 0.5/s2_AB
            G[A,B,:] = np.exp(-f_AB * (R_AB - rs)**2)\
                      -np.exp(-f_AB * (R_AB + rs)**2)           
            G[A,B,:] /= np.sqrt(2*np.pi*s2_AB)
            G[B,A,:] = G[A,B,:]
    return rs, G

def electron_hole_correlation(G, q_elec, q_hole):
    Nat,Nat,Nr = G.shape
    assert len(q_elec) == Nat and len(q_hole) == Nat, "Number of particle/hole charges (%d) does not equal number of atoms (%d)!" % (len(q_elec, Nat))
    g_eh = np.zeros(Nr)
    for ir in range(0, Nr):
        g_eh[ir] = np.dot(q_elec, np.dot(G[:,:,ir], q_hole))
    return g_eh

if __name__ == "__main__":
    import sys
    from os.path import expandvars, expanduser
    from optparse import OptionParser
    usage = "Usage: python %s <xyz file> <files with p-h charges>\n\n" % sys.argv[0]
    usage +="  compute the electron-hole correlation function g_eh(r), which indicates\n"
    usage +="  the likelyhood of finding a positive charge at a distance r from\n"
    usage +="  a negative charge:\n"
    usage += "    g_eh(r) = int rho_e(r') rho_h(r'+r) dr' dOmega\n"
    parser = OptionParser(usage)
    parser.add_option("--save_png", dest="save_png", help="save plot to this file [default: %default]", default="")
    parser.add_option("--Nr", dest="Nr", help="Compute e-h correlation function on a grid of Nr points in the interval r=[0,max atom separation] [default: %default]", default=100, type=int)

    (opts, args) = parser.parse_args()
    if len(args) < 2:
        print usage
        exit(-1)

    xyz_file = args[0]
    ph_charges_files = args[1:]

    atomlist = XYZ.read_xyz(xyz_file)[0]

    rs, Gmat = G_matrix(atomlist, opts.Nr)

    from matplotlib import pyplot as plt
    plt.xlabel("r / bohr")
    plt.ylabel("g$_{e-h}$(r)")
    plt.title("Electron hole correlation function")

    for ph_f in ph_charges_files:
        print "computing p-h correlation function for %s" % ph_f
        ph_charges = np.loadtxt(ph_f)
    # particle charges
        q_elec = ph_charges[:,0]
        q_hole = ph_charges[:,1]
    
        g_eh = electron_hole_correlation(Gmat, q_elec, q_hole)

        plt.plot(rs, g_eh, lw=2)
    plt.show()
