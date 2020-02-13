#!/usr/bin/env python

from DFTB.MolecularIntegrals.BasissetFreeDFT import regular_coulomb_func, irregular_coulomb_func, effective_potential_func, orbital_transformation, residual_func
from DFTB.MolecularIntegrals.MulticenterIntegration import atomlist2arrays, multicenter_operation
from DFTB.MolecularIntegrals.Ints1e import integral, overlap, kinetic, nuclear
from DFTB.MolecularIntegrals import settings

import numpy as np
import numpy.linalg as la

def linear_combination(atomlist, orbitals, coeffs):
    norb = len(orbitals)
    assert len(orbitals) == len(coeffs)
    
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    lcao_orb = multicenter_operation(orbitals,
                              lambda orbs: sum( [coeffs[j] * orbs[j] for j in range(0, norb)] ),
                              atomic_coordinates, atomic_numbers,
                              radial_grid_factor=settings.radial_grid_factor,
                              lebedev_order=settings.lebedev_order)
    return lcao_orb
    
def hmi_variational_kohn():
    # H2^+
    # bond length in bohr
    R = 2.0
    # positions of protons
    posH1 = (0.0, 0.0, -R/2.0)
    posH2 = (0.0, 0.0, +R/2.0)
    # center of charge
    center = (0.0, 0.0, 0.0)

    
    atomlist = [(1, posH1),
                (1, posH2)]
    charge = +2
    
    # choose resolution of multicenter grids for continuum orbitals
    settings.radial_grid_factor = 10      # controls size of radial grid  
    settings.lebedev_order = 21          # controls size of angular grid

    # energy of continuum orbital (in Hartree)
    E = 10.0
    # angular momentum quantum numbers
    l,m = 0,0

    rho = None
    xc = None
    print "effective potential..."
    V = effective_potential_func(atomlist, rho, xc, nelec=0)

    # basis functions
    u1 = regular_coulomb_func(E, +1, l, m, 0.0,
                              center=posH1)
    eta1 = 0.0
    u2 = regular_coulomb_func(E, +1, l, m, 0.0,
                              center=posH2)
    eta2 = 0.0
    
    # asymptotic solutions
    s0 = regular_coulomb_func(E, charge, l, m, 0.0,
                             center=center)
#    c0 = irregular_coulomb_func(E, charge, l, m, 0.0,
#                             center=center)
    c0 = regular_coulomb_func(E, charge, l, m, np.pi/2.0,
                             center=center)


    # switching function
    #  g(0) = 0    g(oo) = 1
    def g(r):
        r0 = 3.0  # radius at which the switching happens
        gr = 0*r
        gr[r > 0] = 1.0/(1.0+np.exp(-(1-r0/r[r > 0])))
        return gr

    def s(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        return g(r) * s0(x,y,z)

    def c(x,y,z):
        r = np.sqrt(x*x+y*y+z*z)
        return g(r) * c0(x,y,z)

    free = [s,c]
    
    # basis potentials
    def V1(x,y,z):
        r2 = (x-posH1[0])**2 + (y-posH1[1])**2 + (z-posH1[2])**2
        return -1.0/np.sqrt(r2)

    def V2(x,y,z):
        r2 = (x-posH2[0])**2 + (y-posH2[1])**2 + (z-posH2[2])**2
        return -1.0/np.sqrt(r2)

    basis = [u1,u2]
    basis_potentials = [V1,V2]

    # number of basis functions
    nb = len(basis)
    
    # matrix elements between basis functions
    Lbb = np.zeros((nb,nb))
    for i in range(0, nb):
        ui = basis[i]
        for j in range(0, nb):
            print "i= %d  j= %d" % (i,j)
            uj = basis[j]
            Vj = basis_potentials[j]
            def integrand(x,y,z):
                return ui(x,y,z) * (V(x,y,z) - Vj(x,y,z)) * uj(x,y,z)
            Lbb[i,j] = integral(atomlist, integrand)

    print Lbb

    # matrix elements between unbound functions
    Lff = np.zeros((2,2))

    for i in range(0, 2):
        fi = free[i]
        for j in range(0, 2):
            print "i= %d  j= %d" % (i,j)
            fj = free[j]
            Lff[i,j] = kinetic(atomlist, fi,fj) + nuclear(atomlist, fi,fj) - E * overlap(atomlist, fi,fj)

    print Lff

    # matrix elements between bound and unbound functions
    Lbf = np.zeros((nb,2))
    Lfb = np.zeros((2,nb))

    for i in range(0, nb):
        ui = basis[i]
        for j in range(0, 2):
            fj = free[j]
            print "i= %d  j= %d" % (i,j)
            Lbf[i,j] = kinetic(atomlist, ui,fj) + nuclear(atomlist, ui,fj) - E * overlap(atomlist, ui,fj)
            Lfb[j,i] = kinetic(atomlist, fj,ui) + nuclear(atomlist, fj,ui) - E * overlap(atomlist, fj,ui)

    print Lbf
    print Lfb

    #
    coeffs = -np.dot(la.inv(Lbb), Lbf)
    
    M = Lff - np.dot(Lfb, np.dot(la.inv(Lbb), Lbf))

    print "M"
    print M

    # The equation
    #    
    #  M. (1) = 0
    #     (t)
    # cannot be fulfilled, usually.
    
    t0 = -M[0,0]/M[0,1]

    tans = np.array([np.tan(eta1), np.tan(eta2)])
    t = (t0 + np.sum(coeffs[:,0]*tans) + t0*np.sum(coeffs[:,1]*tans)) \
         / (1 + np.sum(coeffs[:,0]) + t0*np.sum(coeffs[:,1]))
    eta = np.arctan(t)

    # Because a global phase is irrelevant, the phase shift is only determined
    # modulo pi. sin(pi+eta) = -sin(eta)
    while eta < 0.0:
        eta += np.pi

    print "eta= %s" % eta

    bc0 = linear_combination(atomlist, basis, coeffs[:,0])
    bc1 = linear_combination(atomlist, basis, coeffs[:,1])
    #u = s + np.dot(basis, coeffs[:,0]) + t0 * (c + np.dot(basis, coeffs[:,1]))
    orbitals = [s, bc0, c, bc1]
    nrm = (1 + np.sum(coeffs[:,0]) + t0*np.sum(coeffs[:,1]))
    orb_coeffs = np.array([ 1.0, 1.0, t0, t0 ]) / nrm

    phi = linear_combination(atomlist, orbitals, orb_coeffs)

    residual = residual_func(atomlist, phi, V, E)

    import matplotlib.pyplot as plt
    r = np.linspace(-50.0, 50.0, 5000)

    x = 0.0*r
    y = 0.0*r
    z = r

    plt.plot(r, phi(x,y,z), label=r"$\phi$")
    plt.plot(r, residual(x,y,z), label=r"residual $(H-E)\phi$")
    plt.plot(r, V(x,y,z), label=r"potential")

    phi_lcao = linear_combination(atomlist, basis, [1.0/np.sqrt(2.0), 1.0/np.sqrt(2.0)])
    residual_lcao = residual_func(atomlist, phi_lcao, V, E)

    plt.plot(r, phi_lcao(x,y,z), label=r"$\phi_{LCAO}$")
    plt.plot(r, residual_lcao(x,y,z), label=r"residual $(H-E)\phi_{LCAO}$")

    
    plt.legend()
    plt.show()
    
    
if __name__ == "__main__":
    hmi_variational_kohn()
