#!/usr/bin/env python
"""
check for which exchange-correlation functionals we have the 
correct
   v_eff[rho] = -1/r    for r ---> oo
behaviour

CONCLUSION: All LDA functions have the wrong asymptotics.
  To describe photoionization with DFT we have to solve the
  electronic structure of the cation, for which the asymptotic
  effective potential is v_eff[rho^cation](r) --> -1/r as needed.
  The initial bound orbital is the LUMO of the cation while
  the final continuum orbital is the solution of

               __2           cation
        ( -1/2 \/  + v   [rho       ](r) - E ) phi = 0
                      eff
  
  using the same effective potential. In this way the initial and 
  final states belong to the same Hamiltonian which has the correct
  -1/r limit at long distance.


"""
from DFTB.MolecularIntegrals.BasissetFreeDFT import BasissetFreeDFT, density_func, effective_potential_func
from DFTB.MolecularIntegrals import settings
from DFTB.SlaterKoster import XCFunctionals

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Li^+ atom
    atomlist = [(3, (0.0, 0.0, 0.0))]
    charge = +1

    # choose resolution of multicenter grids for bound orbitals
    settings.radial_grid_factor = 20      # controls size of radial grid  
    settings.lebedev_order = 25          # controls size of angular grid
    # 1s core orbitals for Li+^ atom
    RDFT = BasissetFreeDFT(atomlist, None, charge=charge)
    bound_orbitals = RDFT.getOrbitalGuess()

    print "electron density..."
    # electron density of two electrons in the 1s core orbital
    rho = density_func(bound_orbitals)
    print "effective potential..."

    # List of (exchange, correlation) functionals implemented
    # in libXC
    functionals = [
        ("lda_x", "lda_c_xalpha"),
        ("lda_x_erf", "lda_c_xalpha"),
        ("lda_x_rae", "lda_c_xalpha"),
        ("gga_x_lb", "lda_c_xalpha"),
        ]

    r = np.linspace(0.1, 100.0, 1000)
    plt.xlabel(r"distance r / bohr")
    plt.ylabel(r"potential / Hartree")

    # correct asymptotic HF potential
    plt.plot(r, -2.0/r, "-.", lw=2, label=r"$-2/r$")
    
    for (x_func, c_func) in functionals:
        xc = XCFunctionals.libXCFunctional(x_func, c_func)
        # potential energy for Li nucleus and 2 core electrons
        potential = effective_potential_func(atomlist, rho, xc)
        
        plt.plot(r, potential(r,0*r,0*r), label=r"%s, %s" % (x_func, c_func))

    plt.ylim((-5.0, +0.1))
    plt.legend()
    plt.show()

    
