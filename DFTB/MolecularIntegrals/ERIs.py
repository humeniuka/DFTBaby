#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
compute electron repulsion integrals (ERI) between numerically exact atomic
basis functions using Becke's multicenter integration scheme

The most general ERI has the form

  (ab|1/r12 + f_xc[rho0]|cd)

where rho0 is the superposition of atomic densities of atoms which are individually neutral.
"""
from DFTB.BasisSets import AtomicBasisSet, import_pseudo_atom
from DFTB import XYZ
from DFTB.MolecularIntegrals.MulticenterIntegration import multicenter_poisson, multicenter_integration, atomlist2arrays
from DFTB.MolecularIntegrals import settings
from DFTB.SlaterKoster import XCFunctionals
from DFTB.SlaterKoster.SKIntegrals import spline_radial_density

import numpy as np
from scipy import interpolate

def electron_repulsion_integral(atomlist, bfA, bfB, bfC, bfD, rho, xc):
    """
    compute the electron repulsion integral

      (ab|{1/r12+f_xc[rho]}|cd)

    Parameters
    ==========
    atomlist            :  list of tuples (Z,[x,y,z]) with atomic numbers
                           and positions
    bfA, bfB, bfC, bfD  :  instances of AtomicBasisFunction or callables, e.g. bfA(x,y,z) etc.
    rho                 :  callable, rho(x,y,z) computes the total electron
                           density
    xc                  :  instance of libXCFunctional, only LDA functional will work properly
                           as the density gradient is not available

    Returns
    =======
    Iabcd               :  float, sum of Hartree and exchange correlation part of electron integral
    """
    
    def rhoAB(x,y,z):
        return bfA(x,y,z) * bfB(x,y,z)
    def rhoCD(x,y,z):
        return bfC(x,y,z) * bfD(x,y,z)
    
    Iabcd = electron_repulsion_integral_rho(atomlist, rhoAB, rhoCD, rho, xc)
    return Iabcd

def electron_repulsion_integral_rho(atomlist, rhoAB, rhoCD, rho, xc):
    """
    compute the electron repulsion integral

      (ab|{1/r12+f_xc[rho]}|cd)

    Parameters
    ==========
    atomlist            :  list of tuples (Z,[x,y,z]) with atomic numbers
                           and positions
    rhoAB, rhoCD        :  callables, rhoAB(x,y,z) computes the product a(x,y,z)*b(x,y,z)
                           and rhoCD(x,y,z) computes the product c(x,y,z)*d(x,y,z)
    rho                 :  callable, rho(x,y,z) computes the total electron
                           density
    xc                  :  instance of libXCFunctional, only LDA functional will work properly
                           as the density gradient is not available

    Returns
    =======
    Iabcd               :  float, sum of Hartree and exchange correlation part of electron integral
    """
    # bring data into a form understood by the module MolecularIntegrals
    atomic_numbers, atomic_coordinates = atomlist2arrays(atomlist)
    # Now we compute the integrals numerically on a multicenter grid.
    #
    # compute electrostatic Hartree term
    # (ab|1/r12|cd)
    # 1. solve the Poisson equation to get the electrostatic potential
    #    Vcd(r) due to the charge distribution c(r)*d(r)    
    Vcd = multicenter_poisson(rhoCD, atomic_coordinates, atomic_numbers,
                              radial_grid_factor=settings.radial_grid_factor,
                              lebedev_order=settings.lebedev_order)
    #
    # 2. integrate a(r)*b(r)*Vcd(r)
    def Iabcd_hartree_integrand(x,y,z):
        return rhoAB(x,y,z) * Vcd(x,y,z)

    # Coulomb integral 
    Iabcd_hartree = multicenter_integration(Iabcd_hartree_integrand, atomic_coordinates, atomic_numbers,
                                            radial_grid_factor=settings.radial_grid_factor,
                                            lebedev_order=settings.lebedev_order)

    #
    # compute contribution from exchange-correlation functional
    # (ab|f_xc[rho]|cd)
    def Iabcd_fxc_integrand(x,y,z):
        return rhoAB(x,y,z) * xc.fxc(rho(x,y,z)) * rhoCD(x,y,z)

    Iabcd_xc = multicenter_integration(Iabcd_fxc_integrand, atomic_coordinates, atomic_numbers,
                                       radial_grid_factor=settings.radial_grid_factor,
                                       lebedev_order=settings.lebedev_order)

    Iabcd = Iabcd_hartree + Iabcd_xc

    # check that density integrates to the correct number of electrons
    total_elec_charge = multicenter_integration(rho, atomic_coordinates, atomic_numbers,
                                           radial_grid_factor=settings.radial_grid_factor,
                                           lebedev_order=settings.lebedev_order)
    total_nuc_charge = sum([Zi for (Zi,posi) in atomlist])
    #print "total electronic charge :  %e" % total_elec_charge
    #print "total nuclear charge    :  %e" % total_nuc_charge
    assert abs(total_elec_charge - total_nuc_charge) < 1.0e-3

    #print "Hartree contribution (ab|1/r12|cd)      = %+e" % Iabcd_hartree
    #print "XC-contribution      (ab|f_xc[rho0]|cd) = %+e" % Iabcd_xc       
    
    return Iabcd

class AtomicDensitySuperposition:
    """
    superposition of individually neutral atomic densities

      rho0(r) = sum_A rho_A(r)
    """
    def __init__(self, atomlist, confined=True):
        atomtypes = list(set([Zi for (Zi,posi) in atomlist]))
        atomtypes.sort()
        self.atomlist = atomlist
        self.atomic_densities = {}
        for Zi in atomtypes:
            # load radial density of pseudo atoms
            confined_atom, free_atom = import_pseudo_atom(Zi)
            if confined == True:
                atom = confined_atom
            else:
                atom = free_atom
            # The density at r=0 is not correct, so leave it out
            r = atom.r
            rhoI_spline = spline_radial_density(atom.r[r > 0], atom.radial_density[r > 0])
            self.atomic_densities[Zi] = rhoI_spline
    def __call__(self, x,y,z):
        """ evaluate the superposition of atomic densities on a grid """
        rho0 = 0*x
        for Zi,posi in self.atomlist:
            rhoI_spline = self.atomic_densities[Zi]
            xI,yI,zI = x-posi[0], y-posi[1], z-posi[2]
            # distance to atomic center I
            rI = np.sqrt(xI**2+yI**2+zI**2)
            # add unperturbed density of atom I
            rho0 += rhoI_spline(rI)
        return rho0

# convert angular momentum quantum numbers (l,m) to spectroscopic notation
angmom_to_xyz = {(0,0): "s ",
                 (1,-1): "px", (1,1): "py", (1,0): "pz",
                 (2,-2): "dxy", (2,-1): "dyz", (2,0): "dz2", (2,1): "dzx", (2,2): "dx2y2"}

    
if __name__ == "__main__":
    # set accuracy of multicenter grid
    settings.radial_grid_factor = 3
    settings.lebedev_order = 23
    
    atomlist = [(6, (0,0,0))]
    basis = AtomicBasisSet(atomlist)
    density = AtomicDensitySuperposition(atomlist)
    xc_functional = XCFunctionals.libXCFunctional("lda_x", "lda_c_pw")
    
    for a,bfA in enumerate(basis.bfs):
        for b,bfB in enumerate(basis.bfs):
            for c,bfC in enumerate(basis.bfs):
                for d,bfD in enumerate(basis.bfs):
                    eri = electron_repulsion_integral(atomlist, bfA, bfB, bfC, bfD, density, xc_functional)
                    print "(%d,%d|{1/r12+f_xc[rho0]}|%d,%d)= %e" % (a,b,c,d,eri)
    
