#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

This module implements Johnson's *multiple scattering* (MS) method for bound orbitals [1]_
and Dill and Dehmer's *continuum multiple scattering* (CMS) method [2]_. 


.. rubric:: References
.. [1] K.H. Johnson,  "Scattered-Wave Theory of the Chemical Bond" (1973), Advances in quantum chemistry, volume 7, pages 143-185, Elsevier
.. [2] D.Dill, J.Dehmer, "Electron-molecule scattering and molecular photoionization using the multiple-scattering method", J. Chem. Phys. 61, 692 (1974)
.. [3] M.Danos, L.Maximon, "Multipole Matrix Elements of the Translation Operator", Journal of Mathematical Physics 6, 766 (1965)
.. [4] D.Loomba, S.Wallace, D.Dill, J.Dehmer, "Pictures of unbound molecular electrons, including shape-resonant states. Eigenchannel contour maps", J. Chem. Phys. 75, 4546 (1981)

"""
from __future__ import print_function

from DFTB import AtomicData
from DFTB.AtomicData import atom_names, slater_radii
from DFTB.Analyse import Cube

from DFTB.MolecularIntegrals import settings
from DFTB.MolecularIntegrals.MulticenterIntegration import spherical_average_func, multicenter_grids, join_grids
from DFTB.MolecularIntegrals.BasissetFreeDFT import effective_potential_func, density_func, total_charge, overlap

from DFTB.SlaterKoster import XCFunctionals
from DFTB.SlaterKoster.free_pseudo_atoms import pseudo_atoms_list

from DFTB.Timer import GlobalTimer as T

from DFTB.MultipleScattering.Wavefunctions import MSWavefunction, CMSWavefunction

# Most parts of the MS method are implemented in Fortran in the module `ms.so`.
# The python code just acts as a wrapper. 
from DFTB.MultipleScattering import ms
# for evaluating the muffing tin potential
from DFTB.MultipleScattering import cms
# for Numerov integration in regions I and III
from DFTB.MultipleScattering import numerov
# spherical harmonics
from DFTB.MultipleScattering.sphharm import sphharm, enumerate_lm
# Photoionization cross sections are computed in `photo.so`
from DFTB.MultipleScattering import photo

import numpy as np
import numpy.linalg as la
from scipy import interpolate
from scipy import signal
import mpmath # for gammainc

def close_packed_spheres(atomlist, debug=0):
    """
    In the muffin tin method space is partitioned into three spherical regions:

    -  region I    :  sphere of radius `radii[i]` around each atom `i`
    -  region III  :  sphere of radius `rhoIII` around the origin, which encloses all other spheres
    -  region II   :  interstitial region between I and III

    This function chooses the radii `rhoIII` and `radii[i]` (i=1,...,Nat) based
    on atomic Slater radii, so as to obtain a relatively closely packed assembly. 

    Parameters
    ----------
    atomlist     :  list of tuple
                    list of tuples `(Zat,[xi,yi,zi])` with atomic number
                    Zat and nuclear positions `[xi,yi,zi]` for each atom
    debug        :  int, optional
                    if set to 1, the spheres around regions I and III
                    are plotted in the yz-plane

    Returns
    -------
    rhoIII       :  float
                    radius of sphere around origin separating region II from region III
    radii        :  np.1darray (shape: (Nat,))
                    radii of region I spheres around each atom

    """
    # Find the radii for regions I
    Nat = len(atomlist)
    assert Nat > 1
    
    # None of the atomic centers should lie exactly at the origin,
    # since in this case in N_coeffs(...) one has to evaluate the
    # spherical Bessel functions of the second kind n_l(r) at r=0
    # where n_l(r=0)=inf.
    for i,(Zi,posi) in enumerate(atomlist):
        assert la.norm(posi) > 0.0, "No atom should lie exactly at the origin, \n but atom %d lies at %s !" % (i, str(posi))   
    
    atomic_names = [atom_names[Z-1] for Z,pos in atomlist]
    # A sphere is put around each atom which does not penetrate
    # any other sphere.
    # The distance between two atoms is divided according to the
    # Slater radii.
    radii = np.zeros(Nat)
    for i,(Zi,posi) in enumerate(atomlist):
        radius_i = np.inf
        for j,(Zj,posj) in enumerate(atomlist):
            if i == j:
                continue
            # distance between atoms i and j
            Rij = la.norm(np.array(posi) - np.array(posj))
            # Slater radii
            rsi = slater_radii[atomic_names[i]]
            rsj = slater_radii[atomic_names[j]]
            # s*Rij belongs to atom i, (1-s)*Rij to atom j
            s =  rsi/(rsi+rsj)
            radius_i = min(radius_i, s*Rij)

        radii[i] = radius_i

    # Radius for region III
    rhoIII = 0.0
    for i,(Zi,posi) in enumerate(atomlist):
        rhoIII = max(rhoIII, la.norm(np.array(posi)) + radii[i])

    print("radius of molecular sphere (in bohr)")
    print(" rhoIII= %7.5f" % rhoIII)
    print("radii of atomic spheres (in bohr):")
    for i in range(0, Nat):
        print("  %s-%d    %7.5f" % (atomic_names[i].upper(),i+1, radii[i]))
    print("")
    
    if debug:
        # show partitioning in yz-plane
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.set_aspect(aspect='equal')
        ax.set_axis_off()
        
        ax.set_xlim((-5, 5))
        ax.set_ylim((-5, 5))
        
        # plot circle that separates region II from region III
        regionIII = plt.Circle((0.0, 0,0), rhoIII, color='black', fill=False)
        ax.add_artist(regionIII)

        # plot circles that separate regions I from regions II
        for i,(Zi, posi) in enumerate(atomlist):
            regionII = plt.Circle((posi[1], posi[2]), radii[i], color='blue', fill=False)
            ax.add_artist(regionII)

        plt.show()
        
    return rhoIII, radii

def atomic_superposition_density(atomlist):
    """

    approximate electronic density of a molecule as a superposition of the
    densities of the free atoms

    Parameters
    ----------
    atomlist     :  list of tuple
                    list of tuples `(Zat,(x,y,z))` with the atomic
                    positions that define the multicenter grid

    Returns
    -------
    rho          :  callable `rho(x,y,z)`
                    function for evaluating the superposition density on a grid
                    
    """
    # unique atom types 
    atomtypes = list(set([Zi for (Zi,posi) in atomlist]))
    atomtypes.sort()
    # interpolate atomic radial densities
    atomic_densities = {}
    for Zi in atomtypes:
        # load pseudo atom data
        at = pseudo_atoms_list[Zi-1]
        Nelec = at.Z

        # Remove the r=0 point
        r = at.r[1:]
        rho = at.radial_density[1:]
        # To ensure that extrapolation of rho(r) beyond r.max() gives 0,
        # we explicitly set rho to 0 for the last 4 points.
        rho[-4:] = 0.0
        radial_density = interpolate.interp1d(r, rho,
                                              fill_value="extrapolate")
        atomic_densities[Zi] = radial_density

    def rho(x,y,z):
        dens = 0.0*x
        
        for Zi,(xi,yi,zi) in atomlist:
            # distance from atom i
            ri = np.sqrt( (x-xi)**2 + (y-yi)**2 + (z-zi)**2 )
            # radial density of atom i
            dens_i = atomic_densities[Zi](ri)
            
            dens += dens_i
            
        return dens
            
    return rho

def save_cubefile(atomlist, wfn,
                  filename="/tmp/wavefunction_",
                  ppb=2.0,
                  dbuff=5.0):
    """
    save a wavefunction to a cube file, real and imaginary parts are saved separately

    Parameters
    ----------
    atomlist   : list of tuple
                 list of tuples `(Zat,(x,y,z))` with the atomic
                 positions that define the multicenter grid
    wfn        : callable
                 `wfn(x,y,z)` should evaluate the wavefunction on a grid, 
                 where `x`, `y` and `z` are 1d numpy arrays with the cartesian positions 
                 of the grid points.
    filename   : str, optional
                 path to cube file, the suffixes '_real.cube' and '_imag.cube' are appended 
                 to the filename
    ppb        : float, optional
                 resolution of cube file, points per bohr
    dbuff      : float, optional
                 extra space around molecule in bohr
    
    """
    def real_part(grid,dV):
        x,y,z = grid
        xx,yy,zz = x.ravel(), y.ravel(), z.ravel()
        wfn_xyz = wfn(xx,yy,zz).real
        wfn_grid = np.reshape(wfn_xyz, x.shape)
        return wfn_grid

    def imag_part(grid,dV):
        x,y,z = grid
        xx,yy,zz = x.ravel(), y.ravel(), z.ravel()
        wfn_xyz = wfn(xx,yy,zz).imag
        wfn_grid = np.reshape(wfn_xyz, x.shape)
        return wfn_grid
    
    Cube.function_to_cubefile(atomlist, real_part,
                              filename=filename+"_real.cube", ppb=ppb, dbuff=dbuff)
    Cube.function_to_cubefile(atomlist, imag_part,
                              filename=filename+"_imag.cube", ppb=ppb, dbuff=dbuff)


class Regions(object):
    """
    class for storing potentials of regions I,II and III and for computing
    radial wavefunctions in regions I and III

    Parameters
    ----------
    atomlist     :  list of tuple
                    list of tuples `(Zat,(x,y,z))` with the atomic
                    positions that define the multicenter grid
    charge       :  int, optional
                    total charge of molecule, (total nuclear charge) - (number of electrons)
    chargeIII    :  int, optional
                    non-screened charge felt by an electron in region III, usually `charge+1`
    lmax         :  int, optional
                    highest angular momentum at which the spherical wave expansion is truncated
    NrI          :  int, optional
                    number of radial grid points for interpolating potentials
                    and solving the radial Schroedinger equation in region I
    NrIII        :  int, optional
                    number of radial grid points in region III

    """

    def __init__(self,
                 atomlist, 
                 charge=0, chargeIII=+1,
                 lmax=8, NrI=5000, NrIII=10000):
        self.atomlist = atomlist
        self.nat = len(atomlist)
        # total charge of molecule (not equal chargeIII)
        self.charge = charge
        # nuclear charge felt by an electron in region III
        self.chargeIII = chargeIII
        self.lmax = lmax
        
        # unique atom types 
        self.atomtypes = list(set([Zi for (Zi,posi) in atomlist]))
        self.atomtypes.sort()

        # count total number of electrons
        self.nelec = 0
        for i,(Zi,posi) in enumerate(atomlist):
            self.nelec += Zi
        self.nelec -= self.charge
        assert self.nelec > 0, "number of electrons must be positive"

        print("total molecular charge             : %d" % self.charge)
        print("non-screened charge in region III  : %d" % self.chargeIII)
        print("number of electrons                : %d" % self.nelec)

        # partition space into regions I,II and III
        self.rhoIII, self.radii = close_packed_spheres(atomlist)
        # radii of atomic spheres for each atom type
        self.radii_Z = {}
        for i,(Zi,posi) in enumerate(atomlist):
            self.radii_Z[Zi] = self.radii[i]

        # find number of knots of B-splines in region I
        self.NrI = NrI
        r = np.linspace(1.e-3, 1.0, self.NrI)
        t,c,k = interpolate.splrep(r,r)
        self.nsplI = len(t)
        self.kspl = k

        # find number of knots of B-splines in region III
        self.NrIII = NrIII
        r = np.linspace(1.e-3, 1.0, self.NrIII)
        t,c,k = interpolate.splrep(r,r)
        self.nsplIII = len(t)
        assert k == self.kspl

        # initialize empty arrays for potentials in region I
        
        # functions for evaluating potentials V^I_i(r) for each atom
        self.pot_funcs = [None for i in range(0, self.nat)]
        # assign B-splines for radial potentials (knots and coefficients)
        self.pot_knots  = np.zeros((self.nsplI, self.nat))
        self.pot_coeffs = np.zeros((self.nsplI, self.nat))
        # radial grids and radial potentials at the grid points 
        self.pot_grids = [None for i in range(0, self.nat)]
        self.pot_vals  = [None for i in range(0, self.nat)]

        # initialize empty arrays for wavefunctions in regions I and III
        
        # knots of radial wavefunctions in region I
        self.fIknots  = np.zeros((self.nsplI, self.lmax+1, self.nat))
        # and coefficients
        self.fIcoeffs = np.zeros((self.nsplI, self.lmax+1, self.nat))
        
        # knots and coefficients of wavefunctions in region III
        self.fIIIknots  = np.zeros((self.nsplIII, self.lmax+1))
        self.fIIIcoeffs = np.zeros((self.nsplIII, self.lmax+1))

        # constant potential in region II
        self.vconst = None
        
    ###################################################
    #                                                 #
    #      radial potential V^I  and V^II (const)     #
    #                                                 #
    ###################################################

    def get_potential(self, potential=None):
        """
        determine the spherically symmetric potentials in each atomic sphere.
        If a molecular potential is provided, it is spherically averaged around 
        each atom to obtain atomic potentials that partially account for the 
        molecular environment. Otherwise the region I potentials are assumed to
        be identical to the DFT potentials of free atoms.

        Parameters
        ----------
        potential     :   callable, optional
                          `potential(x,y,z)` should evaluate the molecular potential on a grid

        """
        if potential is None:
            self.potential_type = "atomic"
            self._load_regionI_potentials()
        else:
            self.potential_type = "molecular"
            self.potential = potential
            self._average_regionI_potentials(potential)

        self.vconst = self._average_regionII_potential(self.pot_funcs)
    
    def _average_regionI_potentials(self, potential):
        """
        average molecular potential spherically around each atom

        .. code-block:: none
        
            (sph)               /
           V       =   1/(4 pi) | V(r-R ) dOmega
            i                   /      i        i

        Parameters
        ----------
        potential    :  callable
                        potential(x,y,z) evaluates the potential V on on a grid

        """
        print("average molecular potential spherically around each atom")
                
        for i,(Zi,posi) in enumerate(self.atomlist):
            atom_center = (Zi, posi)
            # average molecular potential spherically around atom I
            potI = spherical_average_func(atom_center, potential,
                                          radial_grid_factor=settings.radial_grid_factor,
                                          lebedev_order=settings.lebedev_order)
            
            # define radial grid for integrating wavefunctions
            rmin = 1.e-3
            # The wavefunction is only needed inside the atomic sphere
            # of radius rho_i
            rmax = 1.2 * self.radii_Z[Zi]
            r = np.linspace(rmin, rmax, self.NrI)
            # evaluate radial potential at the grid points
            pot_r = potI(r)

            # create B-spline
            t,c,k = interpolate.splrep(r, pot_r)

            # save potential for atom i
            self.pot_funcs[i] = potI
            
            self.pot_knots[:,i]  = t
            self.pot_coeffs[:,i] = c
            
            self.pot_grids[i] = r
            self.pot_vals[i]  = pot_r

    def _load_regionI_potentials(self):
        """
        load precalculated atomic potentials for region I

        The potentials in region I are taken to be the atomic DFT potentials.
        :math:`V^I_i(r)` is the self-consistent DFT potential of atom i. The contribution
        from neighbouring atoms to the potential in the atomic sphere i is completely
        neglected.

        """
        # functions for evaluating radial potentials (by atom type)
        pot_funcs_Z = {}
        
        # B-splines of radial potentials (by atom type)
        pot_knots_Z = {}
        pot_coeffs_Z = {}

        # radial grids and potentials on the grid (by atom type)
        pot_grids_Z = {}
        pot_vals_Z = {}

        # First we compute the radial potentials for each atom type. In large
        # molecules, where many atoms are of the same type, this means a huge
        # saving.
        for Zi in self.atomtypes:
            print("load atomic potential for atom type %s" % Zi)
            # load pseudo atom data
            at = pseudo_atoms_list[Zi-1]
            Nelec = at.Z

            # interpolate radial potential
            if self.nelec == 1:
                # If there is only a single electron, the electron-electron repulsion
                # and exchange-correlation part of the atomic potential is not needed,
                # therefore we replace the atomic DFT potential by the bare nuclear
                # attraction, -Z/r
                pot = interpolate.interp1d(at.r[1:], -at.Z/at.r[1:],
                                           fill_value="extrapolate")
                print("   using bare nuclear potential %d/r" % (-at.Z))
            else:
                # load atomic DFT potentials (nuclear attraction + Coulomb repulsion + xc)
                # The first point r=0 is excluded because there the potential
                # diverges.
                pot = interpolate.interp1d(at.r[1:], at.effective_potential[1:],
                                           fill_value="extrapolate")
            pot_funcs_Z[Zi] = pot
            
            # define radial grid for integrating wavefunctions
            rmin = at.r[1]
            # The wavefunction is only needed inside the atomic sphere
            # of radius rho_i
            rmax = 1.2 * self.radii_Z[Zi]
            r = np.linspace(rmin, rmax, self.NrI)
            pot_r = pot(r)
            
            pot_grids_Z[Zi] = r
            pot_vals_Z[Zi] = pot_r

            # create B-spline
            t,c,k = interpolate.splrep(r, pot_r)

            pot_knots_Z[Zi] = t
            pot_coeffs_Z[Zi] = c

        # Next we assign potentials (and their B-spline representation) to atoms,
        # atoms of the same atom type get the same potential irrespective of the
        # molecular environment.
        for i,(Zi,posi) in enumerate(self.atomlist):
            self.pot_funcs[i]    = pot_funcs_Z[Zi]
            
            self.pot_knots[:,i]  = pot_knots_Z[Zi]
            self.pot_coeffs[:,i] = pot_coeffs_Z[Zi]

            self.pot_grids[i] = pot_grids_Z[Zi]
            self.pot_vals[i]  = pot_vals_Z[Zi]
            
    def _average_regionII_potential(self, potentialsI):
        """
        The constant potential in region II (vconst) is obtained by averaging
        the potentials of region I and III over the molecular and atomic spheres

        The potential on the atomic and molecular spheres is weighted by the 
        area of the sphere:

        .. code-block:: none

                         Nat         I                       III                  Nat 
           vconst = [ sum    Area   V (rho )    +   Area    V   (rho   )  ] / (sum     Area   +  Area   )
                         i       i   i    i             III         III           i        i         III

           with 
                           2
           Area  = 4 pi rho 
               i           i                   
          
        In region III, we have a pure Coulomb potential
        
        .. code-block:: none

            III          charge^III
           V   (r)  =  - ----------
                             r

        Parameters
        ----------
        potentialsI  :  list of callables
                        which evaluate the radial potentials around each atom in region I,
                        `potentialsI[i](r)` should give :math:`V^I_i(r)` for the i-th atom

        Returns
        -------
        vconst       :  float
                        constant potential in region II

        """        
        pot_avg = 0.0
        area_tot = 0.0
        # atomic spheres
        for i in range(0, self.nat):
            rho_i = self.radii[i]
            # area of atomic sphere i
            area_i = 4.0 * np.pi * rho_i**2
            # potential I on the sphere
            potI_i = potentialsI[i](rho_i)
            # weight potential by area of atomic sphere
            pot_avg += area_i * potI_i
            area_tot += area_i

        # molecular sphere
        potIII = -self.chargeIII / self.rhoIII
        areaIII = 4.0 * np.pi * self.rhoIII**2

        pot_avg += areaIII * potIII
        area_tot += areaIII

        # average over area
        pot_avg /= area_tot

        vconst = pot_avg
        
        print("constant potential in region II (in Hartree)")
        print(" Vconst= %e" % vconst)

        self.vconst = vconst
        
        return vconst

    ###################################################
    #                                                 #
    #   radial wavefunctions in regions I and III     #
    #                                                 #
    ###################################################

    def _solve_radial_schroedingerI(self, r, vr, energy, rho_i):
        """
        solve the radial Schroedinger equation in an atomic sphare (region I)
        for given energy E and create B-spline for each angular momentum.

        The radial Schroedinger equation

        .. code-block:: none

                d^2                        l*(l+1)
          [ 1/2 ----  +  E  -  V(r)  - 1/2 ------- ] u (r)  =  0
                dr^2                         r^2      l

        is solved for the radial wavefunction :math:`u_l(r) = r R_l(r)`
        on an equidistant grid using Numerov's method.

        Parameters
        ----------
        r           :  np.1darray (shape: (Nr,))
                       array of radial equidistant grid points 
        vr          :  np.1darray (shape: (Nr,))
                       array with radial potential at grid points, `vr[i] = V(r[i])`
        energy      :  float
                       energy E in Hartree
        rho_i       :  float
                       radius of atomic sphere

        Returns
        -------
        fknots      :  np.2darray (shape: (nsplI,lmax+1))
                       knots of B-splines for :math:`R_l(r)`
        fcoeffs     :  np.2darray (shape: (nsplI,lmax+1)) 
                       coefficients of B-spline for :math:`R_l(r)`

        """
        fknots = np.zeros((self.nsplI, self.lmax+1))
        fcoeffs = np.zeros((self.nsplI, self.lmax+1))
        
        # integrate outward, all wavefunctions with l=0,1,...,lmax
        # are generated at once
        #  u[:,l] = u_l(r)
        u = numerov.numerov_out(r,vr, energy, self.lmax)

        # equidistant grid
        dr = np.ediff1d(r[r<=rho_i], to_end=r[2]-r[1])
        for l in range(0, self.lmax+1):
            # normalize u(r) such that
            #   /rho_i           2
            #   |       dr |u(r)|  = 1
            #   /0
            u2int = np.sum(abs(u[r<=rho_i,l])**2 * dr)
            u[:,l] /= np.sqrt(u2int)
        
            # create B-spline for R_l(r) = u_l(r)
            t,c,k = interpolate.splrep(r, u[:,l]/r)
            assert k == self.kspl
            fknots[:,l] = t
            fcoeffs[:,l] = c
            
        return fknots, fcoeffs

    def solve_wavefunctionsI(self, energy):
        """
        find atomic wavefunctions in region I for given energy

        Parameters
        ----------
        energy       :  float
                        energy in Hartree

        """
        if self.potential_type == "molecular":
            # radial wavefunctions are solved for each atom, as the radial potential
            # depends on the molecular environment
            for i,(Zi,posi) in enumerate(self.atomlist):
                # radial grid r and V(r)
                r = self.pot_grids[i]
                vr = self.pot_vals[i]

                # solve for R_l(r) and save wavefunctions in B-spline representation
                fknots, fcoeffs = self._solve_radial_schroedingerI(r, vr, energy, self.radii[i])
                self.fIknots[:,:,i]  = fknots
                self.fIcoeffs[:,:,i] = fcoeffs
                
        elif self.potential_type == "atomic":
            # wavefunctions for atoms of the same type are assumed to be the same

            # knots and coefficients of B-splines for atomic radial wavefunctions in region I
            fIknots_Z = {}
            fIcoeffs_Z = {}
            
            # Radial wavefunctions are computed for each atom type. Since
            # in a molecule many atoms are of the same type, we save a
            # a lot of time by computing the wavefunctions only once
            # per atom type.
            for Zi in self.atomtypes:
                r = self.pot_grids_Z[Zi]
                vr = self.pot_vals_Z[Zi]

                # solve for R_l(r) and save wavefunctions in B-spline representation
                fknots, fcoeffs = self._solve_radial_schroedingerI(r, vr, energy, self.radii[i])
                self.fIknots_Z[Zi]  = fknots
                self.fIcoeffs_Z[Zi] = fcoeffs

            # assign wavefunctions to atoms based on atom type
            for i,(Zi,posi) in enumerate(self.atomlist):
                self.fIknots[:,:,i]  = self.fIknots_Z[Zi]
                self.fIcoeffs[:,:,i] = self.fIcoeffs_Z[Zi]
                
    def solve_wavefunctionsIII(self, energy):
        """
        solve the Schroedinger equation of the Coulomb potential for given energy E

        .. code-block:: none

                __2           chargeIII
           -1/2 \/   phi + (- ---------  - E) phi = 0
                                  r

        in region III.

        Parameters
        ----------
        energy       :  float
                        energy in Hartree

        """
        # Wavefunctions are needed in the region rhoIII < r < oo.
        # rmax should be large enough that the molecular potential has decayed to 0.
        rmin = 0.8 * self.rhoIII
        rmax = self.rhoIII + 30.0
        # uniform radial grid
        r = np.linspace(rmin, rmax, self.NrIII)
        # Coulomb potential
        vr = -self.chargeIII/r
        
        # integrate inward, all wavefunctions with l=0,1,...,lmax
        # are generated
        #  u[:,l] = u_l(r)
        u = numerov.numerov_in(r,vr, energy, self.lmax)
        
        dr = np.ediff1d(r[self.rhoIII<=r], to_end=r[-1]-r[-2])
        for l in range(0, self.lmax+1):
            # normalize u(r) such that
            #   /rmax            2
            #   |       dr |u(r)|  = 1
            #   /rhoIII
            u2int = np.sum(abs(u[self.rhoIII<=r,l])**2 * dr)
            u[:,l] /= np.sqrt(u2int)

            # create B-splines for R_l(r) = u_l(r) / r
            t,c,k = interpolate.splrep(r, u[:,l]/r)
            assert k == self.kspl
            
            self.fIIIknots[:,l] = t
            self.fIIIcoeffs[:,l] = c


            
class MuffinTinPotential(object):
    """
    muffin tin potential for a continuum multiple scattering calculation

    Parameters
    ----------
    atomlist       :  list of tuples `(Zat,[x,y,z])`
                      atomic number and cartesian coordinates (in bohr) of each atom in the molecule
    lmax           :  int, optional
                      highest angular momentum at which the spherical wave expansion is truncated
    charge         :  int, optional 
                      overall charge of the molecule, (total nuclear charge) - (number of electrons)
    chargeIII      :  int, optional
                      non-screened charge felt by an electron in region III
    potential_type :  str ('atomic' or 'molecular'), optional
                      The potential around each atom in region I is spherically symmetric.
                      The radial potential in each atomic sphere is either obtained as the
                      DFT potential of a free atom ('atomic') or by spherically
                      averaging the superpositions of atomic DFT potentials around each center ('molecular'). 
                      'molecular' contains the contribution of one atom to the potential of a
                      neighbouring atom. 'atomic' allows to reuse the same potential
                      and wavefunctions for all atoms of the same type.
    Nr             :  int, optional
                      number of grid points for radial integration using 
                      Gauss-Chebyshev quadrature rule
    debug          :  int, optional
                      if > 0, additional output is generated
    
    """
    def __init__(self, atomlist,
                 lmax=8,
                 charge=0, chargeIII=+1,
                 potential_type="molecular", Nr=2000, debug=0):
        print( " " )
        print( "     ***********************************" )
        print( "     *      Muffin Tin Potential       *" )
        print( "     ***********************************" )
        print( " " )

        print_geometry(atomlist)
        
        self.atomlist = atomlist        
        self.Nr = Nr
        self.debug = debug

        self.R = Regions(self.atomlist, charge=charge, chargeIII=chargeIII, lmax=lmax)

        #
        self.lmax = lmax
        
        # number of angular momentum channels (l,m)
        self.nlm = (self.lmax+1)**2
        # atomic positions
        self.nat = len(self.atomlist)

        # bring geometry into a format that can be passed to the Fortran extensions.
        self.atpos = np.zeros((3,self.nat))
        for i,(Zi,posi) in enumerate(self.atomlist):
            self.atpos[:,i] = posi

        # potentials in region I
        assert potential_type in ['atomic', 'molecular']
        self.potential_type = potential_type
        if self.potential_type == "molecular":
            # superposition of atomic densities averaged spherically around each atom
            if self.R.nelec > 1:
                xc = XCFunctionals.libXCFunctional('lda_x', 'lda_c_xalpha')

                rho = atomic_superposition_density(atomlist)                
                potential = effective_potential_func(atomlist, rho, xc)

                # sanity check, integrated electronic charge
                qnuc, qelec, qtot = total_charge(atomlist, rho)
                assert abs(-qelec - self.R.nelec) < 1.0e-2, "integral of electronic density and number of electrons differ, nelec = %d, integrated charge = %4.4f" % (self.R.nelec, qelec)

            else:
                # single electron system, only nuclear attraction needed
                potential = effective_potential_func(atomlist, None, None, nelec=1)
        else:
            assert charge == 0
            potential = None
        
        self.R.get_potential(potential)
        
    @T.timer
    def _solve_wavefunctions(self, energy):
        self.energy = energy
        self.R.solve_wavefunctionsI(energy)
        if energy < 0.0:
            self.R.solve_wavefunctionsIII(energy)
        
    def potential(self, x,y,z):
        """
        evaluate the muffin tin potential on a grid
                         
        - region I     :   atomic potentials, which are spherically symmetric
        - region II    :   constant potential
        - region III   :   Coulomb potential `-chargeIII/r`
                 
        Parameters
        ----------
        x,y,z      :  np.ndarray (any shape)
                      arrays with cartesian coordinates of grid points

        Returns
        -------
        pot        :  np.ndarray (same shape as `x`) 
                      array with muffin tin potential at the grid points, `V(x,y,z)`

        """
        # call Fortran extension
        pot = cms.muffin_tin_potential(self.R.vconst, self.R.chargeIII, self.R.rhoIII, self.R.radii,
                                       self.atpos, self.R.kspl, self.R.pot_knots, self.R.pot_coeffs,
                                       x,y,z)

        return pot

    ##################################################
    #                                                #
    # bound orbitals (energy E < 0)                  #
    #                                                #
    ##################################################
    
    @T.timer
    def _match_regions_bound(self):
        """
        compute matching matrix M and normalization for E < 0

        Returns
        -------
        mlogcond  : float
                    `-log(cond(M''))`, negative logarithm of condition number,
                    at eigen energies E, M''(E) should be singular and `logcond`=inf
                    
        """
        assert self.energy < 0.0
        M = ms.matching_matrix(self.energy, self.R.vconst, self.R.chargeIII, self.R.rhoIII, self.R.radii,
                               self.atpos, self.R.kspl,
                               self.R.fIknots, self.R.fIcoeffs,
                               self.R.fIIIknots, self.R.fIIIcoeffs)
        normalization = ms.normalization(self.energy, self.R.vconst, self.R.chargeIII, self.R.rhoIII, self.R.radii,
                                         self.atpos, self.R.kspl,
                                         self.R.fIknots, self.R.fIcoeffs,
                                         self.R.fIIIknots, self.R.fIIIcoeffs, self.Nr)
        if self.energy >= self.R.vconst:
            # For E >= V_II the matching matrix should be hermitian, i.e. M = M^+
            # because
            #  J_{L',L}(-R) = [J_{L,L'}(R)]^*
            #  N_{L',L}(-R) = [N_{L,L'}(R)]^*
            #
            err_herm = la.norm(M - M.conjugate().transpose())
            if (self.debug > 0): print( "|M - M^+|= %e" % err_herm )
            assert err_herm < 1.0e-6 * la.norm(M)
        else:
            # For E < V_II the matching matrix is neither symmetric nor hermitian
            # because
            #  I_{L',L}(-R) = (-1)^(l+l') [I_{L,L'}(R)]^*
            #  K_{L',L}(-R) = (-1)^(l+l') [K_{L,L'}(R)]^*
            #
            # We can convert the matching matrix into a Hermitian
            # one by the transformation
            #
            #     M' = D.M.D
            #
            # where D is a diagonal matrix with entries +1,-1,+imag,-imag
            # (with imag^2 = -1). 
            # The roots of the determinant are not affected, since
            #   det(M') = det(D) det(M) det(D) = det(D)^2 det(M)   with |det(D)|=1
            #
            # The diagonal matrix is
            #                     l
            # D             = imag  delta   delta         
            #  (i,L),(j,L')              ij      L,L'
            #
            l_arr, m_arr = enumerate_lm(self.lmax)
            # diagonal elements of D are phases
            phases = np.array( (self.nat+1) * list(1j**l_arr) )
            # M' = D.M.D
            # Because D is diagonal the matrix product becomes
            #
            # M'   = M    (D   D  )           
            #   ij    ij    ii  jj
            #
            # where i and j stand for the multiindices (i,L) and (j,L')
            M *= np.outer(phases, phases)

            # check that M' is indeed hermitian
            err_herm = la.norm(M - M.conjugate().transpose())
            if (self.debug > 0): print( "|M' - M'^+|= %e" % err_herm )
            assert err_herm < 1.0e-6 * la.norm(M)
            
        # At the eigenenergies E the matrix M(E) becomes singular, so M(E) has an
        # eigenvalue e = 0.
        #
        #  M(E) x  = e x
        #        e      e
        #
        # The resulting wavefunction should be normalized, the norm of
        # the wavefunction is given by a quadratic form of the eigenvector
        #
        #                   2         *
        #  |<Psi(E)|Psi(E)>|  = sum  x  S(E)   x
        #                          i  i     ij  j
        #
        # The overlap matrix is only approximate since it accounts only for the probability
        # density in regions I and III. S is diagonal
        #
        #                    2
        #  S (E) = delta    n                n  is the vector returned by `ms.normalization(...)`
        #   ij          ij   i                i
        #
        # To obtain normalized wavefunctions we can solve a generalized eigenvalue problem
        #
        #   M(E) x  = e S(E) x
        #         e           e
        #
        # which can be transformed into an eigenvalue problem by symmetric orthogonalization (Loewdin)
        #             +
        #   M''(E) = X M(E) X              M''(E) x''  = e x''
        #                                          e        e
        #
        # The transformation matrix is X = S^{-1/2}
        #
        #   X   = delta   1/n
        #    ij        ij    i
        #
        #   M'' (E) = M   / (n * n )
        #    ij        ij     i   j
        #
        # References
        # ----------
        # Szabo & Ostlund, section 3.4.5 'Orthogonalization of the Basis'
        #
        M /= np.outer(normalization, normalization)
        
        # At eigenenergies E the matrix M(E) becomes singular and
        # the condition number goes up.
        mlogcond = -np.log(la.cond(M))
        if self.debug > 0:
            print( "energy= %e    -log(cond(M''))= %s" % (self.energy, mlogcond) )
        
        # save symmetrically orthogonalized matching matrix M''
        self._matchM = M
        # save normalization coefficients n_i for constructing S^{-1/2}
        self._normalization = normalization
        
        return mlogcond
        
    def find_eigenstates(self, search_energies,
                         alpha=0.2, kappa_zero=1.0e-4):
        """
        search an energy range for eigen energies where M(E) is singular

        Parameters
        ----------
        search_energies :  np.1darray, (shape (n,))
                           array with energies (in Hartree) that are scanned for minima.
                           The more points are given per energy interval, the less likely it
                           is to miss minima in between. The points can be spaced linearly, e.g.
                           `search_energies=np.linspace(emin, emax, n)`
                           or non-linearly, if you already know the approximate positions of eigenstates.
        alpha           :  float, optional (0 < alpha < 1.0)
                           The approximate positions of the eigen energies
                           are refined by Golden-section search. The bounds for the search are set
                           to E +- alpha*dE, where E is the initial guess of the minimum and dE is
                           the distance to the next lower of higher maximum.
        kappa_zero      :  float, optional
                           At `kappa = sqrt(2*|E-V_II|) = 0`, the matching matrix blows up, therefore
                           energy points, where `kappa < kappa_zero` are excluded from the search.
        
        Returns
        -------
        orbitals     :  list of instances of `MSWavefunction`
                        the energy of an orbital `orb = orbitals[i]` is stored as `orb.energy`

        """
        print( " " )
        print( "     **********************************" )
        print( "     *   searching for bound states   *" )
        print( "     **********************************" )
        print( " " )
        # define objective function of energy that should be minimized, at the minimum
        # there should an eigenvalue
        def func(energy):
            # At kappa = 0, i.e. where E = V_II, the matching matrix
            # always becomes singular because the n_l(kappa*rho) and
            # k_l(kappa*rho) diverge. Therefore eigenstate cannot be
            # found if they lie at kappa = 0
            kappa = np.sqrt(2.0*abs(energy - self.R.vconst))
            if kappa < kappa_zero:
                # abort minimization of -log(cond(M))
                print( "too close to kappa = 0!" )
                raise StopIteration()
            
            self._solve_wavefunctions(energy)
            # f(E) = -log(cond(M(E)))
            f = self._match_regions_bound()
            return f

        print( "The energy range [%12.8f, %12.8f] x Hartree is scanned for bound states." % (search_energies.min(), search_energies.max()) )
        n = len(search_energies)
        # remove points too close to kappa = 0
        kappa = np.sqrt(2.0*abs(search_energies - self.R.vconst))
        search_energies = search_energies[kappa > kappa_zero]
        # measure of singularity of matrix M, -log(cond(M))
        mlogcond_range = np.zeros(n)

        for i,energy in enumerate(search_energies):
            mlogcond_range[i] = func(energy)
            # show progress
            if (((i+1)*100) % n == 0):
                if (self.debug > 0): print( " %6.1d  /%6.1d complete" % (i+1,n) )

        print( "" )
        print( "Minima of f(E)=-log(cond(M(E))) are located and optimized further." )
        print( "Some of the minima are spurious and can be identified because the" )
        print( "norm of the C^I coefficients is close to zero." )
        print( "" )
        # find indices of local minima
        minima = signal.argrelmin(mlogcond_range)[0].tolist()
        # and indices of local maxima
        maxima = signal.argrelmax(mlogcond_range)[0].tolist()

        # save energy scan
        data = np.vstack((search_energies, mlogcond_range)).transpose()
        fh = open("/tmp/mlogcond.dat", "w")
        print ( "# Energy E      -log(cong(M(E)))" , file=fh )
        np.savetxt(fh, data)
        fh.close()

        if self.debug > 0:
            # plot energy scan
            import matplotlib.pyplot as plt
            plt.plot(search_energies, mlogcond_range, label=r"-log(cond($M(E)$))")
            # show minima and minima as dots
            for i in minima:
                plt.plot([search_energies[i]], [mlogcond_range[i]], "o", color="red")
            for i in maxima:
                plt.plot([search_energies[i]], [mlogcond_range[i]], "o", color="blue")
            plt.xlabel("energy / Hartree")
            plt.ylabel("-log(cond($M(E)$))")
            plt.legend()
            plt.show()
            plt.close()

        if len(maxima) == 0:
            # No maximum was found, bracket minimum by end points
            maxima += [0,-1]
        # Each local minimum should be bracketed by two local maxima
        if search_energies[minima[0]] < search_energies[maxima[0]]:
            # first extremum is a minimum, so
            #   maxima[i-1] < minima[i] < maxima[i]
            maxima = [0] + maxima
            # After prepending the first point, we have
            #   maxima[i  ] < minima[i] < maxima[i+1]
        if search_energies[minima[-1]] > search_energies[maxima[-1]]:
            # last extremum is a minimum, which is not bracketed
            # by two maxima
            maxima = maxima + [-1]
            # After appending last point, we have
            #  maxima[i  ] < minima[i] < maxima[i+1]
            # for all minima
        assert len(minima) == len(maxima)-1
        # list of orbital energies
        energies = []
        funcmins = []
        # list of orbitals (instances of `MSWavefunction`)
        orbitals = []

        # Refine energy at each local minimum
        for i in range(0, len(minima)):
            # initial guess
            emin0 = search_energies[minima[i]]
            if (self.debug > 0):
                print( "" )
                print( "minimum at E= %e" % emin0 )
            # maxima that bracket this local minimum
            emax_lower = search_energies[maxima[i]  ]
            emax_upper = search_energies[maxima[i+1]]
            assert emax_lower < emin0 < emax_upper
            # We search for a minimum around emin0 by Golden search
            # (https://en.wikipedia.org/wiki/Golden-section_search)
            # which assumes that there is a single minimum in the interval [l,u]
            l = (1.0-alpha)*emin0 + alpha*emax_lower
            u = (1.0-alpha)*emin0 + alpha*emax_upper
            # find minimum of func(E) = -log(cond(M(E))) in the interval [l,u]
            try:
                emin = minimize_golden(func, l, u)
            except StopIteration:
                continue
            fmin = func(emin)
            assert self.energy == emin

            orbitals += self._get_bound_wavefunctions()

            energies.append(emin)
            funcmins.append(fmin)

        print( " ------------------------------------------------" )
        print( " Orbital                   Energy                " )
        print( "                 Hartree             eV          " )
        print( " ------------------------------------------------" )
        for i,orb in enumerate(orbitals):
            print( "  %4.1d        %12.8f        %12.8f" % (i+1, orb.energy, orb.energy * AtomicData.hartree_to_eV) )
        print( "" )

        if self.debug > 0:
            # plot refined energy scan
            import matplotlib.pyplot as plt
            plt.plot(search_energies, mlogcond_range, label=r"-log(cond($M(E)$))")
            for e,f in zip(energies, funcmins):
                plt.plot([e],[f], "o", color="red")
            plt.xlabel("energy / Hartree")
            plt.ylabel("-log(cond($M(E)$))")
            plt.legend()
            plt.show()
            plt.close()
        
        return orbitals
    
    def _get_bound_wavefunctions(self, singular_thresh=1.0e-8):
        """
        find wavefunctions at the current energy by diagonalizing the symmetrically orthogonalized
        matching matrix.

        Parameters
        ----------
        singular_thresh     :   float, optional
                                threshold for considering an eigenvalue of M(E) to be 0.

        Returns
        -------
        waves               :   list of instances of `MSWavefunction`
                                wavefunctions with energy E,
                                `waves` contains one wavefunction for each eigenvalue `|w| < sing_thresh`.
                                If there are no such eigenvalues, and empty list is returned.
        """
        assert self.energy < 0, "Energy of bound wavefunction has to be negative!"
        ndim,ndim = self._matchM.shape
        # diagonalize symmetrized matrix M' = X^+ M(E) X
        #   M'(E) C' = e C'
        #
        evals,evecs = la.eigh(self._matchM)

        # transform back from orthogonalized basis
        #                          +         + +           +
        # C = X.C'      such that C .S.C = C' X S X C' = C' C'
        #
        #           -1/2
        #   X   = (S    )   = delta   1/n
        #    ij          ij        ij    i
        for j in range(0, ndim):
            # C_ij = (X C')_ij = 1/n_i C'_ij
            evecs[:,j] /= self._normalization
            
        if (self.energy < self.R.vconst):
            # To turn the matching matrix into a hermitian matrix, rows and
            # columns were multiplied with phases. To reverse this, we have
            # to multiply the columns of eigenvectors by the same phases.
            l_arr, m_arr = enumerate_lm(self.lmax)
            # diagonal elements of D are phases
            phases = np.array( (self.nat+1) * list(1j**l_arr) )

            # evecs[i,j] = evecs[i,j] * phases[i]
            for j in range(0, ndim):
                evecs[:,j] *= phases
                
        # The multiple scattering wavefunction is joined smoothly at the boundaries I/II and II/III
        # if M(E) has at least one 0 eigenvalue. If there are more than one such eigenvalue, the
        # wavefunctions are degenerate.

        orbitals = []
        # indices of singular eigenvalues
        indices_singular = np.where(abs(evals) <= singular_thresh)[0]
        #print( "all eigenvalues = %s" % evals )
        #print( "singular values = %s" % evals[indices_singular] )
        for i in indices_singular:
            coeffsA = evecs[:,i]
            # Not all singularities of the matching matrix correspond
            # to smoothly matched solutions of the Schroedinger equation.
            # The spurious solutions have in common that only the coefficients
            # in region III (C^III) are non-zero, while the coefficients in region I
            # (C^I) are exactly zero. It is clear that such a situation is cannot be
            # reconciled with a a continuous wavefunction.

            # The eigenvectors contain the coefficients A^II0 and A^IIi, we need to
            # compute the coefficients C^III and C^Ii.
            coeffsC = ms.convert_coefficients_atoc(self.energy, self.R.vconst, self.R.chargeIII, self.R.rhoIII, self.R.radii,
                                                  self.atpos, self.R.kspl, self.R.fIknots, self.R.fIcoeffs, self.R.fIIIknots, self.R.fIIIcoeffs,
                                                  coeffsA)
            # Spurious solutions are recognized by |C^I|^2 = 0
            C_I = coeffsC[(self.lmax+1)**2:,:]
            nrm2 = np.sum(abs(C_I)**2)
            if (self.debug > 0): print( "norm^2 of C coefficients in region I : %e" % nrm2 )
            if nrm2 < 1.0e-10:
                if (self.debug > 0): print( "skipping spurious solution at energy= %s" % self.energy )
                continue
            
            orbital = MSWavefunction(self.energy, self.R.vconst, self.R.chargeIII, self.R.rhoIII, self.R.radii,
                                     self.atpos, self.R.kspl, self.R.fIknots, self.R.fIcoeffs, self.R.fIIIknots, self.R.fIIIcoeffs,
                                     coeffsA, evals[i])
            orbitals.append(orbital)

        print( "Found %d orbitals at energy E=%s" % (len(orbitals), self.energy) )
        
        return orbitals

    ##################################################
    #                                                #
    # continuum orbitals (energy E > 0)              #
    #                                                #
    ##################################################

    def _match_regions_continuum(self):
        """
        construct the matrix for the inhomogeneous system of linear equations
        (16) and (17) in Ref. [1] and solve 

            M.x = rhs

        for all right-hand sides.
        """
        # call Fortran extension
        # The (S) wavefunctions fulfill the boundary conditions
        #                    -1/2  -1                                   
        #   Psi      ~  (pi k)     r    sum   [ sin(th  ) delta     +  K     cos(th  ) ] Y
        #      III,L                0      L'         l'       L,L'     L,L'       l'     L'
        M, rhs = cms.matching_matrix(self.energy, self.R.vconst, self.R.chargeIII, self.R.rhoIII, self.R.radii,
                                     self.atpos, self.R.kspl, self.R.fIknots, self.R.fIcoeffs)
        
        # check that matching matrix M is hermitian
        err = la.norm(M - M.conjugate().transpose())/la.norm(M)
        assert err < 1.0e-10, "matching matrix is not hermitian, |M-M^+|/|M| = %e" % err
        
        # solve matching equation M.x = b
        self._sol = la.solve(M, rhs)
            
    def get_continuum_wavefunction(self, energy,
                                   boundcond='standing', wavetype='spherical',
                                   l=0,m=0, kth=0.0, kph=0.0):
        """
        find the wavefunction which has the asymptotic form in eqn. (20) in Ref. [2]

        Solutions are either labeled by the quantum numbers (l,m) of the partial
        wave expansion in region III (if wavetype is set to 'spherical') 
        or by the angles kth and kph of the propagation direction of the plane wave
        (if wavetype is set to 'plane').

        Parameters
        ----------
        energy     :  float > 0
                      energy of continuum orbital
        boundcond  :  str, optional
                      impose boundary conditions on wavefunctions.

                        * 'standing': see eqn. (20), wavefunction with K-matrix normalization,
                        * 'incoming': incoming wave (-) normalization according to eqn. (26),
                          wavefunction is complex,
                        * 'outgoing': outgoing wave (+) normalization according to eqn. (25),
                          wavefunction is complex.

        wavetype   :  str, optional
                      specify type of wavefunction.
                      
                        * 'spherical' requests that the wavefuntions are labeled by the quantum numbers `(l,m)`.
                        * 'plane' requests that partial waves are combined to form a plane wave with the 
                          propagation direction determined by the angles `kth` and `kph`.

        l,m        :  int, optional
                      Asymptotic angular momentum in region III of continuum orbital for
                      'spherical' waves
        kth,kph    :  float, optional
                      propagation direction of 'plane' wave, the spherical coordinates
                      of the wave vector are (k,kth,kph)

        Returns
        -------
        wfn        :  callable
                      `wfn(x,y,z)` evaluates the `(l,m)` wavefunction on a grid
        """
        self._solve_wavefunctions(energy)
        self._match_regions_continuum()
        
        assert self.energy >= 0, "Energy of continuum wavefunction has to be positive!"
        # We only need the K-matrix, so set x,y,z to dummy values
        x,y,z = np.ones(3), np.ones(3), np.ones(3)
        # evaluate standing waves
        kmat, wfn = cms.wavefunctions(self.energy, self.R.vconst, self.R.chargeIII, self.R.rhoIII, self.R.radii,
                                      self.atpos, self.R.kspl, self.R.fIknots, self.R.fIcoeffs, self._sol, 0,
                                      x,y,z)
        K = kmat
        Id = np.eye(K.shape[0])
            
        # check that K-matrix is hermitian
        # NOTE: Dill & Dehmer have a symmetric K-matrix, but we use complex spherical harmonics.
        #       I am not sure what the consequences of a hermitian K-matrix are?
        print("|Im(K^(S))|= %e" % la.norm(kmat.imag))
        print("|K - K^T|= %e" % la.norm(kmat - kmat.transpose()))
        err_harm = la.norm(kmat - kmat.conjugate().transpose())
        print("|K - K^dagger|= %e" % err_harm)
        # check that S-matrix is unitary
        # S (1 - i*K) = (1 + i*K)
        S = np.dot(Id+1j*K , la.inv(Id-1j*K))
        err_unitary = la.norm( Id - np.dot(S, S.conjugate().transpose()) )
        print("|S.S^dagger-Id|= %e" % err_unitary)
        
        debug = 0
        if debug == 1:
            # plot K-matrix
            import matplotlib.pyplot as plt
            fig,axes = plt.subplots(1,2)
            axes[0].set_title("K-matrix (real part)")
            axes[0].imshow(kmat.real)
            axes[1].set_title("K-matrix (imag. part)")
            axes[1].imshow(kmat.imag)
            plt.show()
            
        # construct complex wavefunction from K-matrix normalized wavefunctions
        # according to eqn. (28) and (29)
        if boundcond == "incoming":
            print("incoming wave (-) normalization")
            # exp(-i*k*r) = cos(kr) - i * sin(kr)
            # eqn. (29-)
            C = la.inv(Id + 1j*K)
        elif boundcond == "outgoing":
            print("outgoing wave (+) normalization")
            # exp(+i*k*r) = cos(kr) + i * sin(kr)
            # eqn. (29+)
            C = -la.inv(Id - 1j*K)
        elif boundcond == "standing":
            # Wavefunction
            print("standing wave")
            C = np.eye(K.shape[0])
        else:
            raise ValueError("Illegal value '%s' for option 'boundcond'" % boundcond)

        # number of angular momentum channels (l,m)
        nlm = (self.lmax+1)**2
        # The angular momentum quantum numbers (l,m) are enumerated in the order
        # (0,0),
        # (1,0), (1,+1), (1,-1),
        # (2,0), (2,+1), (2,-1), (2,+2), (2,-2),
        # ...., (lmax,-lmax)
        angmoms = []
        for l in range(0, self.lmax+1):
            for m in range(0, l+1):
                angmoms.append( (l,m) )
                if m > 0:
                    angmoms.append( (l,-m) )

        assert len(angmoms) == nlm

        if wavetype == "spherical":
            # Which partial wave solution (l,m) is required?
            # lm is an index into the array of angular momentum channels (l,m)
            #  [(0,0), (1,0),(1,-1),(1,+1), ...]
            for lm1,(l1,m1) in enumerate(angmoms):
                if l1 == l and m1 == m:
                    # (l,m) corresponds to multiindex lm
                    lm = lm1
                    break

            A = np.zeros(nlm, dtype=complex)
            A[lm] = 1.0

        elif wavetype == "plane":
            if boundcond == "outgoing":
                # Partial waves are combined to have the asymptotic form of eqn. (31)
                # a^(+)_L coefficients of eqn. (37) for rotatione matrix R=Id
                A = np.zeros(nlm, dtype=complex)
                for lm1,(l1,m1) in enumerate(angmoms):
                    if m1 == 0:
                        A[lm1] = -1j**l1 * np.sqrt((2*l1+1)/(4.0*np.pi))
            elif boundcond == "incoming":
                # Partial waves are combined to have the asymptotic form of eqn. (32)
                # a^(-)_L coefficients of eqn. (41) for rotatione matrix R=Id
                k = np.sqrt(2*self.energy)
                eta = -self.R.chargeIII/k
                
                A = np.zeros(nlm, dtype=complex)
                # spherical harmonics Y_{l,m}(kth,kph)
                ysph = sphharm(kth,kph, self.lmax)
                for lm1,(l1,m1) in enumerate(angmoms):
                    # Coulomb phase shift, sigma_l = arg Gamma(l+1+i*eta)
                    g = mpmath.gammainc(l1+1.0+1.0j*eta)
                    sigma_l1 = np.angle( complex(g) )
                    # 
                    A[lm1] = 1j**l1 * np.exp(-1j*sigma_l1) * ysph[lm1].conjugate()
            else:
                raise ValueError("Options boundcond='%s' and wavetype='%s' are not compatible." % (boundcond, wavetype))

        else:
            raise ValueError("Illegal value '%s' for option 'wavetype'" % wavetype)        

        coeffs = np.dot(A, C)

        # linear combination of partial waves
        orbital = CMSWavefunction(self.energy, self.R.vconst, self.R.chargeIII, self.R.rhoIII, self.R.radii,
                                  self.atpos, self.R.kspl, self.R.fIknots, self.R.fIcoeffs, self._sol,
                                  coeffs)
        return orbital
    
    ##################################################
    #                                                #
    # photoelectron angular distributions            #
    #                                                #
    ##################################################
    
    def evaluate_orbitals(self, orbitals):
        """

        A set of orbitals are evaluated on a multicenter Becke grid. The resolution
        of the radial and angular grid is set in `MolecularIntegrals.settings`.

        Parameters
        ----------
        orbitals   :  list of `norb` callables
                      `orbitals[b](x,y,z)` should evaluate the `b`-th orbital on a grid

        Returns
        -------
        grid       :  tuple of numpy arrays `(x,y,z, w)`
                      cartesian coordinates and weights for numerical quadrature
        orbs       :  numpy array of shape `(norb,len(x))` 
                      values of orbitals on the grid

        """
        # spherical Becke grids around each atom
        points, weights, volumes = multicenter_grids(self.atomlist,
                                                     radial_grid_factor=settings.radial_grid_factor,
                                                     lebedev_order=settings.lebedev_order)
        # combine multicenter grids into a single grid
        # The quadrature rule for integration on this grid reads:
        #   /
        #   | f(x,y,z) dV = sum  w  f(x ,y ,z ) 
        #   /                  i  i    i  i  i
        x,y,z, w = join_grids(points, weights, volumes)
        # number of grid points
        npts = len(x)
        
        # evaluate all bound orbitals on the grid
        # number of orbitals
        norb = len(orbitals)
        orbs = np.zeros((norb, npts), dtype=complex)
        if (self.debug > 0): print( "evaluating orbitals..." )
        for b,orbital in enumerate(orbitals):
            orbs[b,:] = orbital(x,y,z)

        grid = (x,y,z, w)
        
        return grid, orbs

    def transition_dipoles(self, energy, grid, orbs):
        """

        compute transition dipoles matrix elements between bound orbitals (index by `b`)
        and the CMS continuum wavefunctions (labeled by the asymptotic angular momentum `(l,m)`)

        .. code-block:: none

            ( TDx )                 ( x )
            ( TDy )   = <Psi (k)  | ( y ) | Psi >
            ( TDz )         l,m     ( z )      b

        which are needed to evaluate eqn. (59) in Ref. [2].

        The integration is done numerically using the multicenter Becke grid generated in
        `evaluate_orbitals(...)`.

        Parameters
        ----------
        energy     :  float > 0
                      energy of continuum orbitals
        grid       :  tuple of numpy arrays `(x,y,z, w)`
                      cartesian coordinates and weights for numerical quadrature
        orbs       :  numpy array of shape `(norb,len(x))` 
                      values of bound orbitals on the grid

        Returns
        -------
        tdip       :  np.3darray of shape `(nlm,norb,len(x))`
                      transition dipoles between bound and free orbitals,
                      `tdip[lm,b,:]` is the cartesian transition dipole vector between the bound orbital with index `b`
                      and the K-matrix normalized continuum orbital with asymptotic momentum `lm`.
        norms2     :  np.1darray of shape `(norb,)` 
                      norms**2 of bound orbitals, for checking if the resolution
                      of the multicenter grid is sufficient

        """
        self._solve_wavefunctions(energy)
        self._match_regions_continuum()
        if (self.debug > 0): print( "computing transition amplitudes by numerical integration..." )
        x,y,z, w = grid
        # call Fortran extension for computing wavefunctions and
        # summing over grid points
        kmat, tdip, norms2, projs2 = cms.transition_dipoles(self.energy, self.R.vconst, self.R.chargeIII, self.R.rhoIII, self.R.radii,
                                                            self.atpos, self.R.kspl, self.R.fIknots, self.R.fIcoeffs, self._sol, 0,
                                                            x,y,z, w,
                                                            orbs)

        if self.debug > 0:
            print( "" )
            print( "The bound orbitals should be normalized and approximately orthogonal" )
            print( "to all continuum orbitals, otherwise the transition dipoles depend" )
            print( "on the choice the origin." )
            print( "" )
            print( " Orbital            Norm              Proj. onto Continuum"     )
            print( "    b         |<Psi_b|Psi_b>|^2     sum_lm |<Psi_lm|Psi_b>|^2"  )
            print( " ------------------------------------------------------------"  )
            for b,(nrm,prj) in enumerate(zip(norms2,projs2)):
                print( "  %3.1d           %e           %e" % (b+1,nrm,prj)  )
            print( "" )

        # check that K-matrix is hermitian
        err = la.norm(kmat -kmat.conjugate().transpose())
        #assert err < 1.0e-10, "|K-K^+|= %e" % err
        ### DEBUG
        if self.debug > 0:
            print( "|Re(K)|= %e" % la.norm(kmat.real) )
            print( "|Im(K)|= %e" % la.norm(kmat.imag) )
        ###

        # save K-matrix and cartesian transition dipoles
        self._kmat = kmat
        self._tdip = tdip
        
        return tdip, norms2
        
    def photoelectron_distribution(self, energy, grid, orbs,
                                   pol=0):
        """

        photoelectron angular distribution for isotropically oriented
        emsemble of molecules characterized by the parameters

          sigma (total photoionization cross section),
          beta1 (only non-zero for chiral molecules) and
          beta2 (anisotropy parameter)

        Parameters
        ----------
        energy     :  float > 0
                      energy of continuum orbitals
        grid       :  tuple of numpy arrays `(x,y,z, w)`
                      cartesian coordinates and weights for numerical quadrature
        orbs       :  numpy array of shape `(norb,len(x))` 
                      values of bound orbitals on the grid
        pol        :  int, optional
                      polarization of light, 0 (linear), -1 (left), +1 (right)

        Returns
        -------
        pad      : np.2darray of shape `(norb,3)`
                   PAD parameters  `sigma,beta1,beta2 = pad[b,:]` for each bound orbital `b`
   
        """
        # number of bound orbitals
        norb = orbs.shape[0]
        # transition dipole moments between bound orbitals and continuum orbitals
        # at given energy (creates self._kmat and self._tdip)
        self.transition_dipoles(energy, grid, orbs)
        # convert cartesian transition amplitudes to incoming wave (-) normalization
        tdip_in = cms.transform_tdip_in(self.energy, self.R.chargeIII, self._kmat, self._tdip, self.lmax)
        # The array returned by `transform_tdip_in` has the shape (nlm,norb,3),
        # but `pad_iso` expects an array of shape (norb,3,nlm).
        tdip_in = np.moveaxis(tdip_in, 0, -1)
        if (self.debug > 0): print( "photoelectron angular distribution ..." )
        pad = photo.pad_iso(self.energy, pol, tdip_in, self.lmax)

        if self.debug > 0:
            print( ""                                                               )
            print( "     Photoionization cross section (isotropic averaging)      " )
            print( "                                                              " )
            print( "  dsigma   sigma                                              " )
            print( "  ------ = ----- [1 + beta  P (cos(th)) + beta  P (cos(th)) ] " )
            print( "  dOmega   4 pi           1  1                2  2            " )
            print( "                                                              " )
            pol2str = {0 : "0 (linear)", -1 : "-1 (left)", +1: "+1 (right)"}
            print( "  polarization = %s" % pol2str[pol] )
            print( "  maximum angular momentum l = %d" % self.lmax )
            print( "  photokinetic energy = %e Hartree ( %e eV )" % (self.energy, self.energy * AtomicData.hartree_to_eV) )
            print( " "                                                               )
            print( "  Orbital      sigma             beta1             beta2       " )
            print( "  -------------------------------------------------------------" )
            for b in range(0, norb):
                sigma, beta1, beta2 = pad[b,:]
                print( "    %3.1d       %e      %+e      %+e" % (b+1, sigma, beta1, beta2) )
            print( "" )

        return pad

    def eigenchannel_analysis(self):
        """

        For each transition

           bound orbital b -----> continuum

        the eigenchannel (see Ref. [3]) continuum orbital with the largest
        transition dipole is selected. At shape resonances this orbital should
        resemble a virtual orbital. Because it lies above the ionization threshold,
        it is only metastable.

        Returns
        -------
        orbitals   :  list instance of `CMSWavefunction`
                      eigenchannel continuum orbitals,
                      `orbitals[b]` is the continuum orbital belonging to the transition `b --> continuum`

        """
        assert hasattr(self, '_kmat') and hasattr(self, '_tdip'), "Transition dipoles and K-matrix has to be calculated first!"
        # The hermitian K-matrix is diagonalized to get the eigenvectors U
        # and eigenphases mu.
        #
        # eqn. (2) in Ref. [4]
        #                               *
        # tan(pi mu ) delta     = sum  U     K     U
        #          i       i,j     L,L'  L,i   L,L'  L,j
        #
        eigvals, U = la.eigh(self._kmat)
        # eigenphases
        eigphases = np.arctan(eigvals)/np.pi % (2.0*np.pi)
        
        # eigenchannel wavefunction
        #
        # Phi  = sum U    Psi
        #    i    L   L,i    L
        #
        # The eigen function with the largest transition dipole moment
        # is returned

        # Transition dipole moments between K-matrix normalized continuum orbitals Psi_L
        # and bound states Psi_b are transformed to the basis of eigenchannel functions
        # Phi_i.
        #                       *
        # <Phi |r|Psi > = sum  U    <Psi |r|Psi >
        #     i      b     L    L,i     L      b
        tdip_eig = np.tensordot(U.conjugate(), self._tdip, axes=[0,0])

        tdip2_eig = np.sum(abs(tdip_eig)**2, axis=2)

        # For each bound orbital, we select the continuum eigenstates which has the
        # largest transition dipole moment
        
        nc,norb = tdip2_eig.shape
        # Coefficients of most important (in terms of oscillator strength) continuum
        # eigenstate in the basis of K-matrix normalized states are used to compute
        # eigenchannel orbitals for each transition (for each bound orbital)
        orbitals = []
        for b in range(0, norb):
            # index of continuum eigenstate with largest transition dipole to bound orbital b
            c = np.argmax(tdip2_eig[:,b])

            if self.debug > 0:
                print( "eigenchannel for  %d -> continuum  :   phase shift= %e    |tdip|^2= %e" % (b+1, eigphases[c], tdip2_eig[c,b]) )
            coeffs = U[:,c]
            
            orbital = CMSWavefunction(self.energy, self.R.vconst, self.R.chargeIII, self.R.rhoIII, self.R.radii,
                                      self.atpos, self.R.kspl, self.R.fIknots, self.R.fIcoeffs, self._sol,
                                      coeffs)
            # store eigen phase
            orbital.tags["eigen_phase"] = eigphases[c]
            # These are the same units as in `photo.f90`
            units = 8.0/3.0 * np.pi / AtomicData.speed_of_light * np.sqrt(2*self.energy)
            # store ionization cross section
            orbital.tags["sigma"] = 4.0 * np.pi/3.0 * units * tdip2_eig[c,b]
            
            orbitals.append(orbital)
            
        return orbitals
    
    
def minimize_golden(f, a, b, xtol=1e-16):
    """
    finds the local minimum of a univariate function f(x) 
    in the interval [a,b] by bisection assuming that there is 
    only one minimum.

    The Golden-section search is used to narrow the search interval
    down until `|b-a| < xtol`.

    References
    ----------
    https://en.wikipedia.org/wiki/Golden-section_search

    Parameters
    ----------
    f          : callable
                 strictly unimodal function on [a,b]
    a,b        : float (with a < b) 
                 initial boundaries of search interval
    xtol       : float
                 stop if size of search interval has shrunk below `xtol`

    Returns
    -------
    x0         : float
                 minimum

    Example
    -------
    >>> f = lambda x: (x-2)**2
    >>> x = minimize_golden(f, 1, 5)
    >>> x
    2.000009644875678

    """
    # golden ratio
    gr = (1+np.sqrt(5))/2
    
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    while abs(c - d) > xtol:
        if f(c) < f(d):
            b = d
        else:
            a = c
            
        # we recompute both c and d here to avoid loss
        # of precision which may lead to incorrect results
        # or infinite loop
        c = b - (b - a) / gr
        d = a + (b - a) / gr
        
    return (b + a) / 2


def print_geometry(atomlist):
    """
    print cartesian coordinates in XYZ format
    """
    print( "  Cartesian geometry (in Angstrom)" )
    print( " %d " % len(atomlist))
    print( "    ")
    for Z,pos in atomlist:
        x,y,z = map(lambda xyz: xyz*AtomicData.bohr_to_angs, pos)
        print( " %s   %+10.8e  %+10.8e  %+10.8e" % (AtomicData.atom_names[Z-1], x,y,z) )
    print( "    ")

def geometries_are_equal(atomlistA, atomlistB, eps=1.0e-3):
    """
    compare two molecular geometries
    """
    equal = True
    for (ZatA,posA),(ZatB,posB) in zip(atomlistA, atomlistB):
        equal = equal and (ZatA == ZatB)
        equal = equal and la.norm(np.array(posA)-np.array(posB)) < eps
    return equal

def project_orbital(basis_funcs, orbital):
    """

    project `orbital` onto a basis of orthonormal functions `basis_funcs`

    .. code-block:: none

      |proj> = sum  <b |orbital> |b >
                  i   i            i

    Parameters
    ----------
    basis_funcs   :  list of callables
        `basis_funcs[i](x,y,z)` should evaluate the i-th basis function on a grid, <x,y,z|b_i>.
    orbital       :  callable, instance of `MoldenFile.MolecularOrbital`
        `orbital(x,y,z)` evaluates <x,y,z|orbital>

    Returns
    -------
    projected     :  callable
        `projected(x,y,z)` evaluates <x,y,z|proj>

    """
    nbf = len(basis_funcs)
    # coefficients of orbital in the new basis
    coeffs = np.zeros(nbf, dtype=complex)
    for i,bfunc in enumerate(basis_funcs):
        coeffs[i] = overlap(orbital.atomlist, bfunc, orbital)

    # define wavefunction of projected orbital
    def projected(x,y,z):
        wfn = 0.0j*x
        for i,bfunc in enumerate(basis_funcs):
            wfn += coeffs[i] * bfunc(x,y,z)
        return wfn

    return projected


##################################################
#                                                #
# Testing                                        #
#                                                #
##################################################

def test_hmi_bound_states():
    """
    compute the bound eigenstates of the hydrogen molecular ion (H2+) using the 
    multiple scattering method in the muffin tin approximation.

    The eigen energies should agree with the first column (labelled "MT (L=0)") in
    table 2 of Ref. [5]_. The error in the energy eigenvalues should be on the order 
    of 0.0001 Hartree.

    The energies of the first few bound states should be:

    .. code-block: none
       ------------------------------------------------
        Orbital                   Energy                
                       Hartree             eV          
       ------------------------------------------------
           1         -1.03597887        -28.19043154
           2         -0.64397430        -17.52343984
           3         -0.44433307        -12.09092322
           4         -0.44433307        -12.09092322
           5         -0.35206202         -9.58009908

    .. rubric:: References
    .. [5] B. Wästberg, "Fully numerical solutions of the molecular Schrödinger equation: the MS-X alpha method for non-muffin-tin potentials." (1994) J. Phys. B: At. Mol. Opt. Phys. 27 2105

    """
    # choose resolution of multicenter grids
    settings.radial_grid_factor = 10      # controls size of radial grid 
    settings.lebedev_order = 23           # controls size of angular grid

    # The radius of H2+ is set to R=2 bohr
    atomlist = [
        (1, (0.0, 0.0, -1.0)),
        (1, (0.0, 0.0, +1.0))
        ]

    # The overall charge of the HMI is +1, however the single electron will
    # feel the nuclear attraction from both protons, therefore chargeIII=+2.
    muffin = MuffinTinPotential(atomlist,
                                lmax=8,
                                charge=+1, chargeIII=+2,
                                potential_type="molecular", Nr=2000, debug=0)

    # We search for bound state in the energy range [-2.0, -0.01] Hartree.
    search_energies = np.linspace(-2.0, -0.01, 2000)
    orbitals = muffin.find_eigenstates(search_energies)

    
    
if __name__ == "__main__":
    test_hmi_bound_states()
