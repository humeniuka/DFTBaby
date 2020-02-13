#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
photoangular distributions for isotropically oriented ensembles 
of molecules in the gas phase

Initial bound and final continuum orbitals are obtained by the 
multiple scattering (MS) and the continuum multiple scattering 
(CMS) methods, respectively.
"""
from __future__ import print_function

from DFTB import AtomicData
from DFTB.MolecularIntegrals import settings
from DFTB.MultipleScattering.MuffinTin import MuffinTinPotential, minimize_golden, save_cubefile

import numpy as np
from scipy import signal

class PhotoelectronAngularDistribution:
    """

    Parameters
    ----------
    muffin          :  instance of `MuffinTinPotential`
                       muffin tin potential, which provides the continuum orbitals

    """
    def __init__(self, muffin):
        self.muffin = muffin
        
    def compute_pads(self, bound_orbitals, pke,
                     pol=0, pad_file="/tmp/pad.dat", units="eV-Mb"):
        """

        compute photoelectron angular distributions (PADs) for isotropic ensemble 
        as a function of kinetic energy of the photoelectron using the CMS method
        
        Parameters
        ----------
        bound_orbitals  :  list of callables
                           `bound_orbitals[b](x,y,z)` should evaluate the `b`-th orbital
                           on a cartesian grid
        pke             :  np.1darray (shape (npke,))
                           photokinetic energies (in Hartree) for which the PADs should be calculated
        pol             :  int, optional
                           polarization of light, 0 (linear), -1 (left), +1 (right)
        pad_file        :  str, optional
                           path to file for storing table with PADs
        units           :  str, optional
                           units of energies and cross sections in PAD table
                             *  'eV-Mb' :   PKE in eV, sigma in megabarn
                             *  'au'    :   PKE in Hartree, sigma in bohr^2

        Returns
        -------
        pad             :  np.ndarray (shape (npke,norb,3))
                           `pad[k,b,:]` contains the three parameters
                           `sigma`,`beta1` and `beta2` for ionization from orbital `b` into 
                           the continuum at energy `pke[k]`.
        
        """
        print( " " )
        print( "     *******************" )
        print( "     *      PADs       *" )
        print( "     *******************" )
        print( " " )
        # number of bound orbitals
        norb = len(bound_orbitals)
        # number of energies
        npke = len(pke)

        # find values of bound orbitals on a Becke grid for numerical integration
        grid, orbs = self.muffin.evaluate_orbitals(bound_orbitals)
        
        # compute orientation-averaged PAD for each energy
        pad = np.zeros((npke,norb,3))
        for i,energy in enumerate(pke):
            if (self.muffin.debug > -1):
                print( "%3.1d of %3.1d   PKE = %e Hartree ( %e eV )" % (i+1,npke, energy, energy * AtomicData.hartree_to_eV) )
            pad[i,:,:] = self.muffin.photoelectron_distribution(energy, grid, orbs,
                                                                pol=pol)

        plot_pads(pke, pad, pol)
        save_pads(pke, pad, pol, pad_file, units=units)

        # save intermediate variables for locating resonances
        self._pke = pke
        self._pad = pad
        self._pol = pol
        self._grid = grid
        self._orbs = orbs
        
        return pad

    def find_resonances(self, sigma_thresh=1.0):
        """

        identify resonances as peaks in the photoionization cross section

        Resonances are highly peaked local maxima of sigma(E) which exceed a certain threshold.
        First the local maxima are identified in the curve sigma(E) that was calculated in the
        a previous call to `calculate_pads(...)`. The energetic positions of the maxima are
        refined by bisection. The kinetic energy grid use to compute sigma(E) has to be fine 
        enough to obtain the initial guesses 

        Parameters
        ----------
        sigma_thresh    :  float, optional
            Maxima in the photoionization cross section are considered to be resonances
            if they exceed a threshold, sigma > sigma_thresh (in magebarn)

        Returns
        -------
        resonances         :  dict
            `resonances[i]` contains a list of continuum orbitals at the resonances for ionization
            from initial orbital `i`. The energy of the resonance can be accessed `resonances[i].energy`.

        """
        assert hasattr(self, "_pke"), "`find_resonances(...)` must be preceded by call to `calculate_pads(...)`."
        # retrieve data from previous PAD calculation
        pke = self._pke
        pad = self._pad
        pol = self._pol
        grid = self._grid
        orbs = self._orbs
        
        print( " " )
        print( "     **********************" )
        print( "     *      Resonances    *" )
        print( "     **********************" )
        print( " " )
        npke, norb, dummy = pad.shape

        # `energies[i]` is a list of photoelectron kinetic energy at resonance
        # for ionization from orbital `i`
        energies = {}
        # `sigmas[i]` contains list of values of sigma at resonances
        sigmas = {}
        # `resonances[i] contains list of continuum orbitals at resonances
        # (instances of `CMSWavefunction`)
        resonances = {}
        
        for i in range(0, norb):
            
            # Instead of maximizing sigma we minimize (-1)*sigma.
            
            # find indices of local minima
            minima = signal.argrelmin(-pad[:,i,0])[0].tolist()
            # and indices of local maxima
            maxima = signal.argrelmax(-pad[:,i,0])[0].tolist()

            if len(minima) == 0:
                # No local maximum of sigma, which is a local minimum of (-1)*sigma,
                # so no resonance
                continue
            if len(maxima) == 0:
                # No maximum was found, bracket minimum by end points
                maxima += [0,-1]
            # Each local minimum should be bracketed by two local maxima
            if pke[minima[0]] < pke[maxima[0]]:
                # first extremum is a minimum, so
                #   maxima[j-1] < minima[j] < maxima[j]
                maxima = [0] + maxima
                # After prepending the first point, we have
                #   maxima[j  ] < minima[j] < maxima[j+1]
            if pke[minima[-1]] > pke[maxima[-1]]:
                # last extremum is a minimum, which is not bracketed
                # by two maxima
                maxima = maxima + [-1]
                # After appending last point, we have
                #  maxima[i  ] < minima[i] < maxima[i+1]
                # for all minima
            assert len(minima) == len(maxima)-1
                        
            def func(energy):
                # compute (-1) x photoionization cross section for initial orbital
                # with index `i` at energy `energy`
                pad_i = self.muffin.photoelectron_distribution(energy, grid, orbs[i:i+1,:],
                                                               pol=pol)
                sigma_i = pad_i[0,0]
                return (-1)*sigma_i

            # list of photoelectron kinetic energy at resonance
            energies[i] = []
            # values of sigma at resonances
            sigmas[i] = []
            # list of continuum orbitals at resonances (instances of `CMSWavefunction`)
            resonances[i] = []
            
            # Refine energy at each local minimum of func(E) (= maximum of sigma)
            for j in range(0, len(minima)):
                # initial guess
                emin0 = pke[minima[j]]
                # maxima that bracket this local minimum
                emax_lower = pke[maxima[j]  ]
                emax_upper = pke[maxima[j+1]]
                assert emax_lower < emin0 < emax_upper
                # We search for a minimum around emin0 by Golden search
                # (https://en.wikipedia.org/wiki/Golden-section_search)
                # which assumes that there is a single minimum in the interval [l,u]
                alpha = 0.2
                l = (1.0-alpha)*emin0 + alpha*emax_lower
                u = (1.0-alpha)*emin0 + alpha*emax_upper
                # find minimum of func(E) = -log(cond(M(E))) in the interval [l,u]
                try:
                    emin = minimize_golden(func, l, u)
                except StopIteration:
                    continue
                fmin = func(emin)
                assert self.muffin.energy == emin

                sigma_max = -fmin
                if (sigma_max < sigma_thresh):
                    # sigma at maximum is too small to classify as a resonance
                    continue
                
                resonances[i] += self.muffin.eigenchannel_analysis()

                energies[i].append(emin)
                sigmas[i].append(sigma_max)

        if len(resonances.keys()) > 0:
            print( " -----------------------------------------------------------------------------" )
            print( "   Orbital   Resonance                 Energy                       Sigma     " )
            print( "                             Hartree             eV                  Mb       " )
            print( " -----------------------------------------------------------------------------" )
        else:
            print( " no resonances found with sigma > %e Mb" % sigma_thresh )
        for i in sorted(resonances.keys()):
            for j,res in enumerate(resonances[i]):
                print( "  %4.1d       %4.1d        %12.8f        %12.8f        %6.4e" % 
                   (i+1, j+1,
                    res.energy,
                    res.energy * AtomicData.hartree_to_eV,
                    res.tags["sigma"] * AtomicData.bohr2_to_megabarn) )
        print( "" )

        return resonances
                
    
def save_pads(pke,pad, pol, tbl_file, units="eV-Mb"):
    """
    A table with the PAD is written to `tbl_file`. 
    It contains the 4 columns   PKE   SIGMA  BETA1   BETA_2
    which define the PAD(th) at each energy according to

    .. code-block:: none
                                   
      PAD(th) = SIMGA/(4pi) [ 1 + BETA  P (cos(th)) + BETA  P (cos(th)) ]
                                      1  1                2  2

    For each orbital a block separated by a newline is written.
    """
    npke,norb,dummy = pad.shape
    sigma = pad[:,:,0]
    beta1 = pad[:,:,1]
    beta2 = pad[:,:,2]

    fh = open(tbl_file, "w")
    
    pol2str = {0 : "0 (linear)", -1 : "-1 (left)", +1: "+1 (right)"}

    print( """
#
# photoelectron angular distributions (PAD) for an isotropic ensemble
#
#  PAD(th) = SIMGA/(4pi) [ 1 + BETA  P (cos(th)) + BETA  P (cos(th)) ]
#                                  1  1                2  2
#
# light polarization = %s
#
    """ % pol2str[pol], file=fh)
    
    if units == "eV-Mb":
        pke = pke * AtomicData.hartree_to_eV
        # convert cross section sigma from bohr^2 to Mb
        sigma = sigma * AtomicData.bohr2_to_megabarn

        header = "# PKE/eV       SIGMA/Mb      BETA1     BETA2"
    else:
        header = "# PKE/Hartree  SIGMA/bohr^2  BETA1     BETA2"
    
    for b in range(0, norb):
        print( "# photoionization from orbital %d" % (b+1), file=fh)
        print( header, file=fh)
        block = np.vstack((pke, sigma[:,b], beta1[:,b], beta2[:,b])).transpose()
        np.savetxt(fh, block, fmt=["%e","%e","%+10.7f", "%+10.7f"])
        print( "" , file=fh)

    print( "PAD written to '%s'" % tbl_file )
    fh.close()

def plot_pads(pke,pad, pol, units='eV-Mb'):
    """
    plot PAD parameters sigma, beta1 and beta2 as functions of PKE
    for different orbitals
    """
    npke,norb,dummy = pad.shape
    sigma = pad[:,:,0]
    beta1 = pad[:,:,1]
    beta2 = pad[:,:,2]
    
    import matplotlib
    matplotlib.rc('xtick', labelsize=17)
    matplotlib.rc('ytick', labelsize=17)
    matplotlib.rc('legend', fontsize=17)
    matplotlib.rc('axes', labelsize=17)
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1,3)
    plt.title("polarization = %s" % pol)
    
    if units == "eV-Mb":
        pke = pke * AtomicData.hartree_to_eV
        # convert cross section sigma from bohr^2 to Mb
        sigma = sigma * AtomicData.bohr2_to_megabarn

        for ax in [0,1,2]:
            axes[ax].set_xlabel("PKE / eV")
        axes[0].set_ylabel(r"$\sigma$ / Mb")
    else:
        for ax in [0,1,2]:
            axes[ax].set_xlabel("PKE / Hartree")
        axes[0].set_ylabel(r"$\sigma$ / bohr$^2$")

    # plot sigma in log-scale
    axes[0].set_yscale("log")
    
    axes[1].set_ylabel(r"$\beta_1$")
    axes[2].set_ylabel(r"$\beta_2$")
    
    axes[2].set_ylim((-1.1,2.1))
    
    for b in range(0, norb):
        l, = axes[0].plot(pke, sigma[:,b], lw=2)
        axes[1].plot(pke, beta1[:,b], lw=2, color=l.get_color(), label=r"Orb. %d" % (b+1))
        axes[2].plot(pke, beta2[:,b], lw=2, color=l.get_color())

    axes[1].legend(loc="upper center")

    plt.subplots_adjust(wspace=0.5)
    
    plt.show()
    

##################################################
#                                                #
# Testing                                        #
#                                                #
##################################################

def test_hmi_pads():
    """
    compute PADs for ionization from 1sigma orbital of the hydrogen molecular ion (H2^+)
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
                                lmax=10,
                                charge=+1, chargeIII=+2,
                                potential_type="molecular", Nr=2000, debug=0)

    # We search for bound state in the energy range [-1.2, -0.5] (in Hartree),
    # since we only want to find the 1\sigma_g and 1\sigma_u orbitals.
    search_energies = np.linspace(-1.2, -0.5, 20)
    bound_orbitals = muffin.find_eigenstates(search_energies)

    pad = PhotoelectronAngularDistribution(muffin)

    pke_energies = np.linspace(0.1,400.0,100) / AtomicData.hartree_to_eV
    pad.compute_pads(bound_orbitals, pke_energies,
                     pol=0,
                     pad_file="/tmp/pad.dat")

if __name__ == "__main__":
    test_hmi_pads()
