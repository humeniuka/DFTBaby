#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
classes for representing multiple scattering wavefunctions of bound
and continuum orbitals
"""

from DFTB.MolecularIntegrals.Ints1e import integral
# The modules `ms.so` and `cms.so` are wrappers around Fortran code.
from DFTB.MultipleScattering import ms
from DFTB.MultipleScattering import cms

import numpy as np


class MSWavefunction(object):
    """
    multiple scattering (MS) wavefunction for bound orbital
    """
    # Additional properties of an orbital such as symmetry labels, etc.
    # can be stored in the `tags` dictionary.
    tags = {}
    def __init__(self, energy, Vconst, chargeIII, rhoIII, radii, atpos, kspl,
                 fIknots, fIcoeffs, fIIIknots, fIIIcoeffs,
                 coeffs,
                 singval):
        """
        create MS wavefunction

        This class encapsulates all data necessary for constructing the wavefunction,
        such as atomic coordinates, radial wavefunctions in region I and III and
        their coefficients. 

        The wavefunctions are only approximately normalized to 1, since the contribution
        from the interstitial region II is neglected.

        Parameters
        ----------
        parameters are the same as for `ms.wavefunctions(...)` except for
        singval     :  eigenvalue belonging to the coefficients coeffsA, should be close to 0
        """
        # We have to make copies of numpy arrays, since these might be changed
        # outside this function
        self.energy = energy
        self.Vconst = Vconst
        self.chargeIII = chargeIII
        self.rhoIII = rhoIII
        self.radii = np.copy(radii)
        self.atpos = np.copy(atpos)
        self.kspl = kspl
        self.fIknots = np.copy(fIknots)
        self.fIcoeffs = np.copy(fIcoeffs)
        self.fIIIknots = np.copy(fIIIknots)
        self.fIIIcoeffs = np.copy(fIIIcoeffs)
        self.coeffs = np.copy(coeffs)
        # singular value
        self.singval = singval
        # norm of the orbital, is used to normalize to wavefunction
        # to 1
        self.norm = 1.0
    def __call__(self, x,y,z):
        """
        evaluate wavefunction on a grid

        Parameters
        ----------
        x,y,z      : 1d numpy arrays with coordinates of grid points
                     The cartesian position of the i-th point is given by
                     x[i], y[i], z[i]

        Returns
        -------
        wfn        : complex numpy array, wfn[i] is the wavefunction
                     at the point (x[i],y[i],z[i])
        """
        # The coordinate arrays may have any shape, but the
        # fortran code can only deal with rank 1 arrays
        shape = x.shape
        x,y,z = x.flatten(), y.flatten(), z.flatten()

        # coeffs is a vector with coefficients A^{II0}, A^{IIi}
        wfn = ms.wavefunctions(self.energy, self.Vconst, self.chargeIII, self.rhoIII, self.radii,
                               self.atpos, self.kspl, self.fIknots, self.fIcoeffs, self.fIIIknots, self.fIIIcoeffs,
                               self.coeffs,
                               x,y,z)
        # divide by normalization constant
        wfn /= self.norm

        wfn = wfn[0,:]

        # Convert wavefunction array to the same shape as
        # the input arrays
        wfn = np.reshape(wfn, shape)
        
        return wfn
    
    def normalize(self, atomlist):
        """
        compute and set the normalization constant so that the
        orbital is normalized to 1.

        Parameters
        ----------
        atomlist   :  molecule geometry for defining Becke's grid 
                      for numerical integration numerical 

        Returns
        -------
        norm       :  norm of wavefunction <wfn|wfn>^{1/2} 
                      before normalization
        """
        # norm^2 of wavefunction <wfn|wfn>
        nrm2 = integral(atomlist,
                        lambda x,y,z: abs(self.__call__(x,y,z))**2)
        nrm = np.sqrt(nrm2)
        # update normalization constant
        self.norm *= nrm

        return nrm


class CMSWavefunction(object):
    """
    continuum multiple scattering (CMS) wavefunction for unbound orbital
    """
    # Additional properties of an orbital such as symmetry labels, phase shifts etc.
    # can be stored in the `tags` dictionary.
    tags = {}
    def __init__(self, energy, Vconst, chargeIII, rhoIII, radii, atpos, kspl,
                 fIknots, fIcoeffs, sol,
                 coeffs):
        """
        create CMS wavefunction

        This class encapsulates all data necessary for constructing a continuum wavefunction,
        such as atomic coordinates, radial wavefunctions in region I and their coefficients. 

        The wavefunction is a linear combination of K-matrix normalized standing waves.

             wfn(x,y,z) = sum c  Psi(x,y,z)
                           L   L    L
        Parameters
        ----------
        parameters are the same as for `cms.wavefunctions(...)` except for `coeffs`
        coeffs      :  coefficients `c` for linear combination
        """
        # We have to make copies of numpy arrays, since these might be changed
        # outside this function
        self.energy = energy
        self.Vconst = Vconst
        self.chargeIII = chargeIII
        self.rhoIII = rhoIII
        self.radii = np.copy(radii)
        self.atpos = np.copy(atpos)
        self.kspl = kspl
        self.fIknots = np.copy(fIknots)
        self.fIcoeffs = np.copy(fIcoeffs)
        # solution vector from matching for all L=(l,m) channels
        self.sol = sol
        self.coeffs = np.copy(coeffs)
    def __call__(self, x,y,z):
        """
        evaluate wavefunction on a grid

        Parameters
        ----------
        x,y,z      : 1d numpy arrays with coordinates of grid points
                     The cartesian position of the i-th point is given by
                     x[i], y[i], z[i]

        Returns
        -------
        wfn        : complex numpy array, wfn[i] is the wavefunction
                     at the point (x[i],y[i],z[i])
        """
        # The coordinate arrays may have any shape, but the
        # fortran code can only deal with rank 1 arrays
        shape = x.shape
        x,y,z = x.flatten(), y.flatten(), z.flatten()
                                                        
        kmat, wfn_standing = cms.wavefunctions(self.energy, self.Vconst, self.chargeIII, self.rhoIII, self.radii,
                                               self.atpos, self.kspl, self.fIknots, self.fIcoeffs, self.sol, 0,
                                               x,y,z)
        # form linear combination
        wfn = np.dot(self.coeffs, wfn_standing)

        # Convert wavefunction array to the same shape as
        # the input arrays
        wfn = np.reshape(wfn, shape)
                                            
        return wfn
