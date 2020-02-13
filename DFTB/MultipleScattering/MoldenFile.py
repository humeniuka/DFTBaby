#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Module for reading geometry, spherical GTO basis and MO coefficients from file in molden format.

"""
from __future__ import print_function

from DFTB.MolecularIntegrals.SphericalCoords import cartesian2spherical
from DFTB.MolecularIntegrals.LebedevQuadrature import spherical_harmonics_it
from DFTB import AtomicData

import numpy as np
from scipy import special
from collections import OrderedDict
import re
import itertools

letter2angmom = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6}

class ContractedGaussianOrbital(object):
    """

    Each atomic orbital is the product of a real spherical harmonic function Y~_{l,m}
    and a radial function, which is a contraction of primitive Gaussian functions
    with coefficients c_i, exponents a_i and normalization constants N(a_i,l):

    .. code-block:: none
    
                      l      i=n              -a_i*r^2   ~
      AO(r,th,ph) =  r  [ sum    N(a ,l) c   e         ] Y (th,ph)
                             i=1    i     i               l,m

                 2*(2a)^{l+3/2}  1/2
      N(a,l) = [ -------------- ]
                  Gamma(l+3/2)


    Parameters
    ----------
    center       : tuple (x,y,z)
                   center of atomic basis function (in bohr)
    iatom        : int
                   index into list of atoms, indicates which atom this basis function belongs to
    ifunc        : int
                   index into list of basis functions, used for assigning MO coefficients to AOs
    ishell       : int
                   index into atomic shells, basis functions in the same shell share the same value
                   of `iatom` and `l`.
    l            : int
                   azimuthal quantum number
    m            : int, -l <= m <= l
                   magnetic quantum number
    exponents    : np.1darray (shape (n,))
                   exponents of Gaussian
    coefficients : np.1darray (shape (n,))
                   contraction coefficients

    """

    def __init__(self, center, iatom, ifunc, ishell, l, m, exponents, coefficients):
        assert -l <= m <= l
        assert len(exponents) == len(coefficients), "Number of exponents and coefficients not equal!"
        self.center = center
        self.iatom = iatom
        self.ifunc = ifunc
        self.ishell = ishell
        self.l = l
        self.m = m
        self.exponents = exponents
        self.coefficients = coefficients
        # normalization constants
        self.norms = np.sqrt(2*(2*exponents)**(l+1.5)/special.gamma(l+1.5))
    def __call__(self, x,y,z):
        """
        evaluate basis function on a grid

        Parameters
        ----------
        x,y,z      : np.1darray
                     coordinates of grid points. The cartesian position of the 
                     i-th point is given by `x[i]`, `y[i]`, `z[i]`.

        Returns
        -------
        wfn        : np.1darray
                     wfn[i] is the atomic orbital at the point `(x[i],y[i],z[i])`

        """
        x0,y0,z0 = self.center
        r,th,ph = cartesian2spherical((x-x0,y-y0,z-z0))
        
        # radial part of Gaussian type orbital
        wfn_radial = 0*x
        for a,c,N in zip(self.exponents, self.coefficients, self.norms):
            wfn_radial += N*c*np.exp(-a*r**2)
        wfn_radial *= r**self.l
        
        # angular part of wavefunction, Y^(real)_{l,m}(th,ph)
        sph_it = spherical_harmonics_it(th,ph)
        for Ylm,ll,mm in sph_it:
            # This is very inefficient, since all spherical harmonics with ll<l, mm<m have
            # to be generate first. It is much better to group basis functions on the same
            # atom together and evaluate them in one go.
            if ll == self.l and mm == self.m:
                # real spherical harmonics
                if mm < 0:
                    Ylm_real = -np.sqrt(2.0) * Ylm.imag
                elif mm > 0:
                    Ylm_real =  np.sqrt(2.0) * (-1)**mm * Ylm.real
                else:
                    # mm == 0
                    Ylm_real = Ylm.real
                    
                wfn_angular = Ylm_real
                break

        # combine radial and angular parts
        wfn = wfn_radial * wfn_angular
        
        return wfn

class MolecularOrbital(object):
    """

    LCAO molecular orbital
    
    Parameters
    ----------
    atomlist     : list of tuples `(atomic_number, (x,y,z))`
       molecular geometry in cartesian coordinates (bohr)
    coefficients : np.1darray (shape (nao,))
       MO coefficients with respect to AOs in `basis`
    basis        : list of `ContractedGaussianOrbital`
       basis functions
    tags         : dict, optional
       dictionary containing additiona information about the orbital such as:

          * energy       : float
            MO energy (in Hartree)
          * spin         : str
            spin of orbital ('Alpha' or 'Beta')
          * symmetry     : str
            symmetry or other label characterizing this orbital
          * occup        : float
            occupation or norm of the orbital

    """
    def __init__(self, atomlist, coefficients, basis, tags={}):
        self.atomlist = atomlist
        self.coefficients = coefficients
        # The basis functions are group by atoms and angular momentum shell
        # to be able to reuse intermediate results for all basis function
        # belonging to the same atom or shell.
        self.basis = sorted(basis, key=lambda cgbf: (cgbf.iatom, cgbf.l, abs(cgbf.m), -np.sign(cgbf.m)) )
        self.tags = tags
    def __call__(self, x,y,z):
        """

        evaluate molecular orbital on a grid

        Parameters
        ----------
        x,y,z      : np.1darray
                     coordinates of grid points. The cartesian position of the 
                     i-th point is given by `x[i]`, `y[i]`, `z[i]`.

        Returns
        -------
        wfn        : np.1darray
                     wfn[i] is the molecular orbital at the point `(x[i],y[i],z[i])`

        """
        wfn = 0.0*x
        
        # The following commented section is a slow and stupid way to evaluate the MO,
        # because many primitive Gaussians share the same spherical harmonics.
        # Below is a faster implementation that gives exactly the same result.
        #
        #for cgbf in self.basis:
        #    wfn += self.coefficients[cgbf.ifunc] * cgbf(x,y,z)
        #return wfn
        #

        # group basis functions by atomic center
        atom_group = itertools.groupby(self.basis, key=lambda cgbf: cgbf.iatom)
        for iatom, atomic_basis in atom_group:
            # spherical coordinates around atomic center
            x0, y0, z0 = self.atomlist[iatom][1]
            r,th,ph = cartesian2spherical((x-x0,y-y0,z-z0))
            # spherical harmonics can be shared by all basis functions having the same l,m
            sph_it = spherical_harmonics_it(th,ph)
            # group basis functions on atom by angular momentum quantum numbers (l,m)
            # in the same order in which the spherical harmonics are generated iteratively
            angmom_group = itertools.groupby(atomic_basis, key=lambda cgbf: (cgbf.l, abs(cgbf.m), (-1)*np.sign(cgbf.m)))
            for (l,abs_m,msign_m), shell_basis in angmom_group:
                m = (-1) * msign_m * abs_m

                Ylm,ll,mm = sph_it.next()
                assert ll == l and mm == m
                # real spherical harmonics
                if mm < 0:
                    Ylm_real = -np.sqrt(2.0) * Ylm.imag
                elif mm > 0:
                    Ylm_real =  np.sqrt(2.0) * (-1)**mm * Ylm.real
                else:
                    # mm == 0
                    Ylm_real = Ylm.real
                # radial part is the same for all basis functions in `shell_basis`    
                wfn_angular = Ylm_real
                # radial parts differ
                for cgbf in shell_basis:
                    assert cgbf.l == l and cgbf.m == m
                    # radial part of Gaussian type orbital
                    wfn_radial = 0*x
                    for a,c,N in zip(cgbf.exponents, cgbf.coefficients, cgbf.norms):
                        wfn_radial += N*c*np.exp(-a*r**2)
                    wfn_radial *= r**l

                    wfn += self.coefficients[cgbf.ifunc] * wfn_radial * wfn_angular
        return wfn
                

class MoldenFileInput(object):
    """

    reads geometry, Gaussian basis functions and molecular orbital coefficients from
    a file in the Molden format (see http://cheminf.cmbi.ru.nl/molden/molden_format.html)

    Parameters
    ----------
    moldenfile   :  str
       path to molden file

    """
    atomlist = []
    """list of tuples `(atomic_number, (x,y,z))` with cartesian coordinates (in bohr)"""
    basis = []
    """list of basis functions (instances of class `ContractedGaussianOrbital`)"""
    mos = []
    """list of molecular orbitals (instances of class `MolecularOrbital`)"""
    def __init__(self, moldenfile):
        # read entire molden file
        with open(moldenfile) as fh:
            self.lines = fh.readlines()
        # Matches section headers such as [GTO], [MO], etc.
        regex = re.compile("\[([\ 0-9a-zA-Z]+)\]\s*([\ 0-9a-zA-Z]*)")
        # Split file into different sections, each section can have options associated with
        # it, such as Angs|AU for the [Atoms] section
        self.sections = OrderedDict()
        self.options = OrderedDict()
        curr_section = []
        for il,l in enumerate(self.lines):
            m = regex.search(l.strip())
            if m:
                title = m.group(1).lower()
                option = m.group(2).lower()
                # start a new section with `title`
                curr_section = []
                self.sections[title] = curr_section
                self.options[title] = option
            else:
                # add line to current section
                curr_section.append(l.strip())

        # [Atoms] section is required
        self.read_atoms(self.sections["atoms"], units=self.options.get("atoms", "angs"))
        
        for title, lines in self.sections.iteritems():
            if title == "gto":
                self.read_gto_basis(lines)
            elif title == "mo":
                self.read_mos(lines)
            else:
                pass
                #raise Warning("No reader implemented for section '%s'" % title)
                
    def read_atoms(self, lines, units="angs"):
        """
        parse [Atoms] section, list of tuples (atomic_number, (x,y,z)) for each atom 
        is stored in `self.atomlist`
        """
        #print( "reading atoms ... ")
        self.atomlist = []
        for l in lines:
            words = l.split()
            # atomic number
            Zat = int(words[2])
            # cartesian coordinates
            if units == "au":
                # coordinates are in bohr, already
                c = 1.0
            else:
                # coordinates are in Angstrom, we have to convert them to bohr
                c = 1.0/AtomicData.bohr_to_angs
            pos = map(lambda xyz: float(xyz) * c, words[3:6])
            assert AtomicData.atom_names[Zat-1] == words[0].lower()
            self.atomlist.append( (Zat, pos) )
            
    def read_gto_basis(self, lines):
        """
        parse [GTO] section, list of contracted Gaussian basis functions is stored in `self.basis`
        """
        #print( "reading GTO basis ... ")
        if not (   self.sections.has_key("5d")
                or self.sections.has_key("7f")
                or self.sections.has_key("9g")):
            raise RuntimeError("Only spherical basis (5d, 7f, 9g) are supported!")

        self.basis = []

        nat = len(self.atomlist)
        # count lines
        k = -1
        # count basis functions
        ifunc = 0
        # loop over atoms
        for iatom in range(0, nat):
            k += 1
            # in Molden atom indices are 1-based
            assert int(lines[k]) == iatom+1
            # count atomic shells
            ishell = 0
            # loop over basis functions for atom iatom
            while (k < len(lines)):
                k += 1
                line = lines[k].strip()
                if line == "":
                    # newline marks end of atomic basis
                    break
                words = line.split()
                l, nprim = letter2angmom[words[0].lower()], int(words[1])
                # read exponents and contraction coefficients for `nprim` Gaussians
                exponents    = np.zeros(nprim)
                coefficients = np.zeros(nprim)
                for p in range(0, nprim):
                    k += 1
                    words = lines[k].split()
                    exponents[p] = float(words[0])
                    coefficients[p] = float(words[1])

                center = self.atomlist[iatom][1]

                # create basis functions for l-shell in the correct order
                if l == 1:
                    # p-orbitals are sorted in the order px,py,pz (i.e. m=+1,-1,0)
                    for m in [1,-1,0]:
                        self.basis.append( ContractedGaussianOrbital(center, iatom, ifunc, ishell, 
                                                                 l, m, exponents, coefficients) )
                        ifunc += 1                        
                else:
                    # All other orbitals (l=0 or l > 1) are sorted in the order
                    #  m=0, m=+1, m=-1, m=+2, m=-2, ..., m=+l, m=-l
                    for m in range(0, l+1):
                        self.basis.append( ContractedGaussianOrbital(center, iatom, ifunc, ishell, 
                                                                 l, m, exponents, coefficients) )
                        ifunc += 1
                        if m > 0:
                            self.basis.append( ContractedGaussianOrbital(center, iatom, ifunc, ishell, 
                                                                     l, -m, exponents, coefficients) )
                            ifunc += 1
                ishell += 1

    def read_mos(self, lines):
        """
        parse [MO] section
        """
        #print( "reading MO coefficients ... ")
        # count lines
        k = 0
        # enumerate orbitals
        iorb = 0
        while (k < len(lines)):
            # read energy, spin, symmetry and occupation of orbital
            #   Ene=    0.395606
            #   Spin= Alpha
            #   Sym= 0->1
            #   Occup=     0.852538
            # key-value pairs are stored in the `tags` dictionary
            tags = {}
            line = lines[k]
            while "=" in line:
                key, value = line.split("=")
                if key == "Ene":
                    tags["energy"] = float(value)
                elif key == "Spin":
                    tags["spin"] = value
                elif key == "Sym":
                    tags["symmetry"] = value
                elif key == "Occup":
                    tags["occup"] = float(value)
                else:
                    raise RuntimeError("Unknown key '%s' in [MO] section" % key)
                k += 1
                line = lines[k]
            # read MO coefficients
            nbfs = len(self.basis)
            coefficients = np.zeros(nbfs)
            for k in range(k, k+nbfs):
                words = lines[k].split()
                iao = int(words[0])-1
                assert iao < nbfs, "Number of MO coefficients is larger than number of basis functions (%d), maybe you are using cartesian basis functions instead of spherical ones?" % nbfs
                coefficients[iao] = float(words[1])
                
            self.mos.append( MolecularOrbital(self.atomlist, coefficients, self.basis, tags) )
            iorb += 1
            k += 1

    # access functions
    def getMOs(self):
        """

        Returns
        -------
        mos      :  list of `MolecularOrbital`
            molecular orbitals

        """
        return self.mos
    
    def getGeometry(self):
        """

        Returns
        -------
        atomlist :  list of tuples `(atomic_number, (x,y,z))`
            cartesian coordinates of molecular geometry

        """
        return self.atomlist
            
##########################
#                        #
#       Testing          #
#                        #
##########################

def _check_contraction_normalization():
    """
    check that contracted Gaussian orbital is normalized to 1.
    """
    from DFTB.MolecularIntegrals.Ints1e import overlap

    # D function for carbon from ANO-RCC basis set
    center = np.array([0.0, 0.0, 0.0])
    exponents = np.array([2.0105333, .79109810, .31127870, .12248090, .04899240])
    coefficients = np.array([-0.05130071, 0.33215275, 0.70352371, -0.83718957, -0.33535356])
    l = 2
    m = -1
                  
    cgbf = ContractedGaussianOrbital(center, 0, 0, 0, l, m, exponents, coefficients)

    atomlist = [(6, center)]
    S = overlap(atomlist, cgbf, cgbf)
    print( "overlap = %s" % S )

def _check_mo_normalization(moldenfile):
    """
    check that MOs are properly normalized
    """
    from DFTB.MolecularIntegrals.Ints1e import overlap
    
    mf = MoldenFileInput(moldenfile)
    atomlist = mf.getGeometry()
    mos = mf.getMOs()
    
    nmo = len(mos)
    for i in range(0, nmo):
        nrm = np.sqrt( overlap(atomlist, mos[i], mos[i]) )
        print( "norm of orbital %d  = %e" % (i,nrm) )
    
    
if __name__ == "__main__":
    _check_contraction_normalization()
    #_check_mo_normalization("/local_1tb/net/storage/humeniuka/data/BAGEL_TEST/h2o+_rohf.molden")

    
