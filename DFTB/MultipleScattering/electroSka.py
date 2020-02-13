#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

**********
electroSka
**********

===========
Description
===========
photoangular distributions for isotropically oriented ensembles 
of molecules in the gas phase

Initial bound and final continuum orbitals are obtained by the 
multiple scattering (MS) and the continuum multiple scattering 
(CMS) methods, respectively. Resonances in the continuum can
also be identified.

===============
Format of Input
===============
The input is specified in a JSON file. For an example input file see below.
The calculation is started with:

.. code-block:: bash

   electroSka.py  input.json



**Keywords:**

   | ``lmax``
   |
   | **Description:** Maximum angular momentum l of spherical wave expansion in regions I and III.
   | **Datatype:** int
   | **Default:** 10
   |
   | ``potential_type``
   |
   | **Description:** The potential around each atom in region I is spherically symmetric.
   |     The radial potential in each atomic sphere is either obtained as the
   |     DFT potential of a free atom ('atomic') or by spherically averaging the 
   |     superpositions of atomic DFT potentials around each center ('molecular'). 
   |     'molecular' contains the contribution of one atom to the potential of a
   |     neighbouring atom. 'atomic' allows to reuse the same potential
   |     and wavefunctions for all atoms of the same type.
   | **Datatype:** str
   |       'atomic'
   |       'molecular'
   | **Default:** 'molecular'
   |
   | ``debug``
   |
   | **Description:** Control amount of output
   | **Datatype:** int
   |       -1 (terse)
   |       0  (normal)
   |       >0 (lots of additional output)
   | **Default:** 0
   |  

**Subsections:**

.. topic:: ``molecule``

   | **Description:** Molecular geometry and charge.
   | **Keywords:**
   | 
   |            ``geometry``
   | 
   |            **Description:** List of element symbols and cartesian coordinates.
   |            **Datatype:** List of atoms in the format :code:`{"atom": "symbol", "xyz" : [ x, y, z ] }`
   |            **Default:** no default, has to be specified
   |
   |            ``charge``
   |
   |            **Description:** Total charge of molecule
   |            **Datatype:** int
   |            **Default:** 0
   |
   |            ``units``
   |
   |            **Description:** Units of coordinates
   |            **Datatype:** str
   |            **Default:** 'Angstrom', anything else means bohr.
   |

.. topic:: ``becke_grid``
  
   | **Description:** Multicenter Becke grid for numerical integration.
   | **Keywords:**
   |
   |            ``radial_grid_factor``
   |
   |            **Description:** The number of radial grid point for each atomic grid is increased
   |                by this factor.
   |            **Datatype:** int
   |            **Default:** 10
   |
   |            ``lebedev_order``
   |
   |            **Description:** Order of Lebedev grid for angular integration.
   |            **Datatype:** int
   |            **Default:** 23
   |

.. topic:: ``bound``

   | **Description:** Controls calculation of bound orbitals (E < 0).
   | **Keywords:**
   |
   |            ``search_energies``
   |
   |            **Description:** List of guess energies. Eigenenergies are searched for using the 
   |                shooting method. The energy grid has to be fine enough, so that each interval 
   |                contains at most one eigenenergy. 
   |            **Datatype:** vector<float> or python expression that generates a 1d numpy array
   |                such as 'linspace(-5.0, -0.01, 500)'.
   |            **Default:** no default, has to be specified
   |
   |            ``units``
   |
   |            **Description:** Units for energy
   |            **Datatype:** str ('eV' or 'Hartree')
   |            **Default:** 'Hartree'
   |

.. topic:: ``dyson``

   | **Description:** Read Dyson orbitals from a molden file and use them as initial orbitals 
   |     instead of the bound orbitals if ``orbitals`` is set to 'dyson' in the ``pad`` section.
   | **Keywords:**
   |
   |            ``molden_file``
   |
   |            **Description:** Path to a file in the molden format with molecular orbitals.
   |            **Datatype:** str
   |            **Default:** 'dyson_orbitals.molden'
   |
   |            ``project``
   |
   |            **Description:** Project the Dyson orbitals into the subspace spanned by all bound
   |                orbitals.
   |            **Datatype:** bool
   |            **Default:** false
   |

.. topic:: ``continuum``

   | **Description:** Controls calculation of continuum orbitals (E > 0).
   | **Keywords:**
   |
   |            ``kinetic_energies``
   |
   |            **Description:** Kinetic energies of photoelectron. For each energy in the list
   |                continuum orbitals are obtained using the CMS method and PADs are calculated.
   |            **Datatype:** vector<float> or python expression that generates a 1d numpy array
   |                such as 'linspace(0.1, 5.0, 100)'.
   |            **Default:** no default, has to be specified
   |
   |            ``units``
   |
   |            **Description:** Units for kinetic energy
   |            **Datatype:** str ('eV' or 'Hartree')
   |            **Default:** 'Hartree'
   |
   |            ``resonances``
   |
   |            **Description:** Search for resonances, i.e. for peaks in photoionization cross 
   |                section sigma. After identifying local maxima of sigma(E) among the points
   |                in ``kinetic_energies`` the maxima a optimized by bisection. Only maxima E* where
   |                sigma(E*) exceeds the threshold ``sigma_thresh`` are considered to be resonances.
   |            **Datatype:** bool
   |            **Default:** false
   |
   |            ``sigma_thresh``
   |
   |            **Description:** Threshold (in Mb) for accepting a maximum of sigma(E) as a resonance.
   |            **Datatype:** float
   |            **Default:** 10.0
   |

.. topic:: ``pad``

   | **Description:** Controls computation of photoelectron angular distributions.
   | **Keywords:**
   |
   |            ``orbitals``
   |           
   |            **Description:** Source of the initial orbitals.
   |            **Datatype:** str
   |                   'bound' - use the bound orbitals as initial states
   |                   'dyson' - use Dyson orbitals, requires ``dyson`` section
   |            **Default:** 'bound'
   |
   |            ``initial``
   |           
   |            **Description:** Indices of initial states (starting from 1).
   |            **Datatype:** vector<int>
   |            **Default:** all bound states found in MS calculation
   |
   |            ``pol``
   |             
   |            **Description:** Polarization of ionizing electromagnetic field.
   |            **Datatype:** int
   |                   -1 (left)
   |                   0 (linear)
   |                   +1 (right)
   |            **Default:** 0
   |
   |            ``pad_file``
   |          
   |            **Description:** A table with the photoelectron angular distribution for
   |                ionization from all initial orbitals is written to this file.
   |            **Datatype:** str
   |            **Default:** 'pad.dat'
   |
   |            ``units``
   | 
   |            **Description:** Units of photokinetic energies and cross sections in PAD table.
   |
   |            **Datatype:** str
   |                   'eV-Mb' :   PKE in eV, sigma in megabarn
   |                   'au'    :   PKE in Hartree, sigma in bohr^2
   |            **Default:** 'eV-Mb'
   |

.. topic:: ``cubes``

   | **Description:** Bound orbitals and continuum orbitals at resonances can be saved to cube files.
   | **Keywords:**
   |
   |            ``prefix``
   |           
   |            **Description:** The prefix is prepended to any filename for saving cube files. 
   |                The real and imaginary parts are saved to separate cube files. Bound orbitals are 
   |                saved to 'PREFIX_bound_%d_real.cube' 'PREFIX_bound_%d_imag.cube' where %d is the 
   |                index of the orbital (starting from 1). Resonances in the continuum are saved 
   |                to 'PREFIX_resonance_%d_to_%d_real.cube' and 'PREFIX_resonance_%d_to_%d_imag.cube',
   |                where the first integer indicates the initial orbital and the second index enumerates
   |                the resonances.
   |            **Datatype:** str
   |            **Default:** '/tmp/ms\_'
   | 
   |            ``ppb``
   |
   |            **Description:** Resolution of grid in points per bohr.
   |            **Datatype:** float
   |            **Default:** 2.0
   |
   |            ``dbuff``
   |
   |            **Description:** Additional space around molecule in bohr
   |            **Datatype:** float
   |            **Default:** 5.0
   |
   |            ``export_bound``
   |
   |            **Description:** Should bound orbitals be saved?
   |            **Datatype:** bool
   |            **Default:** false
   |
   |            ``export_projected``
   |
   |            **Description:** Should projected Dyson orbitals be saved?
   |            **Datatype:** bool
   |            **Default:** false
   |
   |            ``export_resonances``
   |
   |            **Description:** Should continuum orbitals at resonances be saved?
   |            **Datatype:** bool
   |            **Default:** false
   |   

=======
Example
=======


Sample input
------------

.. code-block:: javascript

   {
    "electroSka": 	{

	"comment"         : "MS, CMS and PAD calculation at AM1 optimized geometry of water",
	
	"lmax"            :  8,
	"potential_type"  :  "molecular",
	
	"molecule"        : {
	    "units" :  "Angstrom",
	    "geometry" : [
		{"atom": "H", "xyz" : [  0.755131,   -0.475985,    0.000000] },
		{"atom": "O", "xyz" : [  0.000000,    0.118996,    0.000000] },
		{"atom": "H", "xyz" : [ -0.755131,   -0.475985,    0.000000] }
	    ],
	    "charge"   : 0	    
	},
	
	"becke_grid"      : {
	    "radial_grid_factor"   : 10,
	    "lebedev_order"        : 23
	},
	
	"bound"           : {
	    "units"            : "Hartree",
	    "search_energies"  : "linspace(-5.0, -0.01, 500)"
	},
	
	"continuum"       : {
	    "units"            : "Hartree",
	    "kinetic_energies" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            "resonances"       : true,
            "sigma_thresh"     : 10.0
	},
	
	"pad"             : {
	    "initial"         : [1,2,3,4],
	    "pol"             : 0,
            "units"           : "eV-Mb",
	    "pad_file"        : "/tmp/water.pad"
	},
	
	"cubes"           : {
	    "export_bound"      : true,
	    "export_resonances" : true,
	    "prefix"            : "/tmp/water",
	    "ppb"               : 3.0,
	    "dbuff"             : 5.0
	}
    }
   }


==========
References
==========

+----------------------------------------------------+-------------------------------------------------------------+
|          Description of Reference                  |                          Reference                          |
+====================================================+=============================================================+
| Multiple scattering theory (bound states)          | K.\ Johnson, Adv. Quantum Chem. **7**, 143-185 (1973)       |
+----------------------------------------------------+-------------------------------------------------------------+
| Continuum multiple scattering (unbound states)     | D.\ Dill and J.\ Dehmer, J. Chem. Phys. **61**, 692 (1974)  |
+----------------------------------------------------+-------------------------------------------------------------+
| photoelectron angular distributions                | B.\ Ritchie, Phys. Rev. A **13**, 1411-1415, (1976)         |
+----------------------------------------------------+-------------------------------------------------------------+

"""
from __future__ import print_function
from DFTB.MolecularIntegrals import settings
import json
import numpy as np

from DFTB import AtomicData
from DFTB.MultipleScattering.MuffinTin import MuffinTinPotential, save_cubefile, geometries_are_equal, project_orbital
from DFTB.MultipleScattering.PAD import PhotoelectronAngularDistribution
from DFTB.MultipleScattering.MoldenFile import MoldenFileInput

def get_geometry(mol):
    """
    extract atomlist from block with title 'molecule'
    """
    atoms = mol.get("geometry")
    atomlist = []
    for atom in atoms:
        atname = atom.get("atom")
        atnumber = AtomicData.atomic_number(str(atname))
        coords = atom.get("xyz")
        if mol.get("units", "Angstrom") == "Angstrom":
            coords = map(lambda xyz: xyz/AtomicData.bohr_to_angs, coords)
        atomlist.append( (atnumber, coords) )

    return atomlist


def run_from_config(config):
    """

    run MS, CMS and PAD calculation

    Parameters
    ----------
    config      : dict
        input for multiple scattering calculation

    """
    # import linspace into the namespace so that inside the json file
    # we can use 'linspace(...)'
    from numpy import linspace

    becke_opts = config.get("becke_grid", {})
    
    # choose resolution of multicenter grids
    #    controls size of radial grid 
    settings.radial_grid_factor = becke_opts.get("radial_grid_factor", 10)
    #    controls size of angular grid
    settings.lebedev_order = becke_opts.get("lebedev_order", 23)

    try:
        mol = config.get("molecule")
    except KeyError:
        raise KeyError("Required block 'molecule' with geometry and charge not found!")
    
    atomlist = get_geometry(mol)

    bound_opts     = config.get("bound", {})
    dyson_opts     = config.get("dyson", {})
    continuum_opts = config.get("continuum", {})
    pad_opts       = config.get("pad", {})
    cube_opts      = config.get("cubes", {})
    
    # photoelectron kinetic energies
    v = continuum_opts["kinetic_energies"]
    if type(v) == list:
        kinetic_energies = np.array(v)
    elif type(v) == unicode:
        kinetic_energies = eval(v)
    else:
        raise ValueError("'kinetic_energies' should be a list of energies or a python expression")
        
    if continuum_opts.get("units", "Hartree") == "eV":
        kinetic_energies /= AtomicData.hartree_to_eV

    # The overall charge of the HMI is +1, however the single electron will
    # feel the nuclear attraction from both protons, therefore chargeIII=+2.
    muffin = MuffinTinPotential(atomlist,
                                lmax=config.get("lmax", 10),
                                charge=mol.get("charge", 0),
                                chargeIII=mol.get("charge", 0)+1,
                                potential_type=config.get("potential_type", "molecular"),
                                debug=config.get("debug", 0))

    if bound_opts != {}:
        # Energy ranges where to search for bound states
        v = bound_opts["search_energies"]
        if type(v) == list:
            search_energies = np.array(v)
        elif type(v) == unicode:
            search_energies = eval(v)
        else:
            raise ValueError("'search_energies' should be a list of energies or a python expression")

        if bound_opts.get("units", "Hartree") == "eV":
            search_energies /= AtomicData.hartree_to_eV

        # compute bound orbitals
        bound_orbitals = muffin.find_eigenstates(search_energies)

        # save bound orbitals to cube files
        if cube_opts.get("export_bound", False) == True:
            print( "exporting bound orbitals to cube files" )
            prefix = cube_opts.get("prefix", "/tmp/ms_")
            for iorb, orbital in enumerate(bound_orbitals):
                filename = prefix + "_bound_%d" % (iorb+1)
                save_cubefile(atomlist, orbital,
                              filename=filename,
                              ppb=cube_opts.get("ppb", 2.0),
                              dbuff=cube_opts.get("dbuff", 5.0))
                print( "   orbital %d saved to %s_real.cube and %s_imag.cube" % (iorb+1, filename, filename) )

    if dyson_opts != {}:
        # load Dyson orbitals from molden file
        mf = MoldenFileInput(dyson_opts.get("molden_file", "dyson_orbitals.molden"))

        if not geometries_are_equal(atomlist, mf.getGeometry(), eps=5.0e-3):
            raise RuntimeError("Molecular geometry in molden file differs!")
        
        # Molecular geometries have to 
        dyson_orbitals = mf.getMOs()

        if dyson_opts.get("project", False) == True:
            # We project the Dyson orbitals onto the basis of bound orbitals, to make
            # sure the initial orbitals are orthogonal to all continuum orbitals.
            # Otherwise the transition dipole moments are origin-dependent.
            assert bound_opts != {}, "No bound orbitals found for projecting, add a \"bound\" section."

            print( "Dyson orbitals are projected onto bound orbitals to remove their continuum component." )
            dyson_orbitals = [project_orbital(bound_orbitals, d) for d in dyson_orbitals]
            
            if cube_opts.get("export_projected", False) == True:
                print( "exporting projected Dyson orbitals to cube files" )
                prefix = cube_opts.get("prefix", "/tmp/ms_")
                for iorb, orbital in enumerate(dyson_orbitals):
                    filename = prefix + "_dyson_proj_%d" % (iorb+1)
                    save_cubefile(atomlist, orbital,
                                  filename=filename,
                                  ppb=cube_opts.get("ppb", 2.0),
                                  dbuff=cube_opts.get("dbuff", 5.0))
                    print( "   projected Dyson orbital %d saved to %s_real.cube and %s_imag.cube" % (iorb+1, filename, filename) )
                

    # select orbitals for which PAD should be calculated
    initial_orbitals = []
    # Which orbitals should be used as initial states for photoionization,
    # bound orbitals or Dyson orbitals?
    if pad_opts.get("orbitals", "bound") == "dyson":
        # use Dyson orbitals as initial states
        assert dyson_opts != {}, "No Dyson orbitals found, please add a \"dyson\" section."
        for i in pad_opts.get("initial", range(1, len(dyson_orbitals)+1)):
            initial_orbitals.append( dyson_orbitals[i-1] )
    else:
        # use bound orbitals as initial states
        assert bound_opts != {}, "No bound orbitals found, please add a \"bound\" section."
        for i in pad_opts.get("initial", range(1, len(bound_orbitals)+1)):
            initial_orbitals.append( bound_orbitals[i-1] )
    
    pad = PhotoelectronAngularDistribution(muffin)

    pad.compute_pads(initial_orbitals, kinetic_energies,
                     pol=pad_opts.get("pol", 0),
                     pad_file=pad_opts.get("pad_file", "pad.dat"),
                     units=pad_opts.get("units", "eV-Mb"))

    # search for resonances in cross section
    if continuum_opts.get("resonances", False) == True:
        resonances = pad.find_resonances(sigma_thresh=continuum_opts.get("sigma_thresh", 10.0))

        if cube_opts.get("export_resonances", False) == True:
            print( "exporting resonances to cube files" )
            prefix = cube_opts.get("prefix", "/tmp/ms_")
            for iorb in sorted(resonances.keys()):
                for j,res in enumerate(resonances[iorb]):
                    filename = prefix + "_resonance_%d_to_%d" % (iorb+1, j+1)
                    save_cubefile(atomlist, orbital,
                                  filename=filename,
                                  ppb=cube_opts.get("ppb", 2.0),
                                  dbuff=cube_opts.get("dbuff", 5.0))
                    print( "   resonance %d to %d saved to %s_real.cube and %s_imag.cube" % (iorb+1, j+1, filename, filename) )
        
    

if __name__ == "__main__":
    import sys
    args = sys.argv[1:]
    if len(args) > 0:
        filename = args[0]
    else:
        print(" Please provide a JSON input file. " )
        exit(-1)
        
    with open(filename) as json_file:
        json_input = json.load(json_file)
        
    run_from_config(json_input["electroSka"])
