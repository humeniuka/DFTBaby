### Metadynamics Input File ###
{
  ## definition of CVs ##
  # available types: bond, angle, torsion, coordination number, mulliken charge, custom
  "collective variables":  # list of CVs to be used in the metadynamics simulation
  [
    {  # CV1
      "name": "Bond between atom A and atom B",
      "type": "bond",
      "atoms": [0, 1],  # list of atom indices defining the bond
      "width": 0.1,  # width of the added gaussians for this CV (Angstrom)
      "plot":  # Parameters for plotting the potential (Angstrom)
      {
        "min": 0.,
        "max": 1.,
        "npoints": 100,
        "ticks": 0.2
      }
    },
    {  # CV2
      "name": "Angle between atoms A, B and C",
      "type": "angle",
      "atoms": [0, 1, 2],  # list of atom indices defining the angle
      "width": 10.,  # width of the added gaussians for this CV (degrees)
      "plot":  # Parameters for plotting the potential (degrees)
      {
        "min": 0.,
        "max": 180.,
        "npoints": 100,
        "ticks": 30
      }
    },
    {  # CV3
      "name": "Torsion angle defined by atoms A, B, C and D",
      "type": "torsion",
      "atoms": [0, 1, 2, 3],  # list of atom indices defining the torsion angle
      "width": 10.,  # width of the added gaussians for this CV (degrees)
      "plot":  # Parameters for plotting the potential (degrees)
      {
        "min": -180.,
        "max": 180.,
        "npoints": 100,
        "ticks": 45
      }
    },
    {  # CV4
      "name": "Coordination number of atom A wrt other atoms",
      "type": "cn",
      "atom": 0,  # atom index for that the cn should be calculated
      "reference": [2, 3, 4],  # list of atoms that contribute to the cn (default: all)
      "n": 6,  # first exponent
      "m": 12,  # second exponent
      "d": "auto",  # reference length, either float (Angstrom) or "auto"
      "width": 0.1,  # width of the added gaussians for this CV
      "plot":  # Parameters for plotting the potential
      {
        "min": 0.,
        "max": 1.,
        "npoints": 100,
        "ticks": 0.2
      }
    },
    {  # CV5
      "name": "Mulliken charge on one or multiple atoms",
      "type": "mullcharge",
      "atoms": [0, 1, 2, 3],  # single atom index or list of atom indices
      "width": 0.1,  # width of the added gaussians for this CV (e)
      "iface": "dftbaby",  # either dftbaby or turbomole (default: turbomole)
      "method": "numeric",  # either numeric or analytic (only for dftbaby)
      "plot":  # Parameters for plotting the potential (e)
      {
        "min": -1.,
        "max": 1.,
        "npoints": 100,
        "ticks": 0.25
      }
    },
    {  # CV6
      "name": "User-defined CV defined in cvcustom.py",
      "type": "custom",
      "prm":  # paramaters that will be passed to the functions defined in cvcustom.py
      {
          "key1": "value1",
          "key2": "value2"
      },
      "width": 0.1,  # width of the added gaussians for this CV
      "plot":  # Parameters for plotting the potential
      {
        "min": 0.,
        "max": 1.,
        "npoints": 100,
        "ticks": 0.2
      }
    }
  ],
  ## number of dynamics steps until the next gaussian is added ##
  "tau_g": 500,
  ## height of the added gaussians (hartree) ##
  "height": 0.0001
}
