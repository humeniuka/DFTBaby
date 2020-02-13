### Metadynamics Input File ###
{
  ## definition of CVs ##
  # available types: bond, angle, torsion, coordination number, mulliken charge, custom
  "collective variables":  # list of CVs to be used in the metadynamics simulation
  [
    {  # CV2
      "name": "Angle between atoms A, B and C",
      "type": "angle",
      "atoms": [4, 25, 28],  # list of atom indices defining the angle
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
      "atoms": [5, 4, 28, 29],  # list of atom indices defining the torsion angle
      "width": 10.,  # width of the added gaussians for this CV (degrees)
      "plot":  # Parameters for plotting the potential (degrees)
      {
        "min": -180.,
        "max": 180.,
        "npoints": 100,
        "ticks": 45
      }
    },
  ],
  ## number of dynamics steps until the next gaussian is added ##
  "tau_g": 100,
  ## height of the added gaussians (hartree) ##
  "height": 0.0001
}
