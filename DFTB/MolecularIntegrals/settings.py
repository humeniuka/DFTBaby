# -*- coding: utf-8 -*-
"""
settings controlling the resolution of the grid in Becke's 
multicenter integration

You can increase or reduce the resolution of the grid by importing
this module and changing the parameters:

   from DFTB.MolecularIntegrals import settings
   settings.radial_grid_factor = 10
   settings.lebedev_order = 23

"""

# The mesh size is controlled by the following parameters.
radial_grid_factor = 3      # controls size of radial grid 
lebedev_order = 23          # controls size of angular grid

