
Free energy profile for fluorene dimer from Meta-Dynamics simulation
====================================================================

The following input files have to be present:

   -  dynamics.in      with initial geometry and velocities
   -  dftbaby.cfg      configuration file for DFTB
   -  meta-config.py   definition of collective variables

The collective variables are defined in `meta-config.py`:
  - the opening angle of the bridge between the two monomers, angle(C5-C26-C29)
  - the dihedral angle between the monomers,  dihedral(C6-C5-C29-C30)

To request a metadynamics simulation in `dftbaby.cfg` one has to set dyn_mode='M'.

The simulation is run by starting

  SurfaceHopping.py

inside this folder. The reconstructed free energy surface can be viewed with

  reconstruct.py

The image 'vg_2500-g.png' shows a plot of the reconstructed free energy surfaces after
accumulating 2500 Gaussians. 

