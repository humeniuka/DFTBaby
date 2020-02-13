
Non-adiabatic molecular dynamics on fluorene starting in the S1 excited state.
Calling

   SurfaceHopping.py

inside the folder will run the example. The configuration file `dftbaby.cfg` contains additional
key-value pairs related to controlling the dynamics simulation, such as the initial state and the size and number
of nuclear time steps. Output is written to several files, the first column usually contains the time in fs.

  dynamics.xyz         -  geometry along the trajectory
  energy_#.dat         -  adiabatic energies for each state
  state.dat            -  current electronic state
  coeff_#.dat          -  quantum populations for each state


