
 Orbital Localization (Pipek-Mezey)
 ----------------------------------

This example illustrates how to use localized orbitals for defining
initial electronic wavefunctions with well-defined character in MD simulations.

First inspect the localized orbitals by running

   DFTB2.py ethene_dimer.xyz --localize_orbitals=PM

and opening the newly created file 'localized_orbitals.molden'.

The initial electronic wavefunction defined in 'dftbaby.cfg' corresponds
to a local excitation on one ethene molecule. The dynamics is started with

  SurfaceHopping.py

and should run until a conical intersection to the ground state is encountered.


