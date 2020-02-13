
To view the fitted transition charges and transition dipole moment run

  vmd -e show_transition_charges.tcl -args monomer.ff monomer.chromo

To view the unit cell of the tube run

  vmd -e show_unit_cell.tcl

To view the tube with the electric and transition dipoles drawn as blue and red arrows run

  vmd -e show_transition_charges.tcl -args tube_unit_cell.ff tube_unit_cell_ROTDIP.chromo


Parametrization of Frenckel exciton model
-----------------------------------------

* The monomer was optimized with AM1 followed by an optimization with PBEPBE/6-31+G*
* For the optimized geometry the lowest 4 excited states were calculated at the PBEPBE/6-311++G** level of theory.

    g09 < monomer.gjf > monomer.out
    formchk monomer.chk
    
* The transition density for the S0->S1 transition was calculated using the `Multiwfn` program (=> td_S0-S1.cube)

    Multiwfn monomer.fchk

* The electrostatic potential of the transition density was computed by solving the Poisson equation
  using the PSPFFT package.  (=> td_pot_S0-S1.cube)
  
   poisson.py td_S0-S1.cube td_pot_S0-S1.cube --nuclear_potential=0
   
* transition charges were fitted to the electrostatic potential using the CHELPG algorithm. (=> transition_charges_S0-S1.dat)

   chelpg.py td_pot_S0-S1.cube transition_charges_S0-S1.dat

