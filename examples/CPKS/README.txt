
Coupled-Perturbed Kohn-Sham (CPKS) equations
--------------------------------------------

This example explains how to get the gradients of Mulliken charges.
The script `charge_gradients.py` illustrates the sequence of function calls
that have to be made to prepare such a gradient calculation. Running

    ./charge_gradients.py water.xyz

computes the ground state gradient, then solves the coupled-perturbed Kohn-Sham equations
to get the gradients of the MO coefficients and finally prints the gradients of the
Mulliken charges.

