"""
solve Poisson equation for determining the electrostatic potential from the electronic density
and nuclear point charges. 

This module provides two Poisson solvers:

 1) The iterative solver that replaces the second order partial derivatives in
       __2
       \/  V(x) = -4 pi rho(x)

    with finite difference approximations. The resulting system of linear equations is solved by 
    Jacobi iterations.

 2) The Poisson equation is solved in Fourier space using the PSPFFT code [1], which
    can be downloaded from the Computer Physics Communications Program Library.

    The grid with the source data is written to a binary file. Then the external
    program ``poisson_pspfft.x`` is called, which writes the solution to a binary file,
    which is then read back in.


References:

[1] Budiardja,~R.;\ \ Cardall,~Ch.
    Parallel FFT-based Poisson solver for isolated three-dimensional systems.
    Computer Physics Communications 182 (2011) 22652275
"""
