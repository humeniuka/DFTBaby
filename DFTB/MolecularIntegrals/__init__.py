"""
This module contains routines for

  - molecular integrals between Gaussian type orbitals (overlap, kinetic, dipole, angular momentum, electron repulsion integrals)
  - numerical integration on a multicenter grid (Becke's scheme)
  - solving Poisson's equation numerically on a multicenter grid
  - reading formatted Gaussian checkpoint files
  - computing associated Legendre polynomials

DFTB does not need to calculate molecular integrals at runtime. These routines are only used for testing or tabulating integrals.
"""
