Synopsis
========
DFTBaby is a software package for tight-binding DFT calculations on ground and excited states
of molecules and for non-adiabatic molecular dynamics simulations.

The latest version of DFTBaby and a short manual can be found on the following website:

  http://dftbaby.chemie.uni-wuerzburg.de/

Dependencies
============
 - python 2.7
 - numpy 1.8.2 or newer
 - scipy 0.14.0 or newer
 - f2py for compiling the Fortran extensions
 - matplotlib 2.0.0 or newer
 - mpmath 0.19 or newer
 - sympy 1.0 or newer
 - BLAS and LAPACK
 - libxc 3.0.0 (optional, only needed for atomic calculations of pseudoorbitals)


DFTBaby makes extensive use of numpy and scipy for performing linear algebra operations.
It's highly recommended to buy the Intel Math Kernel Library (MKL) and compile numpy/scipy
with MKL. MKL is highly optimized and automatically adds parallelism to many math operations.
Some extensions of DFTBaby are written in Fortran and are parallelised with OpenMP. To turn on parallel
execution on a shared memory machine (with 4 processors), the following environment variables
should be set (for bash shell):
  export OMP_NUM_THREADS=4
  export MKL_NUM_THREADS=4

*Warning*: With OMP_NUM_THREADS and MKL_NUM_THREADS set to 1, some parts of the code
may use up to 300% CPU, instead of the expected 100%. However setting those to
2 or 4 correctly results in 200% or 400% CPU usage, respectively.

Installation
============

  tar -xvf DFTBaby-###.tar.gz
  cd DFTBaby-###/

  python setup.py install --user

This will install the python libraries in $HOME/.local/lib/python2.7 and the executables
in $HOME/.local/bin. These paths should be added to the environment variables PYTHONPATH
and PATH, respectively.
In some cases, you might have to edit the setup.py script. For instance, if the optimal
BLAS or LAPACK library is not in a default location, this can be fixed by editing the
`extra_link_args` options. Also, if you wish to use the Intel compiler instead of the
GNU compiler, you have to uncomment the relevant section in setup.py.

Alternatively, DFTBaby can be installed by adding the folder DFTBaby-#### to the PYTHONPATH
environment variable. The subfolders DFTB/, DFTB/Analyse, DFTB/Modeling, DFTB/Formats, DFTB/Dynamics,
DFTB/Dynamics/Analyse and DFTB/Dynamics/Analyse/Viewer
should be added to the PATH environment variable. This setup has the advantage that you
can modify the code and use it without reinstalling it repeatedly.
The script 'dftbaby_env.sh' contains the necessary lines of code to set those environment
variables and can be sourced from the .bashrc profile.

If you are using environment modules you can use the file 'dftbaby.modulefile' as a template. 

In some cases you will have to compile the extensions for the architecture of your machine
with f2py:

  cd DFTB/extensions/
  make clean
  make
  cd -

If you whish to use the rudimentary DREIDING force field for QM/MM calculation, you have
to compile the code written in C as well:

  cd DFTB/ForceField/src/
  make clean
  make
  cd -

The (continuum) multiple scattering method allows to compute bound and continuum orbitals
of a muffin tin potential and to obtain photoelectron angular distributions. The critical
parts are written in Fortran and you will have to compile the code first:

  cd DFTB/MultipleScattering/src
  make clean
  make
  cd -


Programs
========
DFTBaby provides the module DFTB that can be imported from python. In addition
the following programs (among other scripts) can be called directly from the command line:

 DFTB2.py                 - only ground state calculation, MOs
 LR_TDDFTB.py             - ground and excited state calculations, gradients
 GeometryOptimization.py  - optimizes geometries on ground and excited states
 SurfaceHopping.py        - runs non-adiabatic molecular dynamics

All programs have a help function that shows all available command line options, e.g.:

  LR_TDDFTB.py --help


Features
========
+ Electronic Prametrization
   - Pseudo-orbitals can be generated from atomic DFT calculations
     using non-hybrid functionals from the libxc-library
   - Slater-Koster files for the hamiltonian, the overlap, and dipoles
     can be generated from the pseudo-orbitals
+ Repulsive Potentials
   - scaled repulsive potentials of Hotbit
+ DFTB
   - long-range correction
   - Grimme's dispersion correction (experimental)
+ TD-DFTB
   - excited states with and without long-range correction
   - analytical excited state gradients
+ QM/MM
   - QM part treated with TD-DFTB, MM part with UFF of Gaussian 09
   - QM/MM calculations of molecular crystals (with periodic boundary conditions) using
     rudimentary implementation of the DREIDING force field
+ Analysis
   - cube files for molecular orbitals and transition and difference densities of excited states
+ Non-adiabatic Dynamics
   - surface hopping
   - electronic coefficients are integrated in the locally diabatic basis

Graphical Analysis Tool
=======================
The results of a TD-DFTB calculation can be visualized graphically. To start the graphical user interface
at the end of a calculation it is enough to add the option --graphical=1 when calling LR_TDDFTB.py.
The graphical interface requires Enthought's 'Mayavi' which can be installed via

  pip install mayavi


Limitations
============
The following known limitations should be kept in mind when using (lc)-TD-DFTB:

- DFTB uses a minimal basis set of valence orbitals. Being numerically exact
  atomic orbitals obtained from a DFT calculations, the pseudo-orbitals are much better
  than for instance, STO-3G basis functions. However, Rydberg states, which are frequent
  in organic conjugated molecules, cannot be described in this way.

- The Coulomb interaction between density fluctuations around the neutral reference density
  is calculated using the monompole approximation. In reality, the charge fluctuations around an atom
  are not spherically symmetric. In my experience the monopole approximation leads to occupied Kohn-Sham
  orbitals that are too high in energy. Excitations from these orbitals to unoccupied ones
  produce spurious low-lying states, that are dark.
  Adding the interaction between Mulliken monopoles and dipoles (by setting the option --mulliken_dipoles=1)
  should in principle fix this problem, however the Mulliken dipoles are too large and can
  lead to SCF convergence problems.

- For large molecules, the program will run out of memory during the calculation of the gradients.
  This is a problem of the implementation, since the entire matrix of the gradients of S and H0 is
  constructed and kept in memory.

- Repulsive potentials (from Hotbit) are only available for the elements H,C,N and O. The electronic structure
  does not depend on the repulsive potentials, but the program will abort if it cannot find the repulsive
  potentials for an atom combination present in the input geometry. Therefore, files with dummy potentials for
  other elements such as S, Ag etc. exist in the directory 'reppot_tables', but they are not suitable for
  optimization or dynamics simulations.
  Caution is advised when copying repulsive potentials from other DFTB programs, that were obtained for a
  different electronic parametrization.
  
  

Non-adiabatic Dynamics with DFTB
================================
Performing non-adiabatic dynamics with DFTBaby is quite simple. You have to prepare a folder for the trajectory
with the following content:

  dynamics.in        an xyz-file with the initial geometry in bohr followed by the initial velocities (also in a.u.)
  dftbaby.cfg        a configuration file with options for the DFTB program. The same options can be specified as for
                     the LR_TDDFTB.py program, in addition to options controlling the dynamics simulation such
		     as the initial state, number of time-steps, etc. see DFTB/dftbaby.cfg for a template.

Calling the program 'SurfaceHopping.py' inside the folder will start the dynamics. For each time step
the geometries, excitations energies and currrent states will be stored in the files
dynamics.xyz, energy_#.dat and state.dat. 

An example can be found in

   examples/NAMD/

