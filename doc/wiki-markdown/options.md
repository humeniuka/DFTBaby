Command Line Options
--------------------
The **DFTBaby** programs understand the following options, which can be either appended to
the command line and can be put into a file called `dftbaby.cfg` in the current directory:

<pre>
Usage: LR_TDDFTB.py .xyz-file
   --help option will give more information


Options:
  -h, --help            show this help message and exit
  --nr_unpaired_electrons=NR_UNPAIRED_ELECTRONS
                         number of unpaired electrons for open shell singlet
                        states [default: 0]
  --use_symmetry=USE_SYMMETRY
                         (0 = off, >0 = on) Brings the molecule into standard
                        orientation and tries to classify excited states by
                        the irreducible representations of the symmetry group
                        detected. The assignment only works for non-degenerate
                        irreps. Not all point groups implemented. Note that
                        symmetry adapted orbitals are not used, therefore this
                        option does not speed up the calculation. The symmetry
                        is only determined for the N=use_symmetry lowest
                        excited states. [default: 0]

  Parametrization:
    --parameter_set=PARAMETER_SET
                        "homegrown" - my own parametrization, "hotbit"  -
                        hotbit parameters, "mio" - mio-0-1 DFTB+ parameters
                        set,  [default: homegrown]
    --distance_cutoff=DISTANCE_CUTOFF
                         orbitals whose centers are further away than this
                        cutoff do not interact (in bohr) [default: 30.0]
    --fluctuation_functions=FLUCTUATION_FUNCTIONS
                         Functional form of the spherical charge fluctuations
                        dn(r) on which the gamma-matrix is based, can be
                        'Gaussian' or 'Slater' [default: Gaussian]

  Charges:
    --point_charges_xyz=POINT_CHARGES_XYZ
                         path to file with point charges which interact
                        electrostatically with the molecular partial charges
                        [default: none]
    --initial_charge_guess=INITIAL_CHARGE_GUESS
                         path to file with initial partial charges, on each
                        line the partial charge for one atom should be given
                        in the same order as in the xyz file [default: none]
    --save_converged_charges=SAVE_CONVERGED_CHARGES
                         path to a file where the partial charges of the
                        converged SCF calculation are stored, on each line the
                        partial charge for one atom is written in the same
                        order as in the xyz file. A subsequent calculation can
                        use this as an initial guess via the
                        'initial_charge_guess' option. [default: none]

  Output:
    --verbose=VERBOSE    controls amount of output,0 - silent, 1 - in each
                        iteration write energies and Mulliken charges, >1 -
                        print matrix elements and MO coefficients,  [default:
                        1]
    --spectrum_file=SPECTRUM_FILE
                         write table with excitation energies and oscillator
                        strengths to this file. If the states were classified
                        by symmetry a seperate file is written for each
                        irrep.,  [default: absorption_spectrum.dat]

  Long-Range-Correction:
    --long_range_correction=LONG_RANGE_CORRECTION
                        0 - disabled, 1 - add exact exchange for large
                        distances,  [default: 1]
    --long_range_radius=LONG_RANGE_RADIUS
                        distance R_lr in bohr where the long range correction
                        is switched on,  [default: 3.03]
    --long_range_switching=LONG_RANGE_SWITCHING
                         "erf" or "erfgau", determines how the exact HF
                        exchange is turned on as a function of the distance,
                        "erf" => erf(R_AB/R_lr), "erfgau" => erf(R_AB/R_lr) -
                        2/(sqrt(pi)*R_lr)*R_AB*Exp(-1/3 (R_AB/R_lr)^2)
                        [default: erf]
    --lc_implementation=LC_IMPLEMENTATION
                         The long-range correction can be implemented in
                        slightly different ways. In order to reproduce the
                        calculations from JCP 143 134120 (2015) set this
                        option to "old" [default: none]

  Experimental (not tested properly):
    --long_range_T=LONG_RANGE_T
                        parameter T in Fermi function, see option
                        long_range_switching,  [default: 1.0]
    --tune_range_radius=TUNE_RANGE_RADIUS
                        optimize the long range radius such that the
                        ionization potential equals the negative of the HOMO
                        energy,  [default: 0]
    --save_tuning_curve=SAVE_TUNING_CURVE
                         Save a chart of J(R_lr) plotted agains R_lr to this
                        file [default: none]
    --mulliken_dipoles=MULLIKEN_DIPOLES
                         1 - the electrostatic interaction between monopoles-
                        dipoles and dipoles-dipoles is included in the coulomb
                        part of the hamiltonian, 0 - only monopole-monopole
                        interaction [default: 0]
    --dispersion_correction=DISPERSION_CORRECTION
                         1 - the Grimme dispersion correction is added to the
                        repulsive potential. [default: none]
    --scratch_dir=SCRATCH_DIR
                         Large arrays can be mapped to scratch files to
                        economize memory usage. These files will be created
                        temporarily in the directory provided in
                        'scratch_dir'. Functions using this mapping will run
                        much slower,  [default: none]

  QM/MM:
    --qmmm_partitioning=QMMM_PARTITIONING
                         The system is partitioned into an inner region that
                        is treated with DFTB (QM-part) and an outer region
                        treated with the UFF force field (MM-part). This
                        option should contain the indeces of the atoms in the
                        xyz-file that should belong to the QM part. The
                        counting should start at 1. Indeces can be a list,
                        e.g. "[1,2,3]", or a range, e.g. "range(1,4)", or a
                        combination, e.g. "[1,2,3]+range(4,10)". Since the MM
                        atom types are assigned automatically you need to make
                        sure that the geometry in the xyz-file is optimized
                        and reasonable for a force field. The UFF calculations
                        require a working installation of Gaussian 09
                        [default: none]
    --qmmm_embedding=QMMM_EMBEDDING
                         'mechanical' or 'electrostatic'. The 'mechanical'
                        embedding neglects all electrostatic interactions
                        between MM atoms. In the case of 'electrostatic'
                        embedding, the MM charges are not included into the
                        DFTB-hamiltonian and therefore do not polarize the
                        orbitals, however, the electrostatic interaction
                        between MM charges in the inner and outer region is
                        included. [default: electrostatic]
    --periodic_force_field=PERIODIC_FORCE_FIELD
                         If periodic boundary conditions are require, a
                        rudimentary implementation of the DREIDING force field
                        is used instead of UFF. This parameter should contain
                        the path to a file with the force field definition.
                        The file should specify atom position, atom types and
                        lattice vectors. Examples are given in
                        DFTB/ForceField/examples/. [default: none]

  Cavity:
    --cavity_radius=CAVITY_RADIUS
                         Sets the cavity radius r0 in bohr. A spherical
                        confining cavity with radius r0 holds the geometry
                        together. The walls push atoms back towards the
                        center. The potential is Vcav(r) = 1/2 * k * (r-r0)^2
                        for r > r0. If r0==0, the confinement is turned off.
                        [default: 0.0]
    --cavity_force_constant=CAVITY_FORCE_CONSTANT
                         Sets the force constant k [Hartree/bohr^2]. [default:
                        0.01]

  CPKS:
    --cpks_solver=CPKS_SOLVER
                         Determines how the Coupled-perturbed Kohn-Sham (CPKS)
                        equations should be solved - iteratively with a Krylov
                        subspace method ('iterative') or directly ('direct').
                        The direct solution requires an extremely large amount
                        of memory, whereas the iterative method takes a long
                        time. The CPKS equations are only solved, if gradients
                        of the MO coefficients are required. Gradients of the
                        ground and excited state are computed using Furche's
                        auxiliary function method, so the CPKS equations are
                        not needed in this case.,  [default: direct]

  SCF-Convergence:
    --maxiter=MAXITER    stop SCF calculation after maxiter iterations
                        [default: 1000]
    --scf_conv=SCF_CONV
                         convergence threshold for relative change in SCF-
                        calculation [default: 1e-10]
    --start_from_previous=START_FROM_PREVIOUS
                         use density matrix and partial charges from a
                        previous calculation as initial guess if available (1)
                        or start from reference density and zero partial
                        charges (0),  [default: 1]
    --mixing_threshold=MIXING_THRESHOLD
                         if the relative change drops below this value density
                        mixing is used [default: 0.001]
    --density_mixer=DENSITY_MIXER
                         instance of Mixer object, that determines how the
                        next densitymatrix is obtained from previous guesses,
                        if None a linear mixer is used,  [default:
                        <DFTB.ConvAcceleration.DIIS_80 instance at
                        0x2b536947d368>]
    --fock_interpolator=FOCK_INTERPOLATOR
                         the KS hamiltonian is interpolated based on previous
                        iteration steps [default: none]
    --level_shift=LEVEL_SHIFT
                         shift virtual orbitals up in energy, this shift
                        parameter is gradually reduced to zero as the density
                        matrix converges [default: 0.1]
    --HOMO_LUMO_tol=HOMO_LUMO_TOL
                         level shifting is turned on, as soon as the HOMO-LUMO
                        gap drops below this value [default: 0.05]
    --linear_mixing_coefficient=LINEAR_MIXING_COEFFICIENT
                         is no density mixer object is used
                        (density_mixer=None) the next guess for the density
                        matrix is constructed as P_next = a*P + (1-a)*P_last
                        [default: 0.33]

  Linear-Response TD-DFTB:
    --nr_active_occ=NR_ACTIVE_OCC
                        integer, number of active occupied orbitals
                        if set to None, all occupied orbitals are selected,
                        [default: none]
    --nr_active_virt=NR_ACTIVE_VIRT
                        integer, number of active virtual orbitals
                        if set to None, all virtual orbitals are included,
                        [default: none]
    --select_lm=SELECT_LM
                        Molecular orbitals are included in the active space if
                        they contain atomic orbitals with angular quantum
                        numbers (l,m). This options allows to select only the
                        pz-orbitals in a planar conjugated molecule to mimic a
                        PPP calculation, e.g. --select_lm='[(1,0)]'. If
                        symmetry is enabled, the molecule will be reoriented
                        and the pz-orbitals might actually be the px or py
                        orbitals depending on the point group.,  [default:
                        none]
    --oszis=OSZIS       method by which oscillator strengths are calculated,
                        "Mulliken" -> from Mulliken transition charges,
                        "Dipoles" -> from transition dipole moments (not
                        available with 'mio' or 'hotbit' parameters),
                        [default: Dipoles]
    --response_method=RESPONSE_METHOD
                        "Casida" (for RPA) or "Tamm-Dancoff",  [default:
                        Casida]
    --multiplicity=MULTIPLICITY
                        "S" for Singlets or "T" for Triplets,  [default: S]

  Davidson-like Diagonalization:
    --nstates=NSTATES    Solve iteratively for the lowest 'nstates'
                        eigenvalues of the non-Hermitian TD-DFTB equations. If
                        this options if not set, the full TD-DFTB matrix is
                        constructed in memory and diagonalized. [default:
                        none]
    --diag_ifact=DIAG_IFACT
                         ifact*nstates singly excited states are created as
                        initial guesses. This factor should be increased if
                        the resulting states do not have the desired symmetry.
                        [default: 1]
    --diag_conv=DIAG_CONV
                         convergence threshold for eigenvalues [default:
                        1e-05]
    --diag_maxiter=DIAG_MAXITER
                         maximum number of iterations [default: 50]
    --diag_check=DIAG_CHECK
                         compare iterative solution with the full
                        diagonalization (only for debugging purposes)
                        [default: 0]
    --diag_L2threshold=DIAG_L2THRESHOLD
                         solution vectors whose Lambda2 value is below this
                        threshold are removed. [default: 0.0]
    --diag_selected_irreps=DIAG_SELECTED_IRREPS
                         If symmetry is enabled, expansion vectors belonging
                        to irreps listed in diag_selected_irreps are included
                        preferentially. Note the use of quotes and double
                        quotes, e.g. --diag_selected_irreps="['B1U','B2U']".
                        [default: none]

  Long-Range Charge-Transfer:
    --ct_correction=CT_CORRECTION
                         0 -> no correction, 1 -> shifts low-lying spurious CT
                        states to higher energies as described in J.Chem.Phys.
                        124, 213102 (2006),  [default: 0]

  Charge-Transfer-Analysis:
    --particle_hole_states=PARTICLE_HOLE_STATES
                         list of excited states (starting from 1) for which
                        the positions of the particle and the hole should be
                        determined. [default: []]
    --particle_hole_dir=PARTICLE_HOLE_DIR
                         The charge density of the particle and the hole are
                        written to this file. The density is coarse grained to
                        atomic resolution, so for each atom the portion of the
                        delocalized particle or hole sitting on that atom is
                        save. The atoms are listed in the same order as in the
                        xyz file.,  [default: none]

  Graphical-Analysis:
    --graphical=GRAPHICAL
                         If graphical > 0, a graphical user interface is
                        opened for analysing the results of the calculations.
                        Requires mayavi and a running X-server.,  [default: 0]

  Cube-Files:
    --cubedir=CUBEDIR    directory where cube files of transition or
                        difference densities are stored [default: none]
    --cube_states=CUBE_STATES
                         indeces of excited states for which transition or
                        difference densities should be computed, the first
                        excited state has index 1! [default: [1,2,3]]
    --density_type=DENSITY_TYPE
                         specifies whether transition densities ('transition')
                        or density differences ('difference') between the
                        excited states and the ground state are calculated
                        [default: transition]
    --points_per_bohr=POINTS_PER_BOHR
                         resolution of grid for cube files [default: 2.0]
    --cube_threshold=CUBE_THRESHOLD
                         atomic orbitals with coefficients below
                        threshold*(largest coefficient) are neglected when
                        evaluating MO amplitudes on a cube grid [default:
                        0.001]
    --save_transition_density=SAVE_TRANSITION_DENSITY
                         if 1, the orthogonalized transition density P =
                        S^(1/2).P.S^(1/2)^T in the AO basis is also saved for
                        each state in cube_states as a matrix in a plain text
                        file with suffix '.mat',  [default: 0]

  Molden-File:
    --molden_file=MOLDEN_FILE
                         Save molecular Kohn-Sham orbitals to this file in the
                        Molden format [default: none]
    --molden_nr_occ=MOLDEN_NR_OCC
                         only these highest occupied orbitals will be exported
                        to the Molden file [default: 100]
    --molden_nr_virt=MOLDEN_NR_VIRT
                         only these lowest unoccupied orbitals will be
                        exported to the Molden file,  [default: 100]

  Gradients:
    --gradient_state=GRADIENT_STATE
                         Compute the gradient on this state analytically, 0
                        means ground state, 1 means 1st excited states etc.
                        [default: 2]
    --gradient_file=GRADIENT_FILE
                         save the total gradient (Vrep + E0 + ExcI) to this
                        xyz-file. If this option is not set, the gradient will
                        not be calculated at all! [default: none]
    --gradient_check=GRADIENT_CHECK
                         compute gradient both numerically and analytically
                        and compare (0 -> off, 1 -> perform check). You should
                        increase the SCF-convergence threshold to 1.0e-14 and
                        disable the DIIS mixer.,  [default: 0]
</pre>

The `SurfaceHopping.py` program understands the following options:

<pre>
Usage: SurfaceHopping.py
  run surface hopping dynamics.

  Required Input Files:
    - `dynamics.in` with the initial coordinates and velocities (in bohr)
    - `dftbaby.cfg` with the options controlling the TD-DFTB calculation and dynamics
  Output Files:
    - `dynamics.xyz`: geometries for each time step (in Angstrom)
    - `state.dat`: current electronic state
    - `energy_I.dat`: total energy of electronic state I in Hartree vs. time in fs.
        The ground state energy at time t=0 is subtracted from all later energies.
    - `nonadiabaticIJ.dat`: scalar non-adiabatic coupling between states I and J.


Options:
  -h, --help            show this help message and exit

  Molecular Dynamics:
    --charge=CHARGE      total charge of the molecule, the cation or anion has
                        to be a closed shell species for TD-DFTB to work
                        properly. [default: 0]
    --initial_state=INITIAL_STATE
                         initial electronic state of the trajectory, 0 for
                        ground state. [default: 0]
    --nstates=NSTATES    number of excited states. Only the lowest states are
                        calculated with TD-DFTB. For dynamics in the ground
                        state `nstates` should be set to 0 to avoid the
                        expensive calculation of excitation energies and non-
                        adiabatic couplings. [default: 2]
    --nstep=NSTEP        number of nuclear time steps. [default: 1000]
    --nuclear_step=NUCLEAR_STEP
                         length of nuclear time step for integration of
                        Newton's equations (in fs). [default: 0.1]
    --dyn_mode=DYN_MODE
                         'T' for constant temperature, 'E' for constant
                        energy. To equilibrate a trajectory on the ground
                        state use 'T', non-adiabatic dynamics on excited
                        states should be run at constant energy. [default:
                        "E"]
    --temp=TEMP          temperature in Kelvin, only needed if the dynamics is
                        run at constant temperature [default: 300.0]
    --scalar_coupling_threshold=SCALAR_COUPLING_THRESHOLD
                         Excitation coefficients that are smaller than this
                        threshold are neglected when calculating scalar
                        couplings from the overlap between electronic states
                        at successive time steps. For large molecules this
                        value should be reduced. [default: 0.01]
    --switch_to_groundstate=SWITCH_TO_GROUNDSTATE
                         If set to 1, a hop to the ground state is forced if
                        the S0-S1 energy gap drops below 0.1 eV. In TD-DFT(B)
                        conical intersections to the ground state are not
                        described correctly. If a point of point of degeneracy
                        between S0 and S1 is reached, the TD-DFT(B)
                        calculation usually breaks down.If something goes
                        wrong in the excited state calculation, the trajectory
                        continues with the ground state gradients until the
                        excited state calculation starts working again.,
                        [default: 1]
    --artificial_energy_conservation=ARTIFICIAL_ENERGY_CONSERVATION
                         Energy conservation can be enforced artificially by
                        rescaling the velocities after each time step. This
                        avoids the drift of the total energy observed for long
                        dynamics simulations. This option only takes effect if
                        dyn_mode=="E" [default: 0]
    --time_series=TIME_SERIES
                         To better analyze trajectories it is helpful to write
                        out time series of additional quantities along the
                        trajectory. You can specify a list of the desired
                        quantities (so far only --time_series="['particle-hole
                        charges']" is available).,  [default: []]

</pre>