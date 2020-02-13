Scripts
-------

* `scan_internal_dgf.py`
<pre>
Usage: scan_internal_dgf.py <xyz file> <name of varying degree of freedom>
  The scan variable can be given either by its name (e.g. 'r(c0-c1)') or by its index in the z-matrix (e.g. 0).
  Angles are given in degrees and distances in Angstrom.
  
Options:
  -h, --help            show this help message and exit
  --scan_range=SCAN_RANGE
                        Tuple with (icmin, icmax, Npts) giving the range over
                        which the internal coordinate is varied. [default:
                        (1.0, 2.0, 10)]
  --out_xyz=OUT_XYZ     save the geometries of the scan to this xyz-file
                        [default: scan.xyz]
  --out_var=OUT_VAR     save values of internal coordinate along scan to a
                        table [default: scan_var.dat]

Example:
  $ scan_internal_dgf.py ethene.xyz 0

  will produce the following output:

Internal Coordinates
====================
6
1   1.8771015
1   3.2392549   30.364 deg
6   3.8296891   64.981 deg   0.000 deg
1   1.8771015   145.382 deg   0.000 deg
1   3.2392549   30.364 deg   180.000 deg
                Z-Matrix:
                =========
Index  Coordinate                     Value
---------------------------------------------------
 0   r(c0-h1)                         0.9933194 Ang
 1   r(h1-h2)                         1.7141400 Ang
 2   angle(h1-c0-h2)                 30.3636828 deg
 3   r(h2-c3)                         2.0265844 Ang
 4   angle(h2-h1-c3)                 64.9813302 deg
 5   dihedral(c0-h1-h2-c3)            0.0000000 deg
 6   r(c3-h4)                         0.9933194 Ang
 7   angle(c3-h2-h4)                145.3823526 deg
 8   dihedral(h1-h2-c3-h4)            0.0000000 deg
 9   r(h4-h5)                         1.7141400 Ang
10   angle(h4-c3-h5)                 30.3636828 deg
11   dihedral(h2-c3-h4-h5)          180.0000000 deg
Scan of the variable r(c0-h1) from 1.0 to 2.0 in 10 points
Scan geometries saved to scan.xyz
Values of scan variables saved to scan_var.dat
</pre>

* `scan.py`
<pre>
Usage: scan.py <xyz-file> <nr. states> <dat-file with energies>
  scans the excitation energies along the geometries in the xyz-file
  and writes a column with the energy (in eV) for each state
  (including ground state) to the dat-file
  type --help to see all options
  to reduce the amount of output add the option --verbose=0
</pre>

* `show_molcoords.py`
<pre>
Usage: show_molcoords.py <list of coordinate specifications and xyz-files>

shows bond lengths, angles and dihedrals for selected combinations
of atoms in an xyz-file. 

Coordinates are specified as quoted strings containing the indeces
of
  2 atoms A and B     for the distance A-B (e.g. "1 2")
  3 atoms A,B and C   for the angle B-A-C (e.g. "4 10 3")
  4 atoms A,B,C and D for the dihedral angle between the planes
                      spanned by A->B,B->C and B->C,C->D.

Atom indeces start at 1. Distances are in Angstrom and angles in degrees.

The selected internal coordinates are printed for each geometry in the 
file. If there are several xyz-files the coordinates are averaged over all of them. 
File names and coordinate specifications can be interspersed on the command line.

Examples:
   "1 2" water.xyz "1 2 3"    
          -> shows distance between atoms 1 and 2 and the angle 1-2-3.
   traj1.xyz traj2.xyz traj3.xyz "3 4 5"
          -> shows the bond angle 3-4-5 for each time step averaged over all 3 trajectories.
</pre>

* `interpolate_linearly.py`
<pre>
Usage: interpolate_linearly.py <initial geometry .xyz> <final geometry .xyz> N <path .xyz>
  The cartesian geometries are interpolated linearly between the initial
  and the final structures using N-2 intermediate points.
  The interpolated geometries are written to <path .xyz>
</pre>

* `ff_optimize.py`
<pre>
Usage: ff_optimize.py <.ff force field file>
  optimize geometry using periodic force field
  Force field files can be found in the folder DFTB/ForceField/examples/
  type --help to see all options.


Options:
  -h, --help            show this help message and exit
  --nr_unit_cells=NR_UNIT_CELLS
                        In the output the N neighbouring unit cells are added
                        to the geometry to give an impression of a larger
                        block of the crystal. [default: 0]
</pre>

* `absorption_spectrum.py`
<pre>
Usage: python absorption_spectrum.py <list of files with tabulated absorption spectrum>
  use --help to see all options

Options:
  -h, --help            show this help message and exit
  --save_png=SAVE_PNG   save plot to this file [default: ]
  --broadening=BROADENING
                        convolve stick spectrum with a Lorentzian of
                        FWHM=broadening (in Hartree) [default: 0.005]
  --in_units=IN_UNITS   input units, 'Hartree', 'eV', 'nm', 'cm-1' [default:
                        Hartree]
  --out_units=OUT_UNITS
                        output units, 'Hartree', 'eV', 'nm', 'cm-1' [default:
                        eV]
  --plot_range=PLOT_RANGE
                        plot range in output units, e.g. (0.0, 6.0) to plot
                        from 0 to 6 eV [default: ]
  --plot_labels=PLOT_LABELS
                        show file names as plot labels [default: 1]
  --plot_sticks=PLOT_STICKS
                        plot stick spectrum, too [default: 1]
  --vertical=VERTICAL   use separate subplot for each spectrum [default: 0]
  --invert_xaxis=INVERT_XAXIS
                        Invert energy-axis so that it points from right to
                        left [default: 0
  --save_converted=SAVE_CONVERTED
                        convert the spectrum to out_units and save it to this
                        file [default: ]
  --save_convolved=SAVE_CONVOLVED
                        convert the spectrum to out_units, convole it and save
                        it to this file [default: ]
  --transient_weights=TRANSIENT_WEIGHTS
                        Produce a transient absorption spectrum. Each spectrum
                        is weighted by the time-dependent concentration of the
                        intermediate to which it belongs. The concentrations
                        are read from the data file <transient_weight>, the
                        first column in the file is the reaction time followed
                        by as many columns as there are spectra.
  --transient_difference=TRANSIENT_DIFFERENCE
                        For transient spectra: subtract spectrum at time t=0
                        from all subsequent spectra.
  --labels=LABELS       List of labels to be used for each absorption spectrum
                        [default: ]
  --log_plot=LOG_PLOT   To sho states with very low oscillator strengths plot
                        the y-axis in logarithmic scale [default: 0]
  --normalize           Normalize spectrum so that maximum has intensity of 1
                        [default: none]
</pre>

