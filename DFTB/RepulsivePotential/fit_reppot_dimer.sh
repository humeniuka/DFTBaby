#!/bin/bash

usage="
 This script automatizes the fitting of repulsive potentials for an atom pair A-B
 from potential energy curves of the dimer.

 Usage:
    $(basename $0)  A  B    [charge]

 A,B are the atom symbols (e.g. h,c,br,etc) for atom 1 and atom 2.
 and optionally charge should be chosen such that the dimer is closed shell. For instance, O-Br
 is an open shell system, but [O-Br]^- is closed-shell.

 At the end the repulsive potential is displayed graphically and saved to '<atom 1>_<atom 2>.py'.

 The following steps are performed:
  1. create a sequence of dimer geometries with different bond lengths
  2. compute the electronic DFTB forces (without repulsion between ions and core electrons) for each geometry
  3. compute the full DFT forces with LC-PBE/6-311+G*
  4. fit the repulsive potential from the force differences

"

set -e # abort on first error

# check number of arguments
if [ $# -lt 2 ]
then
    echo "$usage"
    exit -1
fi

atom1=$1
atom2=$2
charge=${3:-0}

# set up minimal configuration for DFTB
cat <<EOF > dftbaby.cfg
[DFTBaby]
long_range_correction=1
dispersion_correction=0

verbose=0
EOF

# create scan geometries -> scan.xyz
create_dimer_trajectory.py $atom1 $atom2   0.5 5.0 100  scan.xyz

# compute DFTB forces along scan
dftb_forces_trajectory.py scan.xyz forces.dftb.xyz --charge=$charge | tee dftb.out

# compute LC-PBE forces along scan using Gaussian 09
gaussian_forces_trajectory.py scan.xyz forces.lc-pbe.xyz --charge=$charge | tee lc-pbe.out

# remove temporary folder created by Gaussian
rm -rf g09_temp

# extract energies along scan
grep "Structure" dftb.out   | awk '{print NR,$4}' > entot.dftb.dat
grep "Structure" lc-pbe.out | awk '{print NR,$4}' > entot.lc-pbe.dat

# fit repulsive potential from dimer curve
FitForces.py $atom1 $atom2  --reppot_dir="./" <<EOF
 $atom1-$atom2-dimer     scan.xyz   forces.dftb.xyz  forces.lc-pbe.xyz  weight=1
EOF

echo "The repulsive potential should have been saved to '${atom1}_${atom2}.py'."
