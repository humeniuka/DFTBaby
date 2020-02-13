#!/bin/bash
#
# optimize all geometries in VALIDATION/STRUCTURES using
# LC-DFTB
#
set -e

module load dftbaby

for xyz in STRUCTURES/*.xyz
do
    name=$(basename $xyz .xyz)
    echo "Optimizing $name"
    # prepare input file
    Nat=$(head -n 1 $xyz)
    cp $xyz DFTB/geom.xyz
    #
    optimize.py DFTB/geom.xyz 0
    # extract optimized geometry
    tail -n $(expr $Nat + 2) DFTB/geom_opt.xyz > DFTB/${name}.xyz

    # delete temporary files
    rm -f DFTB/geom.xyz DFTB/geom_opt.xyz

    echo "optimized geometry written to DFTB/${name}.xyz"
done
