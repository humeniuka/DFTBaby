#!/bin/bash
#
# optimize all geometries in VALIDATION/STRUCTURES using the
# reference method (LC-PBE/6-311+G*)
#
set -e

module load chem/g09

for xyz in STRUCTURES/*.xyz
do
    name=$(basename $xyz .xyz)
    echo "Optimizing $name"
    # prepare input file
    Nat=$(head -n 1 $xyz)
    cat opt_template.gjf > REFERENCE/${name}.gjf
    tail -n $Nat $xyz >> REFERENCE/${name}.gjf
    echo -e "\n\n\n\n" >> REFERENCE/${name}.gjf
    #
    g09 < REFERENCE/${name}.gjf > REFERENCE/${name}.out
    # extract optimized geometry
    newzmat -ichk -oxyz opt opt
    # extract total energy
    formchk opt.chk &> /dev/null
    enTot=$(grep "Total Energy" opt.fchk | awk {'print $4'})
    # save energy and geometry to xyz-file
    echo -e "${Nat}\n energy=${enTot}" > REFERENCE/${name}.xyz
    cat opt.xyz >> REFERENCE/${name}.xyz

    # delete temporary files
    rm -f opt.chk opt.fchk opt.xyz
    rm -f REFERENCE/${name}.out
    rm -f REFERENCE/${name}.gjf

    echo "optimized geometry written to REFERENCE/${name}.xyz"
done
