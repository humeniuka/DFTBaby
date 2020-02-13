#!/bin/bash
#
# compute reference and electronic DFTB energies and forces for all fit paths
#
set -e

if [ $# -eq 0 ]
then
    echo "For which method should the energies and forces along the fit paths be calculated?"
    select target in "reference" "dftb";
    do
	echo "Selected target: $target"
	break
    done
else
    target=$1
fi

function calc_reference() {
    module load chem/g09

    # reference total energies and forces
    echo "Reference"
    for xyz in FIT_PATHS/*.xyz
    do
	name=$(basename $xyz .xyz)
	echo $name
	gaussian_forces_trajectory.py --method="LC-PBEPBE" --basis_set="6-311+G*" $xyz REFERENCE/${name}.energies.dat REFERENCE/${name}.forces.xyz
    done
    # remove temporary folder
    rm -rf g09_temp
}

function calc_dftb() {
    echo "electronic DFTB"
    # electronic DFTB energies and forces
    for xyz in FIT_PATHS/*.xyz
    do
	name=$(basename $xyz .xyz)
	echo $name
	dftb_forces_trajectory.py $xyz DFTB/${name}.energies.dat DFTB/${name}.forces.xyz
    done
}

case "$target" in
    reference)
	calc_reference
	;;
    dftb)
	calc_dftb
	;;
    *)
	calc_reference
	calc_dftb
	;;
esac


