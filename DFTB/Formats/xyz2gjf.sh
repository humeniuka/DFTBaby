#!/bin/bash
#
# convert geometry in xyz-format to a simple Gaussian input file
# that can be opened with GaussView
#
# Usage:
#
#     xyz2gjf.sh    input.xyz   output.gjf
#

if [ ! -f "$1" ]
then
    echo "No input file $1!"
    echo " "
    echo "  Usage: "
    echo "         $(basename $0)  molecule.xyz  script.gjf"
    echo " "
    echo "    creates a Gaussian input file 'script.gjf' using the geometry"
    echo "    in 'molecule.xyz'. "
    echo " "
    exit
fi

# input geometry in xyz format
xyz=$1
name=$(basename $xyz .xyz)
# If no output file is given, we use the same filename but with .gjf suffix.
gjf=${2:-${name}.gjf}
# number of atoms
nat=$(head -n 1 $xyz)

# route section, title header, charge (neutral) and multiplicity (1)
cat > $gjf <<EOF
%chk=${name}.chk
# AM1 

${name}

0,1
EOF
# use last geometry in xyz-file
tail -n $nat $xyz >> $gjf
echo -e "\n\n" >> $gjf

echo "Gaussian input script written to '$gjf'."
