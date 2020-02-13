#!/bin/bash

if [ $# -lt 1 ]
then
    echo "Usage: $0 <xyz-file>"
    echo " optimize the geometry using Gaussian's AM1 method"
    exit 85
fi

xyz_file=$1
# name of molecule
name=$(basename $xyz_file .xyz)
# number of atoms
nat=$(head -n 1 $xyz_file)
# extract cartesian coordinates
cartesian=$(tail -n $nat $xyz_file)
# create input for Gaussian and run g09
echo "running g09..."
g09 &> ${name}.g09.out <<-EOF
%Nproc=4
%Chk=${name}.chk
#T AM1 Opt Freq

optimize $name

0,1
$cartesian


 
EOF
echo "...done"



