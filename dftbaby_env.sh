#!/bin/bash
#
##############################################################
# This bash script sets the environment variables and paths
# necessary for running DFTBaby. You can source this script
# from your .bashrc profile.
############################################################## 

# libxc-library is needed for pseudoatom calculations, you need to modify this path
export LD_LIBRARY_PATH=/local_scratch/humeniuka/software/qc/libxc/lib:$LD_LIBRARY_PATH

export PYTHONUNBUFFERED=1
# increase stack size
export OMP_STACKSIZE=1g
ulimit -s unlimited

# 'root' has to be set to the DFTBaby top directory.
# The following one-liner guesses the installation directory
root="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Alternatively you can hardcode the path
#root=$HOME/DFTB-0.1.0

export PYTHONPATH=${root}:$PYTHONPATH
# 
export PATH=${root}/DFTB/:${root}/DFTB/Analyse:${root}/DFTB/Analyse/mayavi:${root}/DFTB/Modeling:${root}/DFTB/Dynamics:${root}/DFTB/Dynamics/Analyse:${root}/DFTB/Dynamics/Analyse/Viewer:${root}/DFTB/Format:${root}/DFTB/Scattering:${root}/DFTB/ForceField:${root}/DFTB/MetaDynamics:${root}/DFTB/Poisson:${root}/DFTB/MolecularIntegrals:${root}/DFTB/MultipleScattering:${root}/DFTB/Optimize:$PATH

