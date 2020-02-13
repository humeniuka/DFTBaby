#!/bin/bash
#
#

echo   " "
echo   "                                                      E(products) - E(educts) [kcal/mol]  "
echo   " "
echo   "         Reaction                                        Ref.                  DFTB"
echo   " "

# HYDROCARBONS

./reaction_energy.py 1 ethane  1 h2  -2 methane
./reaction_energy.py 1 ethene  2 h2  -2 methane
./reaction_energy.py 1 ethine  3 h2  -2 methane
./reaction_energy.py 1 benzene 9 h2  -6 methane


echo " "
