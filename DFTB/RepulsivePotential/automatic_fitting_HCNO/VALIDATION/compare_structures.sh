#!/bin/bash

module load dftbaby

echo   " "
echo   "                Ref.      DFTB"
echo   " "

# hydrogen molecule
echo   "   hydrogen molecule"
echo   "   ================="
printf "   r(H-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/h2.xyz "1 2") $(show_molcoords.py -n DFTB/h2.xyz "1 2")
echo   " "

# HYDROCARBONS

echo   "   methane"
echo   "   ======="
printf "   r(C-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/methane.xyz "1 2") $(show_molcoords.py -n DFTB/methane.xyz "1 2")
printf "   <(H-C-H):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/methane.xyz "1 2 3") $(show_molcoords.py -n DFTB/methane.xyz "1 2 3")
echo   " "

echo   "   ethane"
echo   "   ======"
printf "   r(C-C):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/ethane.xyz "1 2") $(show_molcoords.py -n DFTB/ethane.xyz "1 2")
printf "   r(C-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/ethane.xyz "1 3") $(show_molcoords.py -n DFTB/ethane.xyz "1 3")
printf "   <(H-C-H):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/ethane.xyz "1 4 5") $(show_molcoords.py -n DFTB/ethane.xyz "1 4 5")
echo   " "

echo   "   ethene"
echo   "   ======"
printf "   r(C=C):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/ethene.xyz "1 2") $(show_molcoords.py -n DFTB/ethene.xyz "1 2")
printf "   r(C-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/ethene.xyz "1 3") $(show_molcoords.py -n DFTB/ethene.xyz "1 3")
printf "   <(H-C-H):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/ethene.xyz "1 3 4") $(show_molcoords.py -n DFTB/ethene.xyz "1 3 4")
echo   " "

echo   "   ethine"
echo   "   ======"
printf "   r(C-C):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/ethine.xyz "1 2") $(show_molcoords.py -n DFTB/ethine.xyz "1 2")
printf "   r(C-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/ethine.xyz "1 3") $(show_molcoords.py -n DFTB/ethine.xyz "1 3")
echo   " "

echo   "   benzene"
echo   "   ======="
printf "   r(C-C):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/benzene.xyz "1 2") $(show_molcoords.py -n DFTB/benzene.xyz "1 2")
printf "   r(C-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/benzene.xyz "1 7") $(show_molcoords.py -n DFTB/benzene.xyz "1 7")
echo   " "

# OXYGEN-CONTAINING COMPOUNDS

echo   "   water"
echo   "   ====="
printf "   r(O-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/water.xyz "1 2"  ) $(show_molcoords.py -n DFTB/water.xyz "1 2")
printf "   <(H-O-H):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/water.xyz "2 1 3") $(show_molcoords.py -n DFTB/water.xyz "2 1 3")
echo   " "

echo   "   methanol"
echo   "   ========"
printf "   r(C-O):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/methanol.xyz "1 2"  ) $(show_molcoords.py -n DFTB/methanol.xyz "1 2")
printf "   r(O-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/methanol.xyz "2 6"  ) $(show_molcoords.py -n DFTB/methanol.xyz "2 6")
printf "   <(C-O-H):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/methanol.xyz "2 1 6") $(show_molcoords.py -n DFTB/methanol.xyz "2 1 6")
echo   " "

echo   "   formaldehyde"
echo   "   ============"
printf "   r(C=O):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/formaldehyde.xyz "2 4"  ) $(show_molcoords.py -n DFTB/formaldehyde.xyz "2 4")
printf "   <(H-C=O):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/formaldehyde.xyz "2 1 4") $(show_molcoords.py -n DFTB/formaldehyde.xyz "2 1 4")
echo   " "

echo   "   acetic acid"
echo   "   ==========="
printf "   r(C-C):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/acetic_acid.xyz "1 2"  ) $(show_molcoords.py -n DFTB/acetic_acid.xyz "1 2")
printf "   r(C=O):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/acetic_acid.xyz "2 3"  ) $(show_molcoords.py -n DFTB/acetic_acid.xyz "2 3")
printf "   r(C-OH):     %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/acetic_acid.xyz "2 4"  ) $(show_molcoords.py -n DFTB/acetic_acid.xyz "2 4")
printf "   r(O-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/acetic_acid.xyz "4 8"  ) $(show_molcoords.py -n DFTB/acetic_acid.xyz "4 8")
printf "   <(C-O-H):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/acetic_acid.xyz "4 2 8") $(show_molcoords.py -n DFTB/acetic_acid.xyz "4 2 8")
echo   " "

echo   "   dimethylether"
echo   "   ============="
printf "   r(C-O):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/dimethylether.xyz "1 2"  ) $(show_molcoords.py -n DFTB/dimethylether.xyz "1 2")
printf "   <(C-O-C):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/dimethylether.xyz "2 1 3") $(show_molcoords.py -n DFTB/dimethylether.xyz "2 1 3")
echo   " "

echo   "   furan"
echo   "   ====="
printf "   r(C-O):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/furan.xyz "1 2"  ) $(show_molcoords.py -n DFTB/furan.xyz "1 2")
printf "   <(C-O-C):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/furan.xyz "1 2 5") $(show_molcoords.py -n DFTB/furan.xyz "1 2 5")
echo   " "


# NITROGEN-CONTAINING COMPOUNDS

echo   "   ammonia"
echo   "   ======="
printf "   r(N-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/ammonia.xyz "1 2"  ) $(show_molcoords.py -n DFTB/ammonia.xyz "1 2")
printf "   <(H-N-H):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/ammonia.xyz "2 1 3") $(show_molcoords.py -n DFTB/ammonia.xyz "2 1 3")
echo   " "

echo   "   dimethylamine"
echo   "   ============="
printf "   r(N-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/dimethylamine.xyz "2 6"  ) $(show_molcoords.py -n DFTB/dimethylamine.xyz "2 6")
printf "   r(C-N):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/dimethylamine.xyz "1 2"  ) $(show_molcoords.py -n DFTB/dimethylamine.xyz "1 2")
printf "   <(C-N-H):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/dimethylamine.xyz "2 1 6") $(show_molcoords.py -n DFTB/dimethylamine.xyz "2 1 6")
echo   " "

echo   "   nitrobenzene"
echo   "   ============"
printf "   r(N=O):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/nitrobenzene.xyz "7 14"   )  $(show_molcoords.py -n DFTB/nitrobenzene.xyz "7 14")
printf "   r(C-N):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/nitrobenzene.xyz "1 7"    )  $(show_molcoords.py -n DFTB/nitrobenzene.xyz "1 7" )
printf "   <(O=N-O):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/nitrobenzene.xyz "7 13 14")  $(show_molcoords.py -n DFTB/nitrobenzene.xyz "7 13 14")
printf "   <(C-C-N-O):  %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/nitrobenzene.xyz "6 1 7 13") $(show_molcoords.py -n DFTB/nitrobenzene.xyz "6 1 7 13")
echo   " "

echo   "   dipeptide (peptide bond)"
echo   "   ========================"
printf "   r(N-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/dipeptide.xyz "3 12"   )  $(show_molcoords.py -n DFTB/dipeptide.xyz "3 12")
printf "   r(C-N):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/dipeptide.xyz "3 2"    )  $(show_molcoords.py -n DFTB/dipeptide.xyz "3 2" )
printf "   r(C=O):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/dipeptide.xyz "2 9")  $(show_molcoords.py -n DFTB/dipeptide.xyz "2 9")
printf "   <(H-N-C=O):  %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/dipeptide.xyz "12 3 2 9") $(show_molcoords.py -n DFTB/dipeptide.xyz "12 3 2 9")
echo   " "

echo   "   acetonitrile"
echo   "   ============"
printf "   r(C-C):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/acetonitrile.xyz "1 2"  ) $(show_molcoords.py -n DFTB/acetonitrile.xyz "1 2")
printf "   r(C~N):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/acetonitrile.xyz "2 3"  ) $(show_molcoords.py -n DFTB/acetonitrile.xyz "2 3")
echo   " "

echo   "   ethanimine"
echo   "   =========="
printf "   r(C=N):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/ethanimine.xyz "2 4"    ) $(show_molcoords.py -n DFTB/ethanimine.xyz "2 4")
printf "   r(N-H):      %6.4f    %6.4f\n" $(show_molcoords.py -n REFERENCE/ethanimine.xyz "4 8"    ) $(show_molcoords.py -n DFTB/ethanimine.xyz "4 8")
printf "   <(C-N-H):    %6.1f    %6.1f\n" $(show_molcoords.py -n REFERENCE/ethanimine.xyz "4 2 8"  ) $(show_molcoords.py -n DFTB/ethanimine.xyz "4 2 8")
echo   " "
