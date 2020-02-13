#!/bin/bash

module load dftbaby

##### HYDROCARBONS ########

# The fit paths are the same as in  JCTC, 2011, 7, 2654-2664

# methane molecule with its central atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/methane.xyz FIT_PATHS/methane_dislocated.xyz dislocate 1 0.75 5 1

# an ethane molecule with one carbon atom displaced on 10 equidistant shells within a 0.75 Ang diameter
perturb_geometry.py STRUCTURES/ethane.xyz FIT_PATHS/ethane_dislocated.xyz dislocate 4 0.75 10 1

#  an ethene molecule with one carbon atom displaced on five shells in a 0.75 Ang diameter sphere
perturb_geometry.py STRUCTURES/ethene.xyz FIT_PATHS/ethene_dislocated.xyz dislocate 3 0.75 5 1

#  an ethine molecule with one carbon atom displaced on five shells in a 0.75 Ang diameter sphere
perturb_geometry.py STRUCTURES/ethine.xyz FIT_PATHS/ethine_dislocated.xyz dislocate 2 0.75 5 1

# a benzene ring with one of its carbon atoms displaced on five equidistant shells within a 0.75 Ang diameter sphere
perturb_geometry.py STRUCTURES/benzene.xyz FIT_PATHS/benzene_dislocated.xyz dislocate 1 0.75 5 1

# hydrogen molecule with the H-H bond stretched in 30 equidistant steps from 0.35 to 2.0 Ang.
scan_internal_dgf.py STRUCTURES/h2.xyz 0 --out_xyz=FIT_PATHS/h2_stretch.xyz --scan_range="(0.35, 2.0, 30)"

#  an ethane molecule with the C-C bond stretched in 30 equidistant steps from 1.0 Ang to 2.5 Ang.
scan_internal_dgf.py STRUCTURES/ethane.xyz 6 --out_xyz=FIT_PATHS/ethane_CC_stretch.xyz --scan_range="(1.0, 2.5, 10)"

#  an ethene molecule with the C-C bond stretched in 30 equidistant steps from 0.4 Ang to 2.0 Ang.
scan_internal_dgf.py STRUCTURES/ethene.xyz 3 --out_xyz=FIT_PATHS/ethene_CC_stretch.xyz --scan_range="(0.4, 2.0, 10)"

#  an ethine molecule with the C-C bond stretched in 30 equidistant steps from 0.2 Ang to 2.0 Ang.
scan_internal_dgf.py STRUCTURES/ethine.xyz 1 --out_xyz=FIT_PATHS/ethine_CC_stretch.xyz --scan_range="(0.4, 2.0, 10)"


# In order to weight the equilibrium geometry more strongly, we generate slightly perturbed
# equilibrium geometries.

# methane molecule with its central atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/methane.xyz FIT_PATHS/methane_equilib_dislocated.xyz dislocate 1 0.2 5 5

# an ethane molecule with one carbon atom displaced on 10 equidistant shells within a 0.2 Ang diameter
perturb_geometry.py STRUCTURES/ethane.xyz FIT_PATHS/ethane_equilib_dislocated.xyz dislocate 4 0.2 10 5

#  an ethene molecule with one carbon atom displaced on five shells in a 0.2 Ang diameter sphere
perturb_geometry.py STRUCTURES/ethene.xyz FIT_PATHS/ethene_equilib_dislocated.xyz dislocate 3 0.2 5 5

#  an ethine molecule with one carbon atom displaced on five shells in a 0.2 Ang diameter sphere
perturb_geometry.py STRUCTURES/ethine.xyz FIT_PATHS/ethine_equilib_dislocated.xyz dislocate 2 0.2 5 5

# a benzene ring with one of its carbon atoms displaced on five equidistant shells within a 0.2 Ang diameter sphere
perturb_geometry.py STRUCTURES/benzene.xyz FIT_PATHS/benzene_equilib_dislocated.xyz dislocate 1 0.2 5 5

##### OXYGEN-CONTAINING CH-COMPOUNDS ####

### O-H, O-O

# water molecule with its central oxygen randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/water.xyz FIT_PATHS/water_equilib_dislocated.xyz dislocate 2 0.2 5 5
# water molecule with its central oxygen randomly displaced on 5 shells within a sphere of diameter 0.55 Ang
perturb_geometry.py STRUCTURES/water.xyz FIT_PATHS/water_dislocated.xyz dislocate 2 0.55 5 1

# hydrogen-peroxide molecule with one oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/hydrogen_peroxide.xyz FIT_PATHS/hydrogen_peroxide_equilib_dislocated.xyz dislocate 2 0.2 5 5
# hydrogen-peroxide molecule with one oxygen atom randomly displaced on 10 shells within a sphere of diameter 0.8 Ang
perturb_geometry.py STRUCTURES/hydrogen_peroxide.xyz FIT_PATHS/hydrogen_peroxide_dislocated.xyz dislocate 2 0.8 10 1

# Zundel cation with central hydrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/zundel_cation.xyz FIT_PATHS/zundel_cation_equilib_dislocated.xyz dislocate 4 0.2 5 5
# hydrogen-peroxide molecule with one oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.55 Ang
perturb_geometry.py STRUCTURES/zundel_cation.xyz FIT_PATHS/zundel_cation_dislocated.xyz dislocate 4 0.55 5 1

### C-O-H

# methanol with oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/methanol.xyz FIT_PATHS/methanol_equilib_dislocated.xyz dislocate 2 0.2 5 5
# methanol with oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.55 Ang
perturb_geometry.py STRUCTURES/methanol.xyz FIT_PATHS/methanol_dislocated.xyz dislocate 2 0.55 5 1

# ethanol with oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/ethanol.xyz FIT_PATHS/ethanol_equilib_dislocated.xyz dislocate 3 0.2 5 5
# ethanol with oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.55 Ang
perturb_geometry.py STRUCTURES/ethanol.xyz FIT_PATHS/ethanol_dislocated.xyz dislocate 3 0.55 5 1

# ethenol with oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/ethenol.xyz FIT_PATHS/ethenol_equilib_dislocated.xyz dislocate 3 0.2 5 5
# ethenol with oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.55 Ang
perturb_geometry.py STRUCTURES/ethenol.xyz FIT_PATHS/ethenol_dislocated.xyz dislocate 3 0.55 5 1

### C=O

# carbon dioxide with one oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/carbon_dioxide.xyz FIT_PATHS/carbon_dioxide_equilib_dislocated.xyz dislocate 1 0.2 5 1

# carbon dioxide with one oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/carbon_dioxide.xyz FIT_PATHS/carbon_dioxide_dislocated.xyz dislocate 1 0.75 5 1


# formaldehyde with oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/formaldehyde.xyz FIT_PATHS/formaldehyde_equilib_dislocated.xyz dislocate 4 0.2 5 5
# formaldehyde with oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.55 Ang
perturb_geometry.py STRUCTURES/formaldehyde.xyz FIT_PATHS/formaldehyde_dislocated.xyz dislocate 4 0.55 5 1

# acetone with oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/acetone.xyz FIT_PATHS/acetone_equilib_dislocated.xyz dislocate 4 0.2 5 5
# acetone with oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.55 Ang
perturb_geometry.py STRUCTURES/acetone.xyz FIT_PATHS/acetone_dislocated.xyz dislocate 4 0.55 5 1

### OH
### | 
### C=O

# carbonic acid with one OH oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/carbonic_acid.xyz FIT_PATHS/carbonic_acid_equilib_dislocated.xyz dislocate 1 0.2 5 1
# carbonic acid with one OH oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/carbonic_acid.xyz FIT_PATHS/carbonic_acid_dislocated.xyz dislocate 1 0.75 5 1


# acetic acid with the OH oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/acetic_acid.xyz FIT_PATHS/acetic_acid_equilib_dislocated.xyz dislocate 4 0.2 5 5
# acetic acid with the OH oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.55 Ang
perturb_geometry.py STRUCTURES/acetic_acid.xyz FIT_PATHS/acetic_acid_dislocated.xyz dislocate 4 0.55 5 1


### Ether C-O-C

# dimethylether with the oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/dimethylether.xyz FIT_PATHS/dimethylether_equilib_dislocated.xyz dislocate 2 0.2 5 1
# dimethylether with the oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/dimethylether.xyz FIT_PATHS/dimethylether_dislocated.xyz dislocate 2 0.75 5 1


### aromatic rings with oxygen

# furan with the oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/furan.xyz FIT_PATHS/furan_equilib_dislocated.xyz dislocate 1 0.2 5 1
# furan with the oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/furan.xyz FIT_PATHS/furan_dislocated.xyz dislocate 1 0.75 5 1

# anisole with the oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/anisole.xyz FIT_PATHS/anisole_equilib_dislocated.xyz dislocate 7 0.2 5 1
# anisole with the oxygen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/anisole.xyz FIT_PATHS/anisole_dislocated.xyz dislocate 7 0.75 5 1


##### NITROGEN-CONTAINING COMPOUNDS ####

# ammonia with the nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/ammonia.xyz FIT_PATHS/ammonia_equilib_dislocated.xyz dislocate 1 0.2 5 1
# ammonia with the nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/ammonia.xyz FIT_PATHS/ammonia_dislocated.xyz dislocate 1 0.75 5 1


# nitrogen molecule with the N-N bond stretched in 10 equidistant steps from 0.35 to 1.6 Ang.
scan_internal_dgf.py STRUCTURES/n2.xyz 0 --out_xyz=FIT_PATHS/n2_stretch.xyz --scan_range="(0.35, 1.6, 10)"


# nitrite with the nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
#perturb_geometry.py STRUCTURES/nitrite.xyz FIT_PATHS/nitrite_equilib_dislocated.xyz dislocate 1 0.2 5 1
# nitrite with the nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
#perturb_geometry.py STRUCTURES/nitrite.xyz FIT_PATHS/nitrite_dislocated.xyz dislocate 1 0.75 5 1

# nitrate with the nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/nitrate.xyz FIT_PATHS/nitrate_equilib_dislocated.xyz dislocate 1 0.2 5 1

# dinitrogen trioxide with one nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/dinitrogen_trioxide.xyz FIT_PATHS/dinitrogen_trioxide_equilib_dislocated.xyz dislocate 1 0.2 5 1
# dinitrogen trioxide with one nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.40 Ang
perturb_geometry.py STRUCTURES/dinitrogen_trioxide.xyz FIT_PATHS/dinitrogen_trioxide_dislocated.xyz dislocate 1 0.40 5 1

# hydrazine with one nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/hydrazine.xyz FIT_PATHS/hydrazine_equilib_dislocated.xyz dislocate 1 0.2 5 1
# hydrazine with one nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/hydrazine.xyz FIT_PATHS/hydrazine_dislocated.xyz dislocate 1 0.75 5 1

# nitrous acid with nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/nitrous_acid.xyz FIT_PATHS/nitrous_acid_equilib_dislocated.xyz dislocate 1 0.2 5 1
# nitrous acid with nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/nitrous_acid.xyz FIT_PATHS/nitrous_acid_dislocated.xyz dislocate 1 0.75 5 1


### Amines

# methylamine with nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/methylamine.xyz FIT_PATHS/methylamine_equilib_dislocated.xyz dislocate 2 0.2 5 1
# methylamine with nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/methylamine.xyz FIT_PATHS/methylamine_dislocated.xyz dislocate 2 0.75 5 1

# dimethylamine with nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/dimethylamine.xyz FIT_PATHS/dimethylamine_equilib_dislocated.xyz dislocate 2 0.2 5 1
# methylamine with nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/dimethylamine.xyz FIT_PATHS/dimethylamine_dislocated.xyz dislocate 2 0.75 5 1

# trimethylamine with nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/trimethylamine.xyz FIT_PATHS/trimethylamine_equilib_dislocated.xyz dislocate 2 0.2 5 1
# methylamine with nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/trimethylamine.xyz FIT_PATHS/trimethylamine_dislocated.xyz dislocate 2 0.75 5 1


### Aminoacids and peptides

# glycine with nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/glycine.xyz FIT_PATHS/glycine_equilib_dislocated.xyz dislocate 3 0.2 5 1
# glycine with nitrogen atom randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/glycine.xyz FIT_PATHS/glycine_dislocated.xyz dislocate 3 0.75 5 1

# nitrogen involved in peptide bond randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/dipeptide.xyz FIT_PATHS/dipeptide_equilib_dislocated.xyz dislocate 3 0.2 5 1
# nitrogen involved in peptide bond randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/dipeptide.xyz FIT_PATHS/dipeptide_dislocated.xyz dislocate 3 0.75 5 1


### N-Heteroaromatics and rings

# central nitrogen in azidobenzene randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/azidobenzene.xyz FIT_PATHS/azidobenzene_equilib_dislocated.xyz dislocate 8 0.2 5 1
# central nitrogen in azidobenzene randomly displaced on 5 shells within a sphere of diameter 0.55 Ang
perturb_geometry.py STRUCTURES/azidobenzene.xyz FIT_PATHS/azidobenzene_dislocated.xyz dislocate 8 0.55 5 1

# nitrogen in nitrobenzene randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/nitrobenzene.xyz FIT_PATHS/nitrobenzene_equilib_dislocated.xyz dislocate 7 0.2 5 1
# nitrogen in nitrobenzene randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/nitrobenzene.xyz FIT_PATHS/nitrobenzene_dislocated.xyz dislocate 7 0.75 5 1

# nitrogen in pyridine randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/pyridine.xyz FIT_PATHS/pyridine_equilib_dislocated.xyz dislocate 1 0.2 5 1
# nitrogen in pyridine randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/pyridine.xyz FIT_PATHS/pyridine_dislocated.xyz dislocate 1 0.75 5 1

# one nitrogen in pyrimidine randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/pyrimidine.xyz FIT_PATHS/pyrimidine_equilib_dislocated.xyz dislocate 1 0.2 5 1
# one nitrogen in pyrimidine randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/pyrimidine.xyz FIT_PATHS/pyrimidine_dislocated.xyz dislocate 1 0.75 5 1

# nitrogen in pyrrole randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/pyrrole.xyz FIT_PATHS/pyrrole_equilib_dislocated.xyz dislocate 1 0.2 5 1
# nitrogen in pyrrole randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/pyrrole.xyz FIT_PATHS/pyrrole_dislocated.xyz dislocate 1 0.75 5 1

# nitrogen in piperidine randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/piperidine.xyz FIT_PATHS/piperidine_equilib_dislocated.xyz dislocate 4 0.2 5 1
# nitrogen in piperidine randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/piperidine.xyz FIT_PATHS/piperidine_dislocated.xyz dislocate 4 0.75 5 1


### C=N double bond

# nitrogen in ethanimine randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/ethanimine.xyz FIT_PATHS/ethanimine_equilib_dislocated.xyz dislocate 4 0.2 5 1
# nitrogen in ethanimine randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/ethanimine.xyz FIT_PATHS/ethanimine_dislocated.xyz dislocate 4 0.75 5 1

# nitrogen in propan-2-imine randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/propan-2-imine.xyz FIT_PATHS/propan-2-imine_equilib_dislocated.xyz dislocate 4 0.2 5 1
# nitrogen in ethanimine randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/propan-2-imine.xyz FIT_PATHS/propan-2-imine_dislocated.xyz dislocate 4 0.75 5 1


### Nitriles

# nitrogen in acetonitrile randomly displaced on 5 shells within a sphere of diameter 0.2 Ang
perturb_geometry.py STRUCTURES/acetonitrile.xyz FIT_PATHS/acetonitrile_equilib_dislocated.xyz dislocate 4 0.2 5 1
# nitrogen in acetonitrile randomly displaced on 5 shells within a sphere of diameter 0.75 Ang
perturb_geometry.py STRUCTURES/acetonitrile.xyz FIT_PATHS/acetonitrile_dislocated.xyz dislocate 4 0.75 5 1
