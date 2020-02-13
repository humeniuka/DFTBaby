/*

 */

#include <stdio.h>

#include "atoms.h"

// conversion factors
#define a2b 1.8897259885789233
#define kcal2Hartree 0.0015936             // Hartree / (kcal/mol)
#define forceConst2Hartree 0.0004462535    // [Hartree / bohr^2] / [(kcal/mol) / Ang^2]

//
#define RE_DELTA (0.01 * a2b)                           // 0.01 Angstrom
#define KE_SINGLE_BOND  (700.0 * forceConst2Hartree)    // 700.00 (kcal/mol)/Angstrom^2
#define K_ANGLE (100.0 * kcal2Hartree)   // 100 (kcal/mol) / rad^2

#define VbarrierD  (25.0 * kcal2Hartree)

#define K_INVERSION (40.0 * kcal2Hartree)
#define PSI0_INVERSION 0.0
#define COS_PSI0_INVERSION cos(PSI0_INVERSION)
#define NONBONDING_CUTOFF (9.0 * a2b)   // 9.0 Angstrom

// The index of the atom type H__HB in AtomType enumeration in `dreiding.h`
#define HYDROGEN_TYPE 1
// parameters for hydrogen bonds from table V (in Mayo et.al., J.Phys.Chem. vol. 94, 26, 1990)
#define DREIDING_HB_RADIUS  (2.75 * a2b)
//#define DREIDING_HB_RADIUS  (4.0 * a2b)
#define DREIDING_HB_STRENGTH (9.0 * kcal2Hartree)

#define HBOND_CUTOFF (5.5 * a2b)

#define EXCITON_CUTOFF (1000.0 * a2b)
  
#define dreiding_nr_atom_types 37

ATOM *read_force_field(FILE *fh, int *nat, CRYSTAL *crystal);
void assign_parameters(ATOM *atoms, int nat, int verbose);
