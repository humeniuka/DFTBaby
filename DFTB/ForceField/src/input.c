/*
Reads the geometry, atom types and lattice vectors
from a file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "input.h"
#include "dreiding.h"

void syntax_error(const char *msg, int line) {
  fprintf(stderr, "Syntax error in line %d\n", line);
  exit(EXIT_FAILURE);
}

#define abort() {  \
      free(atoms); \
      fprintf(stderr, "Syntax error in line %d of force field definition!\n", line); \
      return NULL; \
    }

ATOM *read_force_field(FILE *fh,
		       int *nat, CRYSTAL *crystal)
{
  char buf[2048];
  char atname[20];
  double x,y,z;
  int i, line;
  ATOM *atoms;

  // line number
  line = 1;
  // number of atoms
  if (fscanf(fh, "%d\n", nat) != 1) {
    fprintf(stderr, "Syntax error in line %d of force field definition!\n", line);
    return NULL;
  }
  // read geometries and atom types
  atoms = (ATOM *) malloc(sizeof(ATOM)* *nat);
  // skip comment
  line++;
  if (!fgets(buf, 2048, fh)) abort();

  for(i=0; i<*nat; i++) {
    line++;
    if (fscanf(fh, "%s %lf %lf %lf %d\n", atname, &x, &y, &z, &(atoms[i].atom_type)) != 5) abort();
    atoms[i].x = x * a2b;
    atoms[i].y = y * a2b;
    atoms[i].z = z * a2b;
  }
  // read lattice vectors
  line++;
  if (fscanf(fh, "Tv %lf %lf %lf\n", &x, &y, &z) != 3) abort();
  crystal->ax = x * a2b;
  crystal->ay = y * a2b;
  crystal->az = z * a2b;
  line++;
  if (fscanf(fh, "Tv %lf %lf %lf\n", &x, &y, &z) != 3) abort();
  crystal->bx = x * a2b;
  crystal->by = y * a2b;
  crystal->bz = z * a2b;
  line++;
  if (fscanf(fh, "Tv %lf %lf %lf\n", &x, &y, &z) != 3) abort();
  crystal->cx = x * a2b;
  crystal->cy = y * a2b;
  crystal->cz = z * a2b;

  return atoms;
}

void assign_parameters(ATOM *atoms, int nat, int verbose) {
  int itype;
  int i;

  /* parameters such as `dreiding_valR0` are defined in dreiding.h */
  for(i=0; i<nat; i++) {
    itype = atoms[i].atom_type;
    atoms[i].hyb = dreiding_hyb[itype];
    atoms[i].hDA = dreiding_hDA[itype];
    atoms[i].group = periodic_group[itype];
    atoms[i].covR0 = covalent_radius[itype] * a2b;
    atoms[i].valR0 = dreiding_valR0[itype] * a2b;
    atoms[i].Te = dreiding_angle[itype] * M_PI/180.0;
    atoms[i].cosTe = cos(atoms[i].Te);
    atoms[i].vdWR0 = dreiding_vdWR0[itype] * a2b;
    atoms[i].vdWD0 = dreiding_vdWD0[itype] * kcal2Hartree;
    if (fabs(atoms[i].Te - M_PI) < 1.0e-10) {
      atoms[i].linear = 1;
      atoms[i].Ce = K_ANGLE;
    } else {
      atoms[i].linear = 0;
      atoms[i].Ce = K_ANGLE / pow(sin(atoms[i].Te),2);      
    }
    atoms[i].torsion_barrier = dreiding_torsion_barrier[itype];
    atoms[i].torsion_periodicity = dreiding_torsion_periodicity[itype];
    atoms[i].torsion_phi0 = dreiding_torsion_phi0[itype];
    
  }
}
