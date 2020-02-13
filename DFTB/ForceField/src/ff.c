#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "linked_list.h"
#include "input.h"

typedef struct FORCE_FIELD {
  int nat;
  int nat_repl;
  int ncells; // total number of cells
  CRYSTAL *crystal;
  
  ATOM *atoms;
  REPLICA_ATOM *atoms_repl;
  
  LIST *bond_list;
  LIST *angle_list;
  LIST *torsion_list;
  LIST *inversion_list;
  LIST *nonbonded_list;
  LIST *coulomb_list;
  LIST *hbond_list;

  // exciton couplings
  int state; // current electronic state
  int nchromo; // number of chromophores
  CHROMOPHORE *chromophores;
  LIST *exciton_list;
  int nh;           // dimension of Hamiltonian matrix
  double *Hmatrix;  // hamiltonian matrix, nh*nh elements, row-major order;
  double *Hmatrix_diag; // diagonal elements of Hmatrix, nh elements
  double *exciton_energies; // eigen energies of excitonic Hamiltonian
  double *elec_tdip; // electric transition dipoles (between S0 and exciton states)
  double *magn_tdip; // magnetic transition dipoles (between S0 and exciton states)
  // for LAPACK
  int lwork;
  double *work;
  
} FORCE_FIELD;


typedef struct CONNECTIVITY {
  int ncon; // number of atoms this atom is bonded to
  int partners[8]; // indeces of atoms this atom is connected to,
  // It is assumed that each atom can have at most 8 bonding partners
} CONNECTIVITY;

//----------- Force Field Terms -----------------------

// harmonic bond E = 1/2 Ke (R - Re)^2
typedef struct BOND {
  int I,J;
  double order;   // bond order
  double Re;   // equilibrium bond length
  double Ke;   // harmonic force constant
  // energy due to this bond
  double energy;
  // gradient of energy on atom I, gradient on atom J is -gradI
  double gradI[3];
  // value of bond length
  double r; 
  // derivative of bond length w/r/t atom I, 
  double derI[3]; 
} BOND;

// harmonic cosine angle, E = 1/2 Ce (cos(theta) - cos(Te))^2
typedef struct ANGLE {
  int I,J,K;
  int linear;
  double Ce;  // force constant
  double Te;  // equilibrium angle
  double cosTe;  // cosine of equilibrium angle
  // energy due to angle bend
  double energy;
  // gradients of energy on atoms I, J and K
  double gradI[3];
  double gradJ[3];
  double gradK[3];
  // value of cosine of angle
  double cosTh; 
  // derivative of cosine of angle w/r/t atoms I, J and K
  double derI[3];
  double derJ[3];
  double derK[3];
} ANGLE;

// torsion potential E = 1/2 V_JK (1 - cos(n*(phi - phi^0_JK)))
typedef struct TORSION {
  int I,J,K,L;
  int ntors;  // number of torsions around J-K bond
  // 
  int n;
  double Phi0;
  double cosPhi0;  // cosine of phi0
  double Vbarrier;
  // energy due to this torsion
  double energy;
  // gradients of energy on atoms involved
  double gradI[3];
  double gradJ[3];
  double gradK[3];
  double gradL[3];
  // value of cosine of angle
  double cosPh;
  // derivatives of cosine of torsion angle w/r/t atoms involved
  double derI[3];
  double derJ[3];
  double derK[3];
  double derL[3];
} TORSION;

// spectroscipy inversion 
typedef struct INVERSION {
  int I,J,K,L;  // atom I is bonded to three other atoms J,K,L
  double Ps0;   // inversion angle
  double cosPs0;
  double Kinv;   // force constant
  int planar;   // 1 if Psi=0
  // energy due to this improper torsion
  double energy;
  // gradients on atoms involved
  double gradI[3];
  double gradJ[3];
  double gradK[3];
  double gradL[3];
  // value of sine of angle between bond I-->L and the plane I-J-K
  double sinPs;
  // derivatives of sine of inversion angle w/r/t atoms involved
  double derI[3];
  double derJ[3];
  double derK[3];
  double derL[3];
} INVERSION;

// 12-6 Lennard Jones form   E = A R^(-12) - B R^(-6)
typedef struct VAN_DER_WAALS {
  int I,J;
  double A;   // A = D0 * R0^12
  double B;   // B = D0 * R0^6
  // energy due to this vdW interaction
  double energy;
  // gradient on atom I, gradient on atom J is -gradI
  double gradI[3];

} VAN_DER_WAALS;

// Coulomb interaction E = qq / R
typedef struct COULOMB {
  int I,J;
  // product of partial charges
  double qq;
  // energy due to ionic interaction
  double energy;
  // gradient on atom I, gradient on atom J is -gradJ
  double gradI[3];
} COULOMB;

// Hydrogen bonds
/*

   D-----H
       th \
           \
            \
            A
                 Rhb  12          Rhb   10                     4
  E = Dhb [ 5 (------)     -  6 (------)   ]  cos(theta(D-H-A))
               R(D-A)            R(D-A)

*/
typedef struct HBOND {
  int I,J,K;   // indeces of  donor, hydrogen and acceptor
  double Dhb;
  double Rhb;
  // energy due to hydrogen bond
  double energy;
  // gradients on donor, hydrogen and acceptor
  double gradI[3];
  double gradJ[3];
  double gradK[3];
} HBOND;

// Exciton coupling between transition charges
/*
          qI qJ
   V = -------------
       |r(I) - r(J)|
 */
typedef struct EXCITON {
  int I,J; // indeces of atoms on chromophore A and chromophore B
  double qq;  // product of partial transition charges
  int elem1, elem2; // row and column index into the Hamiltonian matrix
  int elem; // 1D index into the flattened Hamiltonian matrix where this coupling should be added to
  // coupling
  double energy;
  // gradient on atom I
  double gradI[3];
} EXCITON;
//------------------------------------------------------

// This definition is a copy of the one in 'dreiding.h'
enum Hybridization {
  // spR: resonance situation
  // monoval: monovalent atoms such as H, Cl, Fl
  // ionic: ionic bond for atoms such as Na^+, Ca^2+, Fe^2+, Zn^2+
  monoval, sp1, sp2, sp3, spR, ionic
};


//-----------------------------------------------------------------------
/* LAPACK */
extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
		   double* w, double* work, int* lwork, int* info );
//-----------------------------------------------------------------------

inline void cross_product(double *a, double *b, double *cr) {
  // compute the cross product cr = a x b
  cr[0] = a[1]*b[2] - a[2]*b[1];
  cr[1] = a[2]*b[0] - a[0]*b[2];
  cr[2] = a[0]*b[1] - a[1]*b[0];
}

inline double dot_product(double *a, double *b) {
  double dot;
  dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return dot;
}

inline double triple_product(double *a, double *b, double *c) {
  // compute the triple product a . (b x c)
  double cr[3];
  double res;
  cross_product(b, c, cr);

  res = dot_product(a, cr);
  printf("triple product a . (b x c) = %f\n", res);
  return dot_product(a, cr);
}

inline double norm(double *a) {
  double a2;
  a2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
  return sqrt(a2);
}

FORCE_FIELD *build_force_field(int nat, ATOM *atoms, CRYSTAL *crystal,
			       int nchromo, CHROMOPHORE *chromophores,
			       int *conmat, 
			       int verbose)
{
  FORCE_FIELD *ff;
  int n1,n2,n3; // enumerate replicas of unit cell
  int i,j,k,l; // enumerate atoms, also used as general index
  int N;
  int I,J,K,L,tmp; // indeces of replica atoms
  ATOM *atI, *atJ, *atK;
  REPLICA_ATOM *atIr, *atJr, *atKr, *atLr;
  double xIJ, yIJ, zIJ, rIJ;
  double rKI[3], rKJ[3], cosTh;
  CONNECTIVITY *connectivities;
  CONNECTIVITY *conI, *conJ, *conK;
  int is_connected;
  
  NODE *node;
  
  LIST *bond_list;
  BOND *bond;
  
  LIST *nonbonded_list;
  VAN_DER_WAALS *vdW;

  LIST *coulomb_list;
  COULOMB *coul;

  LIST *hbond_list;
  HBOND *hbond;
  
  LIST *angle_list;
  ANGLE *angle;

  LIST *exciton_list;
  EXCITON *exciton;
  CHROMOPHORE *chromo1, *chromo2;
  int c1,c2, s1,s2, i1,i2;
  double q1, q2;
  int elem1, elem2;
  // for LAPACK
  int info;
  double wkopt;
  
  LIST *torsion_list;
  TORSION *torsion;
  int n, ntors;  // number of torsions
  int Is[12];
  int Ls[12];
  
  LIST *inversion_list;
  INVERSION *inversion;
  
  double Rij, Rij6, Dij;
  int continue_flag;
  
  REPLICA_ATOM *at_r;

  int dim;   // enumerate lattice vectors
  int nc[3]; // number of neigbhour cells in a-, b- and c-directions
  int ncells; // total number of cells
  // HACK: If a lattice vectors is zero, this signals that
  //       no periodic calculation is required in that direction (ncX = 0).
  double norm2[3]; // lengths squared of lattice vectors
  norm2[0] = pow(crystal->ax,2) + pow(crystal->ay,2) + pow(crystal->az,2);
  norm2[1] = pow(crystal->bx,2) + pow(crystal->by,2) + pow(crystal->bz,2);
  norm2[2] = pow(crystal->cx,2) + pow(crystal->cy,2) + pow(crystal->cz,2);
  for (dim=0; dim<3; dim++) {
    if (norm2[dim] == 0.0) {  // lattice vector for direction dim is 0
      nc[dim] = 0;
    } else {
      nc[dim] = 1;  // one replica cell in this direction
    }
  }
  if (verbose > 0) {
    printf("number of neighbour cells: %d %d %d\n", nc[0], nc[1], nc[2]);
  }
  
  ff = (FORCE_FIELD *) malloc(sizeof(FORCE_FIELD));
  ff->state = 0;    // default to ground state
  ff->nat = nat;
  // The cells are replicated only in one direction to avoid the double counting of interactions.
  // Otherwise there would be a bond between the lower and the central cell that is identical
  // to a bond between the central and the upper cell.
  ncells = (nc[0]+1)*(nc[1]+1)*(nc[2]+1);
  ff->ncells = ncells;
  ff->nat_repl = nat*ncells;
  ff->crystal = crystal;
  ff->atoms = atoms;
  ff->atoms_repl = (REPLICA_ATOM *) malloc(sizeof(REPLICA_ATOM) * ff->nat_repl);
  // create all replica atoms
  k = 0;
  for (n1=0; n1<=nc[0]; n1++) {
    for (n2=0; n2<=nc[1]; n2++) {
      for (n3=0; n3<=nc[2]; n3++) {
	// enumerate atoms in central unit cell
	for (i=0; i<nat; i++) {
	  at_r = &(ff->atoms_repl[k]);
	  // at_r is a replica of atom i
	  at_r->index = i;
	  if ((n1 == 0) && (n2 == 0) && (n3 == 0)) {
	    at_r->central = 1;
	  } else {
	    at_r->central = 0;
	  }
	  // relative to the i-th atom in the central unit cell the replica
	  // is displaced by the vector [tx,ty,tz]
	  at_r->tx = n1*crystal->ax + n2*crystal->bx + n3*crystal->cx;
	  at_r->ty = n1*crystal->ay + n2*crystal->by + n3*crystal->cy;
	  at_r->tz = n1*crystal->az + n2*crystal->bz + n3*crystal->cz;
	  
	  k++;
	}

      }
    }
  }
  if (verbose > 0) {
    printf("number of atoms: %d\n", nat);
    printf("number of replica atoms: %d\n", k);
  }
  N = ff->nat_repl;
  // connectivities
  connectivities = (CONNECTIVITY *) malloc(sizeof(CONNECTIVITY) * N);
  // BONDS
  bond_list = list_new();
  for (I=0; I<N; I++) {
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    // initialized connectivity of atom I
    conI = &(connectivities[I]);
    conI->ncon = 0;
    
    for (J=0; J<N; J++) {
      if (I == J) continue;
      atJr = &(ff->atoms_repl[J]);
      atJ = &(ff->atoms[atJr->index]);
      // distance vector
      xIJ = atJ->x - atI->x + atIr->tx - atJr->tx;
      yIJ = atJ->y - atI->y + atIr->ty - atJr->ty;
      zIJ = atJ->z - atI->z + atIr->tz - atJr->tz;
      rIJ = sqrt(xIJ*xIJ + yIJ*yIJ + zIJ*zIJ);
      if (conmat != NULL) {
	// If a connectivity matrix has been provided, we use it ...
	is_connected = (conmat[I*N+J] == 1);
      } else { 
	// ... otherwise the connectivity is determined from covalent radii. 
	is_connected = (rIJ < 1.3*(atI->covR0 + atJ->covR0));
      }
      if (is_connected) {
	// atoms I and J are bonded
	if (conI->ncon < 8) {
	  conI->partners[conI->ncon] = J;
	  conI->ncon++;
	} else {
	  continue;
	}
	// at least one atom should lie in the unit cell
	if (!(atIr->central | atJr->central)) continue;
	// add (I,J) bond to list if I<J
	if (I < J) {
	  bond = (BOND *) malloc(sizeof(BOND));
	  bond->I = I;
	  bond->J = J;
	  // default bond order is 1
	  bond->order = 1.0;
	  // see eqn. (6) in the Dreiding paper
	  bond->Re = atI->valR0 + atJ->valR0 - RE_DELTA;
	  bond->Ke = KE_SINGLE_BOND;
	  list_append(bond_list, (void *) bond);
	  if (verbose > 1) {
	    printf("BOND(%d-%d)\n", I+1,J+1);
	  }
	}
      }
    }
  }
  ff->bond_list = bond_list;
  if (verbose > 0) {
    printf("number of bonds: %d\n", (int) list_length(bond_list));
  }
  // connectivity matrix is not needed anymore
  if (conmat != NULL) {
    free(conmat);
  }
  //printf("assign bond orders\n");
  // TODO: better bond order assignment
  //       The force constant bond->Ke has to change with the bond order
  for(node=ff->bond_list->head; node!=NULL; node=node->next) {
    bond = (BOND *) node->element;
    I = bond->I; J = bond->J;
    atIr = &(ff->atoms_repl[I]);
    atJr = &(ff->atoms_repl[J]);
    conI = &(connectivities[I]);
    conJ = &(connectivities[J]);
    atI = &(ff->atoms[atIr->index]);
    atJ = &(ff->atoms[atJr->index]);
    // These assignments can be wrong!
    if ((atI->hyb == sp2) && (atJ->hyb == sp2)) {
      bond->order = 2.0;
      bond->Ke *= 2.0;      // double bond is twice as strong as single bond
    } else if ((atI->hyb == sp1) && (atJ->hyb == sp1)) {
      bond->order = 3.0;
      bond->Ke *= 3.0;      // triple bond is three times as strong as single bond
    } else if ((atI->hyb == spR) && (atJ->hyb == spR)) {
      bond->order = 1.5;
    } else {
      // all other bond orders remain 1
    }
    //printf("BOND(%d-%d,order=%f)\n", I+1,J+1, bond->order);
  }
  // ANGLES 
  angle_list = list_new();
  for (J=0; J<N; J++) {
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    conJ = &(connectivities[J]);
    for (i=0; i<conJ->ncon; i++) {
      I = conJ->partners[i];
      atIr = &(ff->atoms_repl[I]);
      for (k=i+1; k<conJ->ncon; k++) {
	K = conJ->partners[k];
	if (K<I) {
	  // swap atoms
	  tmp = I;
	  I = K;
	  K = tmp;
	}
	atKr = &(ff->atoms_repl[K]);
	// at least one atom should lie in the unit cell
	if (!(atIr->central | atJr->central | atKr->central)) continue;	
	angle = (ANGLE *) malloc(sizeof(ANGLE));
	angle->I = I;
	angle->J = J;
	angle->K = K;
	angle->Ce = atJ->Ce;
	angle->Te = atJ->Te;
	angle->cosTe = atJ->cosTe;
	angle->linear = atJ->linear;
	list_append(angle_list, (void *) angle);
	if (verbose > 1) {
	  printf("ANGLE(%d-%d-%d)\n", I+1,J+1,K+1);
	}
      }
    }
  }
  ff->angle_list = angle_list;
  if (verbose > 0) {
    printf("number of angles: %d\n", (int) list_length(angle_list));
  }
  // TORSIONS
  torsion_list = list_new();
  for(node=ff->bond_list->head; node!=NULL; node=node->next) {
    bond = (BOND *) node->element;
    J = bond->I; K = bond->J;
    atJr = &(ff->atoms_repl[J]);
    atKr = &(ff->atoms_repl[K]);
    conJ = &(connectivities[J]);
    conK = &(connectivities[K]);
    atJ = &(ff->atoms[atJr->index]);
    atK = &(ff->atoms[atKr->index]);
    // count number of torsions around J-K bond and
    // find the neighbouring bond I-J and K-L
    ntors = 0;
    for(i=0;i<conJ->ncon;i++) {
      I = conJ->partners[i];
      if (I == K) continue;
      for(l=0;l<conK->ncon;l++) {
	L = conK->partners[l];
	if (L == J) continue;
	atIr = &(ff->atoms_repl[I]);
	atLr = &(ff->atoms_repl[L]);
	atI = &(ff->atoms[atIr->index]);
	// at least one atom should lie in the unit cell
	if (!(atIr->central | atJr->central | atKr->central | atLr->central)) continue;	
	if (ntors < 12) {
	  Is[ntors] = I;
	  Ls[ntors] = L;
	  ntors++;
	}
      }
    }
    if (ntors > 0) {
      // add all torsions around J-K bond
      for(n=0; n<ntors; n++) {
	torsion = (TORSION *) malloc(sizeof(TORSION));
	torsion->ntors = ntors;
	I = Is[n]; L = Ls[n];
	torsion->I = I;
	torsion->J = J;
	torsion->K = K;
	torsion->L = L;
	// force field parameters
	// distinguish different cases (a)-(j) in the Dreiding paper
	// TODO: not all cases are implemented, also the cases are not mutually exclusive
	//       This 
	if ((bond->order == 1.0) && (atJ->hyb == sp3) && (atK->hyb == sp3) && (atJ->group == 16) && (atK->group == 16)) {
	  // case (h) - A dihedral single bond involving two sp3 atoms of the oxygen column
	  //            (J,K = X_3 of column 16)
	  //printf("torsion case (h)\n");
	  torsion->n = 2;
	  torsion->Phi0 = 90.0;
	  torsion->cosPhi0 = 0.0;
	  torsion->Vbarrier = 2.0 * kcal2Hartree;	  
	} else if ((bond->order == 1.0) && (atJ->hyb == sp3) && (atJ->group == 16) && ((atK->hyb == sp2) || (atK->hyb == spR)) && (atK->group != 16)) {
	  // case (i) - For dihedral bonds involving an sp3 atom of the oxygen column with an sp2 or resonant atom of another column (J = X_3 of column 16, K = X_2, X_R)
	  //printf("torsion case (i)\n");
	  torsion->n = 2;
	  torsion->Phi0 = 180.0;
	  torsion->cosPhi0 = -1.0;
	  torsion->Vbarrier = 2.0 * kcal2Hartree;	  	  
	} else if ((bond->order == 1.0) && (atJ->hyb == sp3) && (atK->hyb == sp3)) {
	  //printf("torsion case (a)\n");
	  // case (a) - dihedral single bond involving two sp3 atoms (J,K = X_3)
	  torsion->n = 3;
	  torsion->Phi0 = 180.0;
	  torsion->cosPhi0 = -1.0;
	  torsion->Vbarrier = 2.0 * kcal2Hartree;
	} else if ((bond->order == 1.0) && ((atJ->hyb == sp2) || (atJ->hyb== spR)) && (atK->hyb == sp3)) {
	  //printf("torsion case (b)\n");
	  // case (b) - dihedral single bond involving one sp2 center and one sp3 center
	  //            (J = X_2,X_R, K=X_3)
	  torsion->n = 6;
	  torsion->Phi0 = 0.0;
	  torsion->cosPhi0 = 1.0;
	  torsion->Vbarrier = 1.0 * kcal2Hartree;	  
	} else if ((bond->order == 2.0) && (atJ->hyb == sp2) && (atK->hyb == sp2)) {
	  //printf("torsion case (c)\n");
	  // case (c) - dihedral double bond involving two sp2 atoms (J,K = X_2)
	  torsion->n = 2;
	  torsion->Phi0 = 0.0;
	  torsion->cosPhi0 = 1.0;
	  torsion->Vbarrier = 45.0 * kcal2Hartree;	  	  
	} else if ((bond->order == 1.5) && (atJ->hyb == spR) && (atK->hyb == spR)) {
	  //printf("torsion case (d)\n");
	  // case (d)  - dihedral resonance bond involving two resonant atoms (J,K = X_R)
	  torsion->n = 2;
	  torsion->Phi0 = 0.0;
	  torsion->cosPhi0 = 1.0;
	  torsion->Vbarrier = 25.0 * kcal2Hartree;
	} else if ((bond->order == 1) && ((atJ->hyb == sp2) || (atJ->hyb == spR)) && ((atK->hyb == sp2) || (atK->hyb == spR))) {
	  //printf("torsions case (e)\n");
	  // case (e)  - dihedral single bond involving two sp2 or resonant atoms
	  //             (e.g. the middle bond of butadiene)
	  torsion->n = 2;
	  torsion->Phi0 = 18.0;
	  torsion->cosPhi0 = -1.0;
	  torsion->Vbarrier = 5.0 * kcal2Hartree;	  
	} else {
	  //printf("torsion case (g)\n");
	  // case (g) For dihedrals involving one or two sp1 atoms (X_1),
	  //          monovalent atoms  (F, Cl, Na, ...) or metals (Fe, Zn, ...)
	  torsion->n = 2;
	  torsion->Phi0 = 0.0; 
	  torsion->cosPhi0 = 1.0; 
	  torsion->Vbarrier = 0.0; 
	}
	// divide by number of torsions
	torsion->Vbarrier /= (double) torsion->ntors;
	list_append(torsion_list, (void *) torsion);
	if (verbose > 1) {
	  printf("TORSION(%d-%d-%d-%d)\n", I+1,J+1,K+1,L+1);
	}
	//
      }
    }
  }
  ff->torsion_list = torsion_list;
  if (verbose > 0) {
    printf("number of torsions: %d\n", (int) list_length(torsion_list));
  }
  // INVERSIONS
  inversion_list = list_new();
  for (I=0; I<N; I++) {
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    conI = &(connectivities[I]);
    if (conI->ncon == 3) {
      // atom I is bonded to exactly three atoms
      J = conI->partners[0]; K = conI->partners[1]; L = conI->partners[2];
      atJr = &(ff->atoms_repl[J]);
      atKr = &(ff->atoms_repl[K]);
      atLr = &(ff->atoms_repl[L]);
      // at least one atom should lie in the unit cell
      if (!(atIr->central | atJr->central | atKr->central | atLr->central)) continue;
      // inversion is only added for I = X_2, X_R
      if (!(atI->hyb == sp2 || atI->hyb == spR)) continue;
      inversion = (INVERSION *) malloc(sizeof(INVERSION));
      inversion->I = I;
      inversion->J = J;
      inversion->K = K;
      inversion->L = L;
      // parameters
      inversion->Ps0 = PSI0_INVERSION;
      inversion->cosPs0 = COS_PSI0_INVERSION;
      inversion->Kinv = K_INVERSION;
      if (fabs(inversion->Ps0 - 0.0) < 1.0e-10) {
	inversion->planar = 1;
      } else {
	inversion->planar = 0;
      }
      list_append(inversion_list, (void *) inversion);
    }
  }
  ff->inversion_list = inversion_list;
  if (verbose > 0) {
    printf("number of inversions: %d\n", (int) list_length(inversion_list));
  }
  // NON-BONDING
  nonbonded_list = list_new();
  for (I=0; I<N; I++) {
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    conI = &(connectivities[I]);
    for (J=I+1; J<N; J++) {
      continue_flag = 0;
      // check if J is bonded to I
      for (i=0; i<conI->ncon; i++) {
	if (J == conI->partners[i]) {
	  continue_flag = 1;
	  break;
	}
      }
      if (continue_flag == 1) continue;
      // check that I and J are not bonded to a common atom
      conJ = &(connectivities[J]);
      for (i=0; i<conI->ncon; i++) {
	for (j=0; j<conJ->ncon; j++) {
	  if (conI->partners[i] == conJ->partners[j]) {
	    continue_flag = 1;
	    break;
	  }
	}
	if (continue_flag == 1) break;
      }
      if (continue_flag == 1) continue;      
      atJr = &(ff->atoms_repl[J]);
      atJ = &(ff->atoms[atJr->index]);
      // at least one atom should lie in the unit cell
      if (!(atIr->central | atJr->central)) continue;
      // distance vector
      xIJ = atJ->x - atI->x + atIr->tx - atJr->tx;
      yIJ = atJ->y - atI->y + atIr->ty - atJr->ty;
      zIJ = atJ->z - atI->z + atIr->tz - atJr->tz;
      rIJ = sqrt(xIJ*xIJ + yIJ*yIJ + zIJ*zIJ);
      if (rIJ > NONBONDING_CUTOFF) continue;
      //
      vdW = (VAN_DER_WAALS *) malloc(sizeof(VAN_DER_WAALS));
      vdW->I = I;
      vdW->J = J;
      Rij = sqrt(atI->vdWR0 * atJ->vdWR0);
      Dij = sqrt(atI->vdWD0 * atJ->vdWD0);
      Rij6 = pow(Rij,6);
      vdW->B = Dij * Rij6;
      vdW->A = Dij * Rij6*Rij6;

      list_append(nonbonded_list, (void *) vdW);
      //printf("VDW(%d,%d)\n", I+1,J+1);
    }
  }

  ff->nonbonded_list = nonbonded_list;
  if (verbose > 0) {
    printf("number of non-bonded interactions: %d\n", (int) list_length(nonbonded_list));
  }

  // COULOMB
  coulomb_list = list_new();
  for (I=0; I<N; I++) {
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    if (atI->charge == 0.0) continue;  // skip neutral atoms I
    conI = &(connectivities[I]);
    for (J=I+1; J<N; J++) {
      // check if J is bonded to I
      continue_flag = 0;
      for (i=0; i<conI->ncon; i++) {
	if (J == conI->partners[i]) {
	  continue_flag = 1;
	  break;
	}
      }
      if (continue_flag == 1) continue;
      
      atJr = &(ff->atoms_repl[J]);
      atJ = &(ff->atoms[atJr->index]);
      if (atJ->charge == 0.0) continue; // skip neutral atoms J
      // at least one atom should lie in the unit cell
      if (!(atIr->central | atJr->central)) continue;
      // distance vector
      xIJ = atJ->x - atI->x + atIr->tx - atJr->tx;
      yIJ = atJ->y - atI->y + atIr->ty - atJr->ty;
      zIJ = atJ->z - atI->z + atIr->tz - atJr->tz;
      rIJ = sqrt(xIJ*xIJ + yIJ*yIJ + zIJ*zIJ);
      if (rIJ > NONBONDING_CUTOFF) continue;
      //
      coul = (COULOMB *) malloc(sizeof(COULOMB));
      coul->I = I;
      coul->J = J;
      coul->qq = atI->charge * atJ->charge;

      list_append(coulomb_list, (void *) coul);
      //printf("COULOMB(%d,%d)\n", I+1,J+1);
    }
  }

  ff->coulomb_list = coulomb_list;
  if (verbose > 0) {
    printf("number of Coulomb interactions: %d\n", (int) list_length(coulomb_list));
  }

  // HYDROGEN BONDS
  hbond_list = list_new();
  for (I=0; I<N; I++) {  // loop over possible hydrogen donors
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    if (atI->hDA == 0) continue;       // skip atoms that cannot donate H-bonds
    conI = &(connectivities[I]);
    // find hydrogen atom K to which donor I is bonded
    continue_flag = 1;
    for (k=0; k<conI->ncon; k++) {
      K = conI->partners[k];
      atKr = &(ff->atoms_repl[K]);
      atK = &(ff->atoms[atKr->index]);
      if (atK->atom_type == HYDROGEN_TYPE) {
	for (J=0; J<N; J++) { // loop over possible hydrogen acceptors
	  if (J == I) continue;
	  atJr = &(ff->atoms_repl[J]);
	  atJ = &(ff->atoms[atJr->index]);
	  if (atJ->hDA == 0) continue;       // skip atoms that cannot accept H-bonds
	  continue_flag = 0;
	  // check if J is bonded to I
	  for (i=0; i<conI->ncon; i++) {
	    if (J == conI->partners[i]) {
	      continue_flag = 1;
	      break;
	    }
	  }
	  if (continue_flag == 1) continue;
	  // check that I and J are not bonded to a common atom
	  conJ = &(connectivities[J]);
	  for (i=0; i<conI->ncon; i++) {
	    for (j=0; j<conJ->ncon; j++) {
	      if (conI->partners[i] == conJ->partners[j]) {
		continue_flag = 1;
		break;
	      }
	    }
	    if (continue_flag == 1) break;
	  }
	  if (continue_flag == 1) continue;      
	  // at least one atom should lie in the unit cell
	  if (!(atIr->central | atJr->central | atKr->central)) continue;
	  // distance vector
	  xIJ = atJ->x - atI->x + atIr->tx - atJr->tx;
	  yIJ = atJ->y - atI->y + atIr->ty - atJr->ty;
	  zIJ = atJ->z - atI->z + atIr->tz - atJr->tz;
	  rIJ = sqrt(xIJ*xIJ + yIJ*yIJ + zIJ*zIJ);
	  if (rIJ > HBOND_CUTOFF) continue;
	  // angle donor--hydrogen--acceptor should be close to 180 deg
	  //   vector  hydrogen->donor
	  rKI[0] = atI->x - atK->x + atKr->tx - atIr->tx;
	  rKI[1] = atI->y - atK->y + atKr->ty - atIr->ty;
	  rKI[2] = atI->z - atK->z + atKr->tz - atIr->tz;
	  //   vector  hydrogen->acceptor
	  rKJ[0] = atJ->x - atK->x + atKr->tx - atJr->tx;
	  rKJ[1] = atJ->y - atK->y + atKr->ty - atJr->ty;
	  rKJ[2] = atJ->z - atK->z + atKr->tz - atJr->tz;
	  // cos(D-H-A)
	  cosTh = dot_product(rKI, rKJ)/( norm(rKI)*norm(rKJ) );
	  // The angle(D-H--A) has to be larger than ~ 143 deg, tolerance 37 deg
	  //if (cosTh > -0.8) continue; //  cosTh = -1  ->  theta = 180 deg
	  // The angle(D-H--A) has to be larger than 120 deg, tolerance 60 deg
	  if (cosTh > -0.5) continue;  // cosTh = -0.5 -> theta = 120
	  // The angle(D-H--A) has to be larger than 90 deg, tolerance 90 deg
	  //if (cosTh > 0.0) continue;  // cosTh = 0 -> theta = 90
	  //
	  hbond = (HBOND *) malloc(sizeof(HBOND));
	  hbond->I = I;
	  hbond->J = K;   // the hydrogen atom is at the central vertex of the angle D-H-A
	  hbond->K = J;
	  hbond->Rhb = DREIDING_HB_RADIUS;
	  hbond->Dhb = DREIDING_HB_STRENGTH;
	  
	  list_append(hbond_list, (void *) hbond);
	  //printf("H-BOND(%d-%d--%d)\n", I+1,K+1,J+1);
	}
      }
    }
  }
  
  ff->hbond_list = hbond_list;
  if (verbose > 0) {
    printf("number of hydrogen bonds: %d\n", (int) list_length(hbond_list));
  }
    
  free(connectivities);

  // EXCITONIC COUPLINGS
  exciton_list = list_new();
  ff->nchromo = nchromo;
  ff->chromophores = chromophores;
  if (nchromo > 0) {
    // determine number of basis set = size of Hamiltonian matrix
    ff->nh = 0;
    for(n1=0; n1<ncells; n1++) { // enumerate cells
      for(c1=0; c1<nchromo; c1++) { // enumerate chromophores
	chromo1 = &(chromophores[c1]);
	for(s1=0; s1<chromo1->nst; s1++) { // enumerate excited state
	  ff->nh += 1;
	}
      }
    }
    // Hamiltonian matrix
    ff->Hmatrix = (double *) malloc(sizeof(double) * ff->nh*ff->nh);
    ff->Hmatrix_diag = (double *) malloc(sizeof(double) * ff->nh);
    ff->exciton_energies = (double *) malloc(sizeof(double) * ff->nh);
    // transition dipoles
    ff->elec_tdip = (double *) malloc(sizeof(double) * ff->nh * 3);
    ff->magn_tdip = (double *) malloc(sizeof(double) * ff->nh * 3);
    
    // query and allocate optimal workspace for Lapack
    ff->lwork = -1;
    dsyev_( "Vectors", "Upper", &ff->nh, ff->Hmatrix, &ff->nh, ff->exciton_energies,
	   &wkopt, &ff->lwork, &info);
    ff->lwork = (int) wkopt;
    ff->work = (double *) malloc(sizeof(double) * ff->lwork);
    
    // build list of excitonic couplings
    elem1 = -1;   // enumerate rows of Hamiltonian matrix
    for(n1=0; n1<ncells; n1++) { // enumerate cells 1
      for(c1=0; c1<nchromo; c1++) { // enumerate chromophores 1
	chromo1 = &(chromophores[c1]);
	for(s1=0; s1<chromo1->nst; s1++) { // enumerate excited states 1
	  elem1 += 1;
	  // set diagonal elements of Hamiltonian matrix
	  ff->Hmatrix_diag[elem1] = chromo1->excitation_energies[s1];
	  for(i1=0; i1<chromo1->nat; i1++) { // enumerate atoms of chromophore 1
	    q1 = chromo1->transition_charges[i1*4*chromo1->nst+4*s1+0];
	    if (q1 == 0.0) continue;
	    // replica atom I on replica chromophore
	    I = n1*nat + chromo1->indeces[i1];
	    atIr = &(ff->atoms_repl[I]);
	    atI = &(ff->atoms[atIr->index]);

	    elem2 = -1;  // enumerate columns of Hamiltonian matrix
	    for(n2=0; n2<ncells; n2++) { // enumerate cells 2
	      for(c2=0; c2<nchromo; c2++) { // enumerate chromophores 2
		chromo2 = &(chromophores[c2]);
		for(s2=0; s2<chromo2->nst; s2++) { // enumerate excited states 2
		  elem2 += 1;
		  // Since the Hamiltonian matrix is symmetric only the LOWER half
		  // has to be filled in
		  if (elem2 >= elem1) continue;
		  for(i2=0; i2<chromo2->nat; i2++) { // enumerate atoms of chromophore 2
		    q2 = chromo2->transition_charges[i2*4*chromo2->nst+4*s2+0];
		    if (q2 == 0.0) continue;
		    // replica atom J on replica chromophore
		    J = n2*nat + chromo2->indeces[i2];
		    atJr = &(ff->atoms_repl[J]);
		    atJ = &(ff->atoms[atJr->index]);
		    // distance vector
		    xIJ = atJ->x - atI->x + atIr->tx - atJr->tx;
		    yIJ = atJ->y - atI->y + atIr->ty - atJr->ty;
		    zIJ = atJ->z - atI->z + atIr->tz - atJr->tz;
		    rIJ = sqrt(xIJ*xIJ + yIJ*yIJ + zIJ*zIJ);
		    if (rIJ > EXCITON_CUTOFF) continue;
		    //
		    exciton = (EXCITON *) malloc(sizeof(EXCITON));
		    exciton->qq = q1*q2;
		    exciton->I = I;
		    exciton->J = J;
		    exciton->elem1 = elem1;
		    exciton->elem2 = elem2;
		    exciton->elem = elem1*ff->nh + elem2;

		    list_append(exciton_list, (void *) exciton);
		    //printf("EXCITON COUPLING V(cell=%d chromophore=%d elec.state=%d atom=%d; cell=%d chromophore=%d elec.state=%d atom=%d ; q1*q2=%f)\n", n1,c1,s1,I, n2,c2,s2,J, q1*q2);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    if (verbose > 0) {
      printf("number of exciton couplings: %d\n", (int) list_length(exciton_list));
      printf("size of Hamiltonian matrix: %d x %d\n", ff->nh, ff->nh);
    }
  } else {
    ff->nh = 0; // only ground state
    ff->Hmatrix = NULL;
    ff->Hmatrix_diag = NULL;
    ff->exciton_energies = NULL;
    ff->work = NULL;
  }
  ff->exciton_list = exciton_list;

  return ff;
}

// The following formulae are taken from Guanghua Gao' PhD thesis
// "Large Scale Molecular Simulations with Application to Polymers and Nano-scale Materials"

inline double angle_gradient(double *A, double *B,
			     double *gradI, double *gradJ, double *gradK) {
  double cosTh;
  double nA, nB;
  double iAB, iA2, iB2;
  int i;
  nA = norm(A);
  nB = norm(B);
  cosTh = dot_product(A,B)/(nA*nB);
  // gradients
  iAB = 1.0/(nA*nB);
  iA2 = 1.0/(nA*nA);
  iB2 = 1.0/(nB*nB);

  for(i=0;i<3;i++) {
    gradI[i] =  iAB * B[i]        - cosTh *  iA2 * A[i];
    gradJ[i] = -iAB * (A[i]+B[i]) + cosTh * (iA2 * A[i] + iB2 * B[i]);
    gradK[i] =  iAB * A[i]        - cosTh               * iB2 * B[i];
  }
  return cosTh;
}

inline double torsion_gradient(double *F, double *G, double *H,
			       double *gradI, double *gradJ, double *gradK, double *gradL) {
  double A[3], B[3];
  double nA, nB;
  double cosPh;
  double AxF[3], BxF[3], AxG[3], BxG[3], AxH[3], BxH[3];
  int i;
  double iAB, iA2, iB2;
  // A = F x G
  cross_product(F,G, A);
  // B = H x G
  cross_product(H,G, B);
  nA = norm(A);
  nB = norm(B);
  cosPh = dot_product(A,B)/(nA*nB);
  // gradients
  cross_product(A,F, AxF);
  cross_product(B,F, BxF);
  cross_product(A,G, AxG);
  cross_product(B,G, BxG);
  cross_product(A,H, AxH);
  cross_product(B,H, BxH);
  iAB = 1.0/(nA*nB);
  iA2 = 1.0/(nA*nA);
  iB2 = 1.0/(nB*nB);
  
  for(i=0;i<3;i++) {
    gradI[i] = - iAB * BxG[i]                      + cosPh * iA2 * AxG[i];
    gradJ[i] =   iAB * (-AxH[i] + BxG[i] - BxF[i]) - cosPh * (iA2 * (AxG[i]-AxF[i]) - iB2*BxH[i]);
    gradK[i] =   iAB * ( AxG[i] + AxH[i] + BxF[i]) - cosPh * (iA2 * AxF[i] + iB2*(BxG[i] + BxH[i]));
    gradL[i] =  -iAB * AxG[i]                      + cosPh * iB2 * BxG[i];
  }
  return cosPh;
}

inline double inversion_gradient(double *U, double *V, double *B,
				 double *gradI, double *gradJ, double *gradK, double *gradL) {
  double A[3];
  double nA, nB;
  double sinPs;
  double AxV[3], BxV[3], AxU[3], BxU[3];
  int i;
  double iAB, iA2, iB2;
  // A = U x V
  cross_product(U,V, A);
  nA = norm(A);
  nB = norm(B);
  sinPs = dot_product(A,B)/(nA*nB);
  // gradients
  cross_product(A,V, AxV);
  cross_product(B,V, BxV);
  cross_product(A,U, AxU);
  cross_product(B,U, BxU);
  iAB = 1.0/(nA*nB);
  iA2 = 1.0/(nA*nA);
  iB2 = 1.0/(nB*nB);
  
  for(i=0;i<3;i++) {
    gradI[i] =  iAB * (-A[i] + BxV[i] - BxU[i]) - sinPs * (iA2*(AxV[i] - AxU[i]) - iB2 * B[i]);
    gradJ[i] = -iAB * BxV[i]                    + sinPs * iA2 * AxV[i];
    gradK[i] =  iAB * BxU[i]                    - sinPs * iA2 * AxU[i];
    gradL[i] =  iAB * A[i]                      - sinPs                          * iB2 * B[i];
  }
  return sinPs;
}

void set_current_state(FORCE_FIELD *ff, int state) {
  if (state > ff->nh) {
    printf("WARNING: Attempting to set state=%d, but there are only %d electronic states!\n", state, ff->nh);
    ff->state = ff->nh;
  } else {
    ff->state = state;
  }
}

void evaluate_force_field(FORCE_FIELD *ff) {
  int i;
  int I,J,K,L;
  NODE *node;
  BOND *bond;
  ANGLE *angle;
  TORSION *torsion;
  INVERSION *inversion;
  VAN_DER_WAALS *vdW;
  COULOMB *coul;
  HBOND *hbond;
  EXCITON *exciton;
  
  ATOM *atI, *atJ, *atK, *atL;
  REPLICA_ATOM *atIr, *atJr, *atKr, *atLr;
  double gradI[3], gradJ[3], gradK[3], gradL[3];
  double rJI[3], rJK[3], rKL[3];
  // for bond and vdW interactions
  double xIJ, yIJ, zIJ, nIJ, dR;
  double rinv, ri2, ri6, ri8, ri12, ri14, fac;
  // for angle interactions
  double cosTh;
  double dcos;
  // for torsion interactions
  double cosPh, x, x2, x3, x4, x5, x6, cosn, dcosn;
  // for inversion interactions
  double sinPs, cosPs, dcosPs;
  double rIJ[3], rIK[3], rIL[3];
  // for Coulomb interactions
  double ri3;
  // for H-bonds
  double nIK;
  double x10, x12, f, fprime;
  double cosTh3, cosTh4;
  double fac1, fac2;
  // for exciton couplings
  int elem1, elem2;
  // for Lapack
  int info;

  // evaluate bonded interaction
  for(node=ff->bond_list->head; node!=NULL; node=node->next) {
    bond = (BOND *) node->element;
    I = bond->I; J = bond->J;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    // distance vector
    xIJ = atJ->x - atI->x + atIr->tx - atJr->tx;
    yIJ = atJ->y - atI->y + atIr->ty - atJr->ty;
    zIJ = atJ->z - atI->z + atIr->tz - atJr->tz;
    nIJ = sqrt(xIJ*xIJ + yIJ*yIJ + zIJ*zIJ);

    dR = nIJ - bond->Re; 
    // energy
    bond->energy = 0.5*bond->Ke*dR*dR;
    // gradient
    fac = -bond->Ke*dR/nIJ;
    bond->gradI[0] = fac * xIJ;
    bond->gradI[1] = fac * yIJ;
    bond->gradI[2] = fac * zIJ;

    bond->r = nIJ;
    // derivative of bond length w/r/t to atom I
    bond->derI[0] = -xIJ/nIJ;
    bond->derI[1] = -yIJ/nIJ;
    bond->derI[2] = -zIJ/nIJ;
  }
  //printf("bonds evaluated\n");
  // angles
  for(node=ff->angle_list->head; node!=NULL; node=node->next) {
    angle = (ANGLE *) node->element;
    angle->energy = 0.0;
    I = angle->I; J = angle->J; K = angle->K;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    atKr = &(ff->atoms_repl[K]);
    atK = &(ff->atoms[atKr->index]);

    // vector J-->I
    rJI[0] = atI->x - atJ->x + atJr->tx - atIr->tx;
    rJI[1] = atI->y - atJ->y + atJr->ty - atIr->ty;
    rJI[2] = atI->z - atJ->z + atJr->tz - atIr->tz;
    // vector J-->K
    rJK[0] = atK->x - atJ->x + atJr->tx - atKr->tx;
    rJK[1] = atK->y - atJ->y + atJr->ty - atKr->ty;
    rJK[2] = atK->z - atJ->z + atJr->tz - atKr->tz;

    cosTh = angle_gradient(rJI, rJK, gradI, gradJ, gradK);

    angle->cosTh = cosTh;
    // derivative of cos(theta_IJK) w/r/t to atom positions
    for(i=0;i<3;i++) {
      angle->derI[i] = gradI[i];
      angle->derJ[i] = gradJ[i];
      angle->derK[i] = gradK[i];
    }
    
    if (angle->linear == 1) {
      // E = Ce * (1 + cos(theta_IJK))
      fac = angle->Ce;
      angle->energy = fac * (1.0 + cosTh);
      for(i=0;i<3;i++) {
	angle->gradI[i] = fac * gradI[i];
	angle->gradJ[i] = fac * gradJ[i];
	angle->gradK[i] = fac * gradK[i];
      }
    } else {
      // E = 1/2 * Ce * [cos(theta_IJK) - cos(Te)]^2
      dcos = cosTh - angle->cosTe;
      angle->energy = 0.5 * angle->Ce * dcos*dcos;
      fac = angle->Ce * dcos;
      for(i=0;i<3;i++) {
	angle->gradI[i] = fac * gradI[i];
	angle->gradJ[i] = fac * gradJ[i];
	angle->gradK[i] = fac * gradK[i];
      }      
    }
  }
  //printf("angles evaluated\n");
  // torsions
  for(node=ff->torsion_list->head; node!=NULL; node=node->next) {
    torsion = (TORSION *) node->element;
    torsion->energy = 0.0;
    I = torsion->I; J = torsion->J; K = torsion->K; L = torsion->L;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    atKr = &(ff->atoms_repl[K]);
    atK = &(ff->atoms[atKr->index]);
    atLr = &(ff->atoms_repl[L]);
    atL = &(ff->atoms[atLr->index]);
    // vector J-->I
    rJI[0] = atI->x - atJ->x + atJr->tx - atIr->tx;
    rJI[1] = atI->y - atJ->y + atJr->ty - atIr->ty;
    rJI[2] = atI->z - atJ->z + atJr->tz - atIr->tz;
    // vector J-->K
    rJK[0] = atK->x - atJ->x + atJr->tx - atKr->tx;
    rJK[1] = atK->y - atJ->y + atJr->ty - atKr->ty;
    rJK[2] = atK->z - atJ->z + atJr->tz - atKr->tz;
    // vector K-->L
    rKL[0] = atL->x - atK->x + atKr->tx - atLr->tx;
    rKL[1] = atL->y - atK->y + atKr->ty - atLr->ty;
    rKL[2] = atL->z - atK->z + atKr->tz - atLr->tz;

    cosPh = torsion_gradient(rJI, rJK, rKL,
			     gradI, gradJ, gradK, gradL);

    torsion->cosPh = cosPh;
    // derivative of cos(phi_IJKL) w/r/t to atom positions
    for(i=0;i<3;i++) {
      torsion->derI[i] = gradI[i];
      torsion->derJ[i] = gradJ[i];
      torsion->derK[i] = gradK[i];
      torsion->derL[i] = gradL[i];
    }

    // use Chebyshev polynomials to evaluate cos(n*phi)
    // cos(n*phi) = sum_m T_n,m * cos(phi)^m
    //
    x = cosPh;
    if (torsion->n == 2) {
      x2 = x*x;
      cosn = 2.0*x2-1.0;
      // derivative
      dcosn = 4.0*x;
    } else if (torsion->n == 3) {
      x2 = x*x;
      x3 = x2*x;
      cosn = 4.0*x3 - 3.0*x;
      dcosn = 12.0*x2 - 3.0;
    } else if (torsion->n == 6) {
      x2 = x*x;
      x3 = x2*x;
      x4 = x2*x2;
      x5 = x3*x2;
      x6 = x3*x3;
      cosn = 32.0*x6 - 48.0*x4 + 18.0*x2 - 1.0;
      dcosn = 192.0*x5 - 192.0*x3 + 36.0*x;
    } else {
      cosn = 0.0;
      dcosn = 0.0;
      printf("WARNING: periodicity of torsion should be n=2,3 or 6!");
    }
    // cos(n*(phi-phi0)) = cos(n*phi)*cos(n*phi0)   provided that phi0 = k*pi/n
    // E = V * [1 - cos(n*(phi-phi0))]
    torsion->energy = torsion->Vbarrier * (1.0 - torsion->cosPhi0 * cosn);
    // compute gradients of torsion energy
    fac = -torsion->Vbarrier * torsion->cosPhi0 * dcosn;
    for(i=0;i<3;i++) {
      torsion->gradI[i] = fac * gradI[i];
      torsion->gradJ[i] = fac * gradJ[i];
      torsion->gradK[i] = fac * gradK[i];
      torsion->gradL[i] = fac * gradL[i];
    }
  }
  //printf("torsions evaluated\n");
  // inversions
  for(node=ff->inversion_list->head; node!=NULL; node=node->next) {
    inversion = (INVERSION *) node->element;
    inversion->energy = 0.0;
    I = inversion->I; J = inversion->J; K = inversion->K; L = inversion->L;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    atKr = &(ff->atoms_repl[K]);
    atK = &(ff->atoms[atKr->index]);
    atLr = &(ff->atoms_repl[L]);
    atL = &(ff->atoms[atLr->index]);
    // vector I-->J
    rIJ[0] = atJ->x - atI->x + atIr->tx - atJr->tx;
    rIJ[1] = atJ->y - atI->y + atIr->ty - atJr->ty;
    rIJ[2] = atJ->z - atI->z + atIr->tz - atJr->tz;
    // vector I-->K
    rIK[0] = atK->x - atI->x + atIr->tx - atKr->tx;
    rIK[1] = atK->y - atI->y + atIr->ty - atKr->ty;
    rIK[2] = atK->z - atI->z + atIr->tz - atKr->tz;
    // vector I-->L
    rIL[0] = atL->x - atI->x + atIr->tx - atLr->tx;
    rIL[1] = atL->y - atI->y + atIr->ty - atLr->ty;
    rIL[2] = atL->z - atI->z + atIr->tz - atLr->tz;

    // angle between I-->L  and the plane I-J-K
    sinPs = inversion_gradient(rIJ, rIK,  rIL,
			       gradI, gradJ, gradK, gradL);
    cosPs = sqrt(1.0 - sinPs*sinPs);
    dcosPs = cosPs - inversion->cosPs0;
    inversion->energy += 1.0/6.0 * inversion->Kinv * dcosPs*dcosPs;
    // gradient
    fac = -1.0/3.0 * inversion->Kinv * dcosPs * sinPs/cosPs;
    for(i=0;i<3;i++) {
      inversion->gradI[i] = fac * gradI[i];
      inversion->gradJ[i] = fac * gradJ[i];
      inversion->gradK[i] = fac * gradK[i];
      inversion->gradL[i] = fac * gradL[i];
    }
    
    // derivatives for Wilson's B-matrix, only the angle between
    // bond I-->L and the plane I-J-K is used
    inversion->sinPs = sinPs;
    for(i=0;i<3;i++) {
      inversion->derI[i] = gradI[i];
      inversion->derJ[i] = gradJ[i];
      inversion->derK[i] = gradK[i];
      inversion->derL[i] = gradL[i];
    }
    
    // angle between I-->J  and the plane I-K-L
    sinPs = inversion_gradient(rIK, rIL,  rIJ,
			       gradI, gradK, gradL, gradJ);
    cosPs = sqrt(1.0 - sinPs*sinPs);
    dcosPs = cosPs - inversion->cosPs0;
    inversion->energy += 1.0/6.0 * inversion->Kinv * dcosPs*dcosPs;
    // gradient
    fac = -1.0/3.0 * inversion->Kinv * dcosPs * sinPs/cosPs;
    for(i=0;i<3;i++) {
      inversion->gradI[i] += fac * gradI[i];
      inversion->gradJ[i] += fac * gradJ[i];
      inversion->gradK[i] += fac * gradK[i];
      inversion->gradL[i] += fac * gradL[i];
    }
    // angle between I-->K  and the plane I-L-J
    sinPs = inversion_gradient(rIL, rIJ,  rIK,
			       gradI, gradL, gradJ, gradK);
    cosPs = sqrt(1.0 - sinPs*sinPs);
    dcosPs = cosPs - inversion->cosPs0;
    inversion->energy += 1.0/6.0 * inversion->Kinv * dcosPs*dcosPs;
    // gradient
    fac = -1.0/3.0 * inversion->Kinv * dcosPs * sinPs/cosPs;
    for(i=0;i<3;i++) {
      inversion->gradI[i] += fac * gradI[i];
      inversion->gradJ[i] += fac * gradJ[i];
      inversion->gradK[i] += fac * gradK[i];
      inversion->gradL[i] += fac * gradL[i];
    }
    
    
  }
  //printf("inversions evaluated\n");  
  // evaluate non-bonding interaction
  for(node=ff->nonbonded_list->head; node!=NULL; node=node->next) {
    vdW = (VAN_DER_WAALS *) node->element;
    I = vdW->I; J = vdW->J;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    // distance vector
    xIJ = atJ->x - atI->x + atIr->tx - atJr->tx;
    yIJ = atJ->y - atI->y + atIr->ty - atJr->ty;
    zIJ = atJ->z - atI->z + atIr->tz - atJr->tz;
    nIJ = sqrt(xIJ*xIJ + yIJ*yIJ + zIJ*zIJ);

    rinv = 1.0/nIJ;
    ri2 = rinv*rinv;
    ri6 = pow(rinv,6);
    ri8 = ri6*ri2;
    ri12 = ri6*ri6;
    ri14 = ri8*ri6;
    // energy
    vdW->energy = vdW->A * ri12 - vdW->B*ri6;
    // gradient
    fac = 12.0 * vdW->A * ri14 - 6.0 * vdW->B * ri8;
    vdW->gradI[0] = fac * xIJ;
    vdW->gradI[1] = fac * yIJ;
    vdW->gradI[2] = fac * zIJ;
  }
  //printf("vdWs evaluated\n");  
  // evaluate Coulomb interaction
  for(node=ff->coulomb_list->head; node!=NULL; node=node->next) {
    coul = (COULOMB *) node->element;
    I = coul->I; J = coul->J;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    // distance vector
    xIJ = atJ->x - atI->x + atIr->tx - atJr->tx;
    yIJ = atJ->y - atI->y + atIr->ty - atJr->ty;
    zIJ = atJ->z - atI->z + atIr->tz - atJr->tz;
    nIJ = sqrt(xIJ*xIJ + yIJ*yIJ + zIJ*zIJ);

    rinv = 1.0/nIJ;
    ri3 = rinv*rinv*rinv;
    // energy
    coul->energy = coul->qq * rinv;
    // gradient
    fac = coul->qq * ri3;
    coul->gradI[0] = fac * xIJ;
    coul->gradI[1] = fac * yIJ;
    coul->gradI[2] = fac * zIJ;
  }
  //printf("Coulomb interactions evaluated\n");
  // evaluate hydrogen bonds
  for(node=ff->hbond_list->head; node!=NULL; node=node->next) {
    hbond = (HBOND *) node->element;
    // donor         hydrogen    acceptor
    I = hbond->I; J = hbond->J; K = hbond->K;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    atKr = &(ff->atoms_repl[K]);
    atK = &(ff->atoms[atKr->index]);

    // vector J-->I
    rJI[0] = atI->x - atJ->x + atJr->tx - atIr->tx;
    rJI[1] = atI->y - atJ->y + atJr->ty - atIr->ty;
    rJI[2] = atI->z - atJ->z + atJr->tz - atIr->tz;
    // vector J-->K
    rJK[0] = atK->x - atJ->x + atJr->tx - atKr->tx;
    rJK[1] = atK->y - atJ->y + atJr->ty - atKr->ty;
    rJK[2] = atK->z - atJ->z + atJr->tz - atKr->tz;

    cosTh = angle_gradient(rJI, rJK, gradI, gradJ, gradK);

    // distance vector I-->K  between donor and acceptor
    rIK[0] = atK->x - atI->x + atIr->tx - atKr->tx;
    rIK[1] = atK->y - atI->y + atIr->ty - atKr->ty;
    rIK[2] = atK->z - atI->z + atIr->tz - atKr->tz;
    nIK = sqrt(rIK[0]*rIK[0] + rIK[1]*rIK[1] + rIK[2]*rIK[2]);
    rinv = 1.0/nIK;
    x = hbond->Rhb * rinv;
    x2 = x*x;
    x10 = pow(x,10);
    x12 = x10*x2;

    f = hbond->Dhb * (5.0*x12 - 6.0*x10);
    fprime = - hbond->Dhb * 60.0*(x12 - x10) * rinv;

    cosTh3 = pow(cosTh, 3);
    cosTh4 = cosTh*cosTh3;
    
    hbond->energy = f * cosTh4;

    fac1 = fprime*cosTh4*rinv;
    fac2 = 4.0*f*cosTh3;

    for(i=0;i<3;i++) {
      hbond->gradI[i] = -fac1 * rIK[i] + fac2 * gradI[i];
      hbond->gradJ[i] =                  fac2 * gradJ[i];
      hbond->gradK[i] =  fac1 * rIK[i] + fac2 * gradK[i];
    }
    //printf("energy HBOND(%d-%d--%d)=%f\n", I+1, J+1, K+1, hbond->energy);
  }
  //printf("H-bonds evaluated\n");

  // evaluate excitonic couplings
  if (ff->nchromo > 0) {
    // zero all elements of Hmatrix
    for(elem1=0; elem1<ff->nh; elem1++) {
      for(elem2=0; elem2<ff->nh; elem2++) {
	ff->Hmatrix[elem1*ff->nh + elem2] = 0.0;
      }
      // diagonal elements are just the monomer excitation energies
      ff->Hmatrix[elem1*ff->nh + elem1] = ff->Hmatrix_diag[elem1];
    }
    // compute couplings and fill in Hamiltonian matrix
    for(node=ff->exciton_list->head; node!=NULL; node=node->next) {
      exciton = (EXCITON *) node->element;
      I = exciton->I; J = exciton->J;
      // same as Coulomb coupling
      atIr = &(ff->atoms_repl[I]);
      atI = &(ff->atoms[atIr->index]);
      atJr = &(ff->atoms_repl[J]);
      atJ = &(ff->atoms[atJr->index]);
      // distance vector
      xIJ = atJ->x - atI->x + atIr->tx - atJr->tx;
      yIJ = atJ->y - atI->y + atIr->ty - atJr->ty;
      zIJ = atJ->z - atI->z + atIr->tz - atJr->tz;
      nIJ = sqrt(xIJ*xIJ + yIJ*yIJ + zIJ*zIJ);
      
      rinv = 1.0/nIJ;
      ri3 = rinv*rinv*rinv;
      // energy
      exciton->energy = exciton->qq * rinv;
      // gradient
      fac = exciton->qq * ri3;
      exciton->gradI[0] = fac * xIJ;
      exciton->gradI[1] = fac * yIJ;
      exciton->gradI[2] = fac * zIJ;
      
      // add coupling to the respective matrix element
      ff->Hmatrix[exciton->elem] += exciton->energy;
    }
    /*
    // DEBUG
    // print Hamiltonian
    printf("excitonic Hamiltonian (only lower triangle of symmetric matrix is filled)\n");
    for(I=0; I<ff->nh; I++) {
      for(J=0; J<ff->nh; J++) {
	printf("%+4.3f  ", ff->Hmatrix[I*ff->nh+J]);
      }
      printf("\n");
    }
    */
    // diagonalize excitonic Hamiltonian
    // In row major order the lower half of H is set, but Fortran uses column major
    // order, so we have to use "Upper"
    dsyev_( "Vectors", "Upper", &ff->nh, ff->Hmatrix, &ff->nh, ff->exciton_energies,
	   ff->work, &ff->lwork, &info);
    /*
    // DEBUG
    // show eigenvalues
    for(I=0; I<ff->nh; I++) {
      printf("E(%d) = %f\n", I, ff->exciton_energies[I]);
    }
    */
    // The eigenvectors are stored columnwise in Hmatrix
    if (info > 0) {
      printf("ERROR: DSYEV failed to diagonalize exciton Hamiltonian!\n");
      exit(info);
    }
  //printf("exciton couplings evaluated\n");
  }
}

double collect_force_field(FORCE_FIELD *ff)
{
  int i,n;
  int I,J,K,L;
  NODE *node;
  BOND *bond;
  ANGLE *angle;
  TORSION *torsion;
  INVERSION *inversion;
  VAN_DER_WAALS *vdW;
  COULOMB *coul;
  HBOND *hbond;
  EXCITON *exciton;
  int st; // current exciton state
  double ci,cj, cc;
  
  ATOM *atI, *atJ, *atK, *atL;
  REPLICA_ATOM *atIr, *atJr, *atKr, *atLr;

  double energy = 0.0;

  // for testing
  double nrm2;
  double energy_expec = 0.0;
  
  // set all gradients to 0
  for(n=0;n<ff->nat;n++) {
    atI = &(ff->atoms[n]);
    for(i=0;i<3;i++) {
      atI->grad[i] = 0.0;
    }
  }

  // bonds
  for(node=ff->bond_list->head; node!=NULL; node=node->next) {
    bond = (BOND *) node->element;
    I = bond->I; J = bond->J;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    // add energies
    energy += bond->energy;
    // add gradients
    for(i=0;i<3;i++) {
      atI->grad[i] += bond->gradI[i];
      atJ->grad[i] -= bond->gradI[i];
    }
  }
  //  angles
  for(node=ff->angle_list->head; node!=NULL; node=node->next) {
    angle = (ANGLE *) node->element;
    I = angle->I; J = angle->J; K = angle->K;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    atKr = &(ff->atoms_repl[K]);
    atK = &(ff->atoms[atKr->index]);
    // add energies
    energy += angle->energy;
    // add gradients
    for(i=0;i<3;i++) {
      atI->grad[i] += angle->gradI[i];
      atJ->grad[i] += angle->gradJ[i];
      atK->grad[i] += angle->gradK[i];
    }
  }
  // torsions
  for(node=ff->torsion_list->head; node!=NULL; node=node->next) {
    torsion = (TORSION *) node->element;
    I = torsion->I; J = torsion->J; K = torsion->K; L = torsion->L;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    atKr = &(ff->atoms_repl[K]);
    atK = &(ff->atoms[atKr->index]);
    atLr = &(ff->atoms_repl[L]);
    atL = &(ff->atoms[atLr->index]);
    // add energies
    energy += torsion->energy;
    // add gradients
    for(i=0;i<3;i++) {
      atI->grad[i] += torsion->gradI[i];
      atJ->grad[i] += torsion->gradJ[i];
      atK->grad[i] += torsion->gradK[i];
      atL->grad[i] += torsion->gradL[i];
    }
  }
  // inversions
  for(node=ff->inversion_list->head; node!=NULL; node=node->next) {
    inversion = (INVERSION *) node->element;
    I = inversion->I; J = inversion->J; K = inversion->K; L = inversion->L;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    atKr = &(ff->atoms_repl[K]);
    atK = &(ff->atoms[atKr->index]);
    atLr = &(ff->atoms_repl[L]);
    atL = &(ff->atoms[atLr->index]);
    // add energies
    energy += inversion->energy;
    // add gradients
    for(i=0;i<3;i++) {
      atI->grad[i] += inversion->gradI[i];
      atJ->grad[i] += inversion->gradJ[i];
      atK->grad[i] += inversion->gradK[i];
      atL->grad[i] += inversion->gradL[i];
    }
  }
  // non-bonding interactions
  for(node=ff->nonbonded_list->head; node!=NULL; node=node->next) {
    vdW = (VAN_DER_WAALS *) node->element;
    I = vdW->I; J = vdW->J;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    // add energies
    energy += vdW->energy;
    // add gradients
    for(i=0;i<3;i++) {
      atI->grad[i] += vdW->gradI[i];
      atJ->grad[i] -= vdW->gradI[i];
    }
  }
  // Coulomb interactions
  for(node=ff->coulomb_list->head; node!=NULL; node=node->next) {
    coul = (COULOMB *) node->element;
    I = coul->I; J = coul->J;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    // add energies
    energy += coul->energy;
    // add gradients
    for(i=0;i<3;i++) {
      atI->grad[i] += coul->gradI[i];
      atJ->grad[i] -= coul->gradI[i];
    }
  }
  // hydrogen bonds
  for(node=ff->hbond_list->head; node!=NULL; node=node->next) {
    hbond = (HBOND *) node->element;
    I = hbond->I; J = hbond->J; K = hbond->K;
    atIr = &(ff->atoms_repl[I]);
    atI = &(ff->atoms[atIr->index]);
    atJr = &(ff->atoms_repl[J]);
    atJ = &(ff->atoms[atJr->index]);
    atKr = &(ff->atoms_repl[K]);
    atK = &(ff->atoms[atKr->index]);
    // add energies
    energy += hbond->energy;
    // add gradients
    for(i=0;i<3;i++) {
      atI->grad[i] += hbond->gradI[i];
      atJ->grad[i] += hbond->gradJ[i];
      atK->grad[i] += hbond->gradK[i];
    }
  }

  // excitonic coupling
  if (ff->state > 0) {    
    st = ff->state-1;
    // DEBUB
    // check normalization of eigenfunction and compute expectation value of energy
    nrm2 = 0.0;
    energy_expec = 0.0;
    for(i=0; i<ff->nh; i++) {
      ci = ff->Hmatrix[i*ff->nh + st];
      nrm2 += ci*ci;
      energy_expec += ci*ci * ff->Hmatrix_diag[i];
    }
    //printf("norm of eigenvector = %f\n", nrm2);
    //

    // add energies
    //printf("exciton energy: %f\n", ff->exciton_energies[st]);
    energy += ff->exciton_energies[st];
    // compute Hellmann-Feynmann gradients
    // dE(I)/dR = sum_i sum_j C(i,I) dH(i,j)/dR C(j,J)
    for(node=ff->exciton_list->head; node!=NULL; node=node->next) {
      exciton = (EXCITON *) node->element;
      I = exciton->I; J = exciton->J;
      atIr = &(ff->atoms_repl[I]);
      atI = &(ff->atoms[atIr->index]);
      atJr = &(ff->atoms_repl[J]);
      atJ = &(ff->atoms[atJr->index]);
      
      // coefficients of eigenfunctions are stored columnwise in H
      ci = ff->Hmatrix[exciton->elem1*ff->nh+st];
      cj = ff->Hmatrix[exciton->elem2*ff->nh+st];
      
      cc = 2.0*ci*cj;

      // add gradients
      for(i=0;i<3;i++) {
	atI->grad[i] += cc * exciton->gradI[i];
	atJ->grad[i] -= cc * exciton->gradI[i];
      }

      // DEBUG
      energy_expec += cc * exciton->energy;
    }
    // compare excitation energy with expectation value E(I) = <I|H|I>
    //printf("exciton energy E(I)= %f    <I|H|I>= %f\n", ff->exciton_energies[st], energy_expec);
  }
  return energy;
}

void delete_force_field(FORCE_FIELD *ff)
{
  int n;
  CHROMOPHORE *chromo;
  
  free(ff->atoms);
  free(ff->atoms_repl);
  list_delete(ff->bond_list);
  list_delete(ff->angle_list);
  list_delete(ff->torsion_list);
  list_delete(ff->inversion_list);
  list_delete(ff->nonbonded_list);
  list_delete(ff->coulomb_list);
  list_delete(ff->hbond_list);

  // free chromophores
  for(n=0; n<ff->nchromo; n++) {
      chromo = &(ff->chromophores[n]);
      free(chromo->indeces);
      free(chromo->excitation_energies);
      free(chromo->transition_charges);
    }
  if (ff->chromophores != NULL) {
    free(ff->chromophores);
  }
  // free exciton Hamiltonian
  if (ff->Hmatrix != NULL) {
    free(ff->Hmatrix);
    free(ff->Hmatrix_diag);
    free(ff->exciton_energies);
    free(ff->elec_tdip);
    free(ff->magn_tdip);
    free(ff->work);
  }
  list_delete(ff->exciton_list);
  free(ff);
}

void local_chromophore_axes(FORCE_FIELD *ff)
{
  /*
    compute the local axes for each chromophore by diagonalizing the tensor of
    inertia (with all masses set to 1 !)
    The magnetic transition dipoles moments are defined relative to the local axes,
    so that they rotate with the molecule.
   */
  int n;
  CHROMOPHORE *chromo;
  int i, I;
  ATOM *atI;
  double x,y,z;
  double cm[3];  // center of mass
  double Ixx, Iyy, Izz, Ixy, Iyz, Ixz;
  double moments[3]; // principle moments of inertia, eigenvalues of tensor I
  int info;
  double work[100];
  int lwork=100;
  int lda = 3;
  double pos[3], proj[3];
  int j, k;
  double xaxis[3], yaxis[3], zaxis[3];
  
  for(n=0; n<ff->nchromo; n++) {  // enumerate chromophores
    chromo = &(ff->chromophores[n]);
    // compute center of mass (all masses set to 1)
    cm[0] = 0.0; cm[1] = 0.0; cm[2] = 0.0;
    for(i=0; i<chromo->nat; i++) { // enumerate atoms belonging to chromophore i
      I = chromo->indeces[i];
      atI = &(ff->atoms[I]);
      x = atI->x; y = atI->y; z = atI->z;
      cm[0] += x;
      cm[1] += y;
      cm[2] += z;
    }
    for(k=0; k<3; k++) {
      cm[k] /= (float) chromo->nat;
    }
    
    // compute tensor of inertia
    Ixx = Iyy = Izz = Ixy = Iyz = Ixz = 0.0;
    for(i=0; i<chromo->nat; i++) { // enumerate atoms belonging to chromophore i
      I = chromo->indeces[i];
      atI = &(ff->atoms[I]);
      // shift chromophore to center of mass
      x = atI->x - cm[0]; y = atI->y - cm[1]; z = atI->z - cm[2];
      Ixx += y*y+z*z;
      Iyy += x*x+z*z;
      Izz += x*x+y*y;
      Ixy += -x*y;
      Iyz += -y*z;
      Ixz += -x*z;
    }
    
    // diagonalize symmetric tensor
    chromo->axes[0][0] = Ixx;
    chromo->axes[1][1] = Iyy;
    chromo->axes[2][2] = Izz;
    
    chromo->axes[1][0] = Ixy;
    chromo->axes[2][0] = Ixz;
    chromo->axes[2][1] = Iyz;

    // In row major order the lower half of I is set, but Fortran uses column major
    // order, so we have to use "Upper"
    dsyev_("Vectors", "Upper", &lda, (double *) chromo->axes, &lda, moments, work, &lwork, &info);

    // Phases of eigenvectors are not unique, fix them by some convention
    // find projection of atom positions on axes
    proj[0] = 0.0; proj[1] = 0.0; proj[2] = 0.0;
    for(i=0; i<chromo->nat; i++) { // enumerate atoms belonging to chromophore i
      I = chromo->indeces[i];
      atI = &(ff->atoms[I]);
      pos[0] = atI->x - cm[0]; pos[1] = atI->y - cm[1]; pos[2] = atI->z - cm[2];
      for(j=0; j<3; j++) {
	for(k=0; k<3; k++) {
	  proj[j] += (i+1) * pos[k] * chromo->axes[j][k];
	}
      }
    }
    printf("proj = %+7.5f %7.5f %7.5f\n", proj[0], proj[1], proj[2]);
    // find right-handed coordinate system
    if (proj[2] < 0.0) {
      printf("flip z-axis\n");
      for(k=0; k<3; k++) {
	chromo->axes[2][k] *= -1;
      }
    }
    if (proj[1] < 0.0) {
      printf("flip y-axis\n");
      for(k=0; k<3; k++) {
	chromo->axes[1][k] *= -1;
      }
    }
    // triple product
    for(k=0; k<3; k++) {
      xaxis[k] = chromo->axes[0][k];
      yaxis[k] = chromo->axes[1][k];
      zaxis[k] = chromo->axes[2][k];
    }
    if (triple_product(zaxis, xaxis, yaxis) < 0.0) {
      printf("flip x-axis\n");
      for(k=0; k<3; k++) {
	chromo->axes[0][k] *= -1;
      }
    }
    
    printf("Local coordinates, chromophore %d\n", n);
    printf("  center of mass (all masses = 1): %+7.5f %+7.5f %+7.5f\n", cm[0], cm[1], cm[2]);
    printf("  tensor of inertia (all masses = 1)\n");
    printf("    %+7.5f \n", Ixx);
    printf("    %+7.5f %+7.5f \n", Ixy, Iyy);
    printf("    %+7.5f %+7.5f %+7.5f \n", Ixz, Iyz, Izz);
    printf("  principal moments: %+7.5f %+7.5f %+7.5f\n", moments[0], moments[1], moments[2]);
    printf("  x-axis: %+7.5f %+7.5f %+7.5f \n", chromo->axes[0][0], chromo->axes[0][1], chromo->axes[0][2]);
    printf("  y-axis: %+7.5f %+7.5f %+7.5f \n", chromo->axes[1][0], chromo->axes[1][1], chromo->axes[1][2]);
    printf("  z-axis: %+7.5f %+7.5f %+7.5f \n", chromo->axes[2][0], chromo->axes[2][1], chromo->axes[2][2]);
    printf("\n");
    
  }
}

void evaluate_transition_dipoles(FORCE_FIELD *ff)
{
  /*
    compute electronic and magnetic transition dipoles between the ground
    state and the excitonic basis states, where 1 chromophore is excited and all others
    are in their respective ground states.
   */
  int n, c, s; // enumerate cells, chromophores and states
  int b; // enumerate excitonic basis states
  int i, I; // enumerate atoms 
  int xyz;
  CHROMOPHORE *chromo;
  ATOM *atI;
  REPLICA_ATOM *atIr;
  double q; // electric transition charge
  double mloc[3], m[3]; // magnetic transition dipoles in local and global frame
  double xI, yI, zI;
  
  // need local axes for orientation of magnetic dipoles
  local_chromophore_axes(ff);

  b = 0;
  for(n=0; n<ff->ncells; n++) {  
    for(c=0; c<ff->nchromo; c++) {
      chromo = &(ff->chromophores[c]);
      for(s=0; s<chromo->nst; s++) {
	// electric transition dipole
	ff->elec_tdip[3*b+0] = 0.0;
	ff->elec_tdip[3*b+1] = 0.0;
	ff->elec_tdip[3*b+2] = 0.0;
	// magnetic transition dipole
	ff->magn_tdip[3*b+0] = 0.0;
	ff->magn_tdip[3*b+1] = 0.0;
	ff->magn_tdip[3*b+2] = 0.0;
	for(i=0; i<chromo->nat; i++) {	  
	  I = n*ff->nat + chromo->indeces[i];
	  atIr = &(ff->atoms_repl[I]);
	  atI = &(ff->atoms[atIr->index]);
	  // position vector
	  xI = atI->x - atIr->tx;
	  yI = atI->y - atIr->ty;
	  zI = atI->z - atIr->tz;
	  // transition charge
	  q = chromo->transition_charges[4*i*chromo->nst+4*s+0];
	  
	  ff->elec_tdip[3*b+0] += q*xI;
	  ff->elec_tdip[3*b+1] += q*yI;
	  ff->elec_tdip[3*b+2] += q*zI;

	  // magnetic dipoles in local frame
	  mloc[0] = chromo->transition_charges[4*i*chromo->nst+4*s+1];
	  mloc[1] = chromo->transition_charges[4*i*chromo->nst+4*s+2];
	  mloc[2] = chromo->transition_charges[4*i*chromo->nst+4*s+3];

	  // transform magnetic dipoles from local molecular frame into global frame
	  for(xyz=0; xyz<3; xyz++) {
	    m[xyz] =  mloc[0] * chromo->axes[0][xyz]
	            + mloc[1] * chromo->axes[1][xyz]
	            + mloc[2] * chromo->axes[2][xyz];
	    // add to total magnetic dipole moment
	    ff->magn_tdip[3*b+xyz] += m[xyz];
	  }
	}
	b += 1;
      }
    }
  }
}

void set_force_field_coordinates(FORCE_FIELD *ff, int n, double x, double y, double z)
{
  ATOM *atom;
  atom = &(ff->atoms[n]);
  atom->x = x;
  atom->y = y;
  atom->z = z;
}

void get_force_field_gradient(FORCE_FIELD *ff, int n, double *gradx, double *grady, double *gradz)
{
  ATOM *atom;
  atom = &(ff->atoms[n]);
  *gradx = atom->grad[0];
  *grady = atom->grad[1];
  *gradz = atom->grad[2];
}

int get_number_of_atoms(FORCE_FIELD *ff)
{
  return ff->nat;
}

int get_number_of_excitons(FORCE_FIELD *ff)
{
  return ff->nh;
}

void get_exciton_energy(FORCE_FIELD *ff, double *energy, int state)
{
  *energy = ff->exciton_energies[state];
}

void get_exciton_coefficient(FORCE_FIELD *ff, double *coef, int i, int state)
{
  *coef = ff->Hmatrix[i*ff->nh + state];
}

void get_transition_dipoles(FORCE_FIELD *ff,
			    double *tx, double *ty, double *tz,
			    double *mx, double *my, double *mz, int b)
{
  /*
    fetch the electronic and magnetic transition dipoles for the transition from the ground
    state to the b-th excitonic basis state (index b starts from 0)
   */
  if (b >= ff->nh) {
    printf("WARNING: There are only %d excitonic basis states, got b=%d!\n", ff->nh, b);
    b = ff->nh-1;
  }
  
  *tx = ff->elec_tdip[3*b+0];
  *ty = ff->elec_tdip[3*b+1];
  *tz = ff->elec_tdip[3*b+2];

  *mx = ff->magn_tdip[3*b+0];
  *my = ff->magn_tdip[3*b+1];
  *mz = ff->magn_tdip[3*b+2];
}

int count_internal_coordinates(FORCE_FIELD *ff)
{
  // count the number of bonds, angles and torsions
  int n=0;

  n += (int) list_length(ff->bond_list);
  n += (int) list_length(ff->angle_list);
  n += (int) list_length(ff->torsion_list);
  n += (int) list_length(ff->inversion_list);
  
  return n;
}

void get_coordinate_definitions(FORCE_FIELD *ff,
				char **types_ptr, int **atoms_ptr, int *length)
{
  /*
  build list of characters which describe the type of each
  internal coordinate:
     'B' - bond
     'A' - valence angle
     'D' - dihedral angle
     'I' - inversion angle
  types[k] contains the character label for the k-th internal coordinate

  Each internal coordinate is defined by at most 4 atom indices
     I J      -  bond
     I J K    -  valence angle
     I J K L  -  dihedral angle
     I J K L  -  inversion, angle between bond I--L and plane I-J-K

  I = atoms[k*4], J = atoms[k*4+1], K = atoms[k*4+2], L = atoms[k*4+3]
  are the atom indices defining the k-th internal coordinate. For bonds
  and angles the indices K and L should not be used.
  */
  int I,J,K,L;
  NODE *node;
  BOND *bond;
  ANGLE *angle;
  TORSION *torsion;
  INVERSION *inversion;
  
  int nred;
  char *types;
  int *atoms;
  int k;
  
  // count number of redundant internal coordinates
  nred = count_internal_coordinates(ff);
  // allocate memory for labels and atom indices
  types = (char *) malloc(sizeof(char) * nred);
  atoms = (int *) malloc(sizeof(int) * 4 * nred);
  *length = nred;
  
  k=0;  // enumerates redundant internal coordinates
  // bonds
  for(node=ff->bond_list->head; node!=NULL; node=node->next) {
    bond = (BOND *) node->element;
    I = bond->I; J = bond->J;
    types[k] = 'B';
    atoms[4*k  ] = I;
    atoms[4*k+1] = J;
    atoms[4*k+2] = -1;
    atoms[4*k+3] = -1;
    k++;
  }
  // valence angles
  for(node=ff->angle_list->head; node!=NULL; node=node->next) {
    angle = (ANGLE *) node->element;
    I = angle->I; J = angle->J; K = angle->K;
    types[k] = 'A';
    atoms[4*k  ] = I;
    atoms[4*k+1] = J;
    atoms[4*k+2] = K;
    atoms[4*k+3] = -1;
    k++;
  }
  // torsions
  for(node=ff->torsion_list->head; node!=NULL; node=node->next) {
    torsion = (TORSION *) node->element;
    I = torsion->I; J = torsion->J; K = torsion->K; L = torsion->L;
    types[k] = 'D';
    atoms[4*k  ] = I;
    atoms[4*k+1] = J;
    atoms[4*k+2] = K;
    atoms[4*k+3] = L;
    k++;
  }
  // inversions
  for(node=ff->inversion_list->head; node!=NULL; node=node->next) {
    inversion = (INVERSION *) node->element;
    I = inversion->I; J = inversion->J; K = inversion->K; L = inversion->L;
    types[k] = 'I';
    atoms[4*k  ] = I;
    atoms[4*k+1] = J;
    atoms[4*k+2] = K;
    atoms[4*k+3] = L;
    k++;
  }

  if (k != nred) {
    printf("BUG: k= %d  !=  nred= %d !\n", k, nred);
  }

  *types_ptr = types;
  *atoms_ptr = atoms;
}

// internal coordinates q and Wilson B-matrix B_ij = d(q_i)/d(x_j)
void collect_internal_coordinates(FORCE_FIELD *ff,
				  double **internal_ptr, double **bmatrix_ptr, int *dimensions)
{
  int i;
  int I,J,K,L;
  NODE *node;
  BOND *bond;
  ANGLE *angle;
  TORSION *torsion;
  INVERSION *inversion;
  
  double *internal;
  double *bmatrix;
  int nred, ncar;
  int k;
  double fchain;
  double scale;
  
  // count number of redundant internal coordinates
  nred = count_internal_coordinates(ff);
  // number of cartesian coordinates
  ncar = 3*ff->nat;
  //
  dimensions[0] = nred;
  dimensions[1] = ncar;

  // allocate memory for redundant internal coordinates and
  // gradients of internal coordinates w/r/t cartesian coordinates
  internal = (double *) malloc(sizeof(double) * nred);
  bmatrix = (double *) malloc(sizeof(double) * ncar * nred);

  // initialize memory to zero
  memset(internal, 0, sizeof(double) * nred);
  memset(bmatrix, 0, sizeof(double) * ncar * nred);
  
  k=0;  // enumerates redundant internal coordinates
  // bonds
  for(node=ff->bond_list->head; node!=NULL; node=node->next) {
    bond = (BOND *) node->element;
    I = bond->I; J = bond->J;
    // k-th internal coordinate is bond length
    internal[k] = bond->r;
    // derivatives of bond length w/r/t cartesian coordinates
    // of atoms I and J
    for(i=0;i<3;i++) {
      bmatrix[k*ncar+3*I+i]  = bond->derI[i];
      bmatrix[k*ncar+3*J+i] -= bond->derI[i];
    }
    k++;
  }

  // To avoid dividing by 0 in the chain rule for angles of 0 and 180 degrees
  // the definition of the valence and torsion angles is modified as
  //     theta = arccos(scale*cos(theta))
  scale = 1.0 - 1.0e-12; //1.0 - 1.0e-10;

  //  angles
  for(node=ff->angle_list->head; node!=NULL; node=node->next) {
    angle = (ANGLE *) node->element;
    I = angle->I; J = angle->J; K = angle->K;

    // use angle in radians as internal coordinate
    internal[k] = acos(scale*angle->cosTh);
    // factor from chain rule
    //   d(theta)/dx = d(arccos(cos(theta)))/dx = arccos'(cos(theta)) d(cos(theta))/dx
    //               = fchain * d(cos(theta))/dx
    fchain = -1.0/sqrt(1.0-pow(scale*angle->cosTh,2));

    /*
    // use cosine of angle as internal coordinate
    internal[k] = angle->cosTh;
    fchain = 1.0;
    */
    for(i=0;i<3;i++) {
      bmatrix[k*ncar+3*I+i] = fchain * angle->derI[i];
      bmatrix[k*ncar+3*J+i] = fchain * angle->derJ[i];
      bmatrix[k*ncar+3*K+i] = fchain * angle->derK[i];
    }
    k++;
  }
  // torsions
  for(node=ff->torsion_list->head; node!=NULL; node=node->next) {
    torsion = (TORSION *) node->element;
    I = torsion->I; J = torsion->J; K = torsion->K; L = torsion->L;

    /*
       Using the torsion angle as internal coordinate is more
       intuitive but we get coordinate singularities at 
       phi = 0 or 180 deg. Maybe one can also use
       cos(phi) directly.
    */

    // use angle 
    internal[k] = acos(scale*torsion->cosPh);
    // factor from chain rule
    //   d(theta)/dx = d(arccos(cos(theta)))/dx = arccos'(cos(theta)) d(cos(theta))/dx
    //               = fchain * d(cos(theta))/dx
    fchain = -1.0/sqrt(1.0-pow(scale*torsion->cosPh,2));

    /*
    // use cosine 
    internal[k] = torsion->cosPh;
    fchain = 1.0;
    */
    for(i=0;i<3;i++) {
      bmatrix[k*ncar+3*I+i] = fchain * torsion->derI[i];
      bmatrix[k*ncar+3*J+i] = fchain * torsion->derJ[i];
      bmatrix[k*ncar+3*K+i] = fchain * torsion->derK[i];
      bmatrix[k*ncar+3*L+i] = fchain * torsion->derL[i];
    }
    k++;
  }

  // inversions
  for(node=ff->inversion_list->head; node!=NULL; node=node->next) {
    inversion = (INVERSION *) node->element;
    I = inversion->I; J = inversion->J; K = inversion->K; L = inversion->L;

    // use angle 
    internal[k] = asin(scale*inversion->sinPs);
    // factor from chain rule
    //   d(theta)/dx = d(arcsin(sin(theta)))/dx = arcsin'(sin(theta)) d(sin(theta))/dx
    //               = fchain * d(sin(theta))/dx
    fchain = 1.0/sqrt(1.0-pow(scale*inversion->sinPs,2));

    /*
    // use sine
    internal[k] = inversion->sinPs;
    fchain = 1.0;
    */
    for(i=0;i<3;i++) {
      bmatrix[k*ncar+3*I+i] = fchain * inversion->derI[i];
      bmatrix[k*ncar+3*J+i] = fchain * inversion->derJ[i];
      bmatrix[k*ncar+3*K+i] = fchain * inversion->derK[i];
      bmatrix[k*ncar+3*L+i] = fchain * inversion->derL[i];
    }
    k++;
  }

  if (k != nred) {
    printf("BUG: k= %d  !=  nred= %d !\n", k, nred);
  }
  
  *internal_ptr = internal;
  *bmatrix_ptr = bmatrix;
}
