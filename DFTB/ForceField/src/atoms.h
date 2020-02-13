
typedef struct ATOM {
  int atom_type;
  double x,y,z;
  // partial charge
  double charge;
  // covalent radius
  double covR0;
  // hybridization
  int hyb;
  // Is this atom capable of donating or accepting hydrogen bonds?
  int hDA;
  // periodic group
  int group;
  // Dreiding force field parameters
  // bond
  double valR0;
  // angle
  double Te;    // equilibrium angle
  double cosTe; // cosine of equilibrium angle
  double Ce;    // force constant for angle bend
  int linear;  // 1 if angle = 180 deg
  // torsion
  double torsion_barrier;
  double torsion_periodicity;
  double torsion_phi0;
  // van der Waals
  double vdWR0;
  double vdWD0;
  // gradient on atom
  double grad[3];
} ATOM;

typedef struct
REPLICA_ATOM {
  int index; // index into the atoms inside central unit cell
  int central; // 1: atom lies in central unit cell, 0: atom lies in one of the neighbouring cells
  double tx,ty,tz; // translation vector
} REPLICA_ATOM;

typedef struct CRYSTAL {
  // components of lattice vectors
  double ax,ay,az; 
  double bx,by,bz;
  double cx,cy,cz;
} CRYSTAL;

// A chromophore is a collection of atoms which have transition charges associated
// with them for different excited states. Only the indeces to the atoms are stored.
typedef struct CHROMOPHORE {
  int nat;  // number of atoms
  int nst; // number of excited states
  int *indeces; // indeces of atoms, shape (nat)
  double *transition_charges; // transition charges (q) and magnetic dipoles (mx, my, mz)
                              // for each atom and state, shape (nat,4*nst)
                              //    q(1) mx(1) my(1) mz(1)
                              //    q(2) mx(2) my(2) mz(2)
                              //    ...
                              //    q(nat) mx(nat) my(nat) mz(nat)
  
  double *excitation_energies; // excitation energies for each state in Hartree, shape (nst)
  // axis for local molecular coordinate system
  double axes[3][3];   // axes[i][:] is the i-th axis
} CHROMOPHORE;

