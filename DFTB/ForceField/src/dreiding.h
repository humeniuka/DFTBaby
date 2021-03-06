// Parameters of Dreiding force field for H_ and C_R
/*
    Atom Types:
    ===========

*/

enum AtomType {
  H_, H__HB, H_b,
  B_3, B_2,
  C_3, C_R, C_2, C_1,
  N_3, N_R, N_2, N_1,
  O_3, O_R, O_2, O_1, F_,
  Al3, Si3, P_3, S_3, Cl,
  Ga3, Ge3, As3, Se3, Br,
  In3, Sn3, Sb3, Te3, I_,
  Na, Ca, Fe, Zn
};

enum Hybridization {
  // spR: resonance situation
  // monoval: monovalent atoms such as H, Cl, Fl
  // ionic: ionic bond for atoms such as Na^+, Ca^2+, Fe^2+, Zn^2+
  monoval, sp1, sp2, sp3, spR, ionic
};

// The parameters are listed in the same order as the atom types
const int dreiding_hyb[] = {monoval, monoval, monoval,
				     sp3, sp2,
				     sp3, spR, sp2, sp1,
				     sp3, spR, sp2, sp1,
				     sp3, spR, sp2, sp1, monoval,
				     sp3, sp3, sp3, sp3, monoval,
				     sp3, sp3, sp3, sp3, monoval,
				     sp3, sp3, sp3, sp3, monoval,
				     ionic, ionic, ionic, ionic};

const int periodic_group[] = { 1,  1, 1,
			       13, 13,
			       14, 14, 14, 14,
			       15, 15, 15, 15,
			       16, 16, 16, 16, 17,
			       13, 14, 15, 16, 17,
			       13, 14, 15, 16, 17,
			       13, 14, 15, 16, 17,
			       1, 2, 8, 12 };

// lengths in Angstrom, angles in degrees, energies in kcal/mol
const double covalent_radius[] = {0.31, 0.31, 0.31,
				  0.84, 0.84,
				  0.76, 0.73, 0.73, 0.69,
				  0.71, 0.71, 0.71, 0.71,
				  0.66, 0.66, 0.66, 0.66, 0.57,
				  1.21, 1.11, 1.07, 1.05, 1.02,
				  1.22, 1.20, 1.19, 1.20, 1.20,
				  1.42, 1.39, 1.39, 1.38, 1.39,
				  1.66, 1.76, 1.32, 1.22};

const double dreiding_valR0[]  = {0.330, 0.330, 0.510,
				  0.880, 0.790,
				  0.770, 0.700, 0.670, 0.602,
				  0.702, 0.650, 0.615, 0.556,
				  0.660, 0.660, 0.560, 0.528, 0.611,
				  1.047, 0.937, 0.890, 1.040, 0.997,
				  1.210, 1.210, 1.210, 1.210, 1.167,
				  1.390, 1.373, 1.432, 1.280, 1.360,
				  1.860, 1.940, 1.285, 1.330};

const double dreiding_angle[]  = {180.0, 180.0, 90.0,
				  109.471, 120.0,
				  109.471, 120.0, 120.0, 180.0,
				  106.7  , 120.0, 120.0, 180.0,
				  104.51 , 120.0, 120.0, 180.0, 180.0,
				  109.471, 109.471, 93.3, 92.1, 180.0,
				  109.471, 109.471, 92.1, 90.6, 180.0,
				  109.471, 109.471, 91.6, 90.3, 180.0,
				  90.0, 90.0, 90.0, 109.471};

const double dreiding_torsion_barrier[] = {0.0, 0.0, 0.0,
					   2.0, 0.0,
					   2.0, 25.0, 45.0, 0.0,
					   2.0, 25.0, 45.0, 0.0,
					   2.0, 25.0, 45.0, 0.0, 0.0,
					   2.0,  2.0,  2.0, 2.0, 0.0,
					   2.0,  2.0,  2.0, 2.0, 0.0,
					   2.0,  2.0,  2.0, 2.0, 0.0,
					   0.0,  0.0,  0.0, 0.0, 0.0};
const int dreiding_torsion_periodicity[] = {0, 0, 0,
					    3, 0,
					    3, 2, 2, 0,
					    3, 2, 2, 0,
					    2, 2, 2, 0, 0,
					    3, 3, 3, 2, 0,
					    3, 3, 3, 2, 0,
					    3, 3, 3, 2, 0,
					    0, 0, 0, 0};
const double dreiding_torsion_phi0[] = {0.0, 0.0, 0.0,
					180.0, 0.0,
					180.0, 180.0, 180.0, 0.0,
					180.0, 180.0, 180.0, 0.0,
					90.0,  180.0, 180.0, 0.0,  0.0,
					180.0, 180.0, 180.0, 90.0, 0.0,
					180.0, 180.0, 180.0, 90.0, 0.0,
					180.0, 180.0, 180.0, 90.0, 0.0,
					0.0,     0.0,   0.0,  0.0};

// non-bonding parameters
const double dreiding_vdWR0[]  = {3.195, 3.195, 3.195,
				  4.02, 4.02,
				  3.8983, 3.8983, 3.8983, 3.8983,
				  3.6621, 3.6621, 3.6621, 3.6621,
				  3.4046, 3.4046, 3.4046, 3.4046, 3.4720,
				  4.39, 4.27, 4.1500, 4.030, 3.9503,
				  4.39, 4.27, 4.15, 4.03, 3.95,
				  4.59, 4.47, 4.35, 4.23, 4.15,
				  3.144, 3.472, 4.54, 4.54};

const double dreiding_vdWD0[]  = {0.0152, 0.0152, 0.0001,
				  0.095, 0.095,
				  0.0951, 0.0951, 0.0951, 0.0951,
				  0.0774, 0.0774, 0.0774, 0.0774,
				  0.0957, 0.0957, 0.0957, 0.0957, 0.0725,
				  0.31, 0.31, 0.3200, 0.3440, 0.2833,
				  0.40, 0.40, 0.41, 0.43, 0.37,
				  0.55, 0.55, 0.55, 0.57, 0.51,
				  0.5, 0.05, 0.055, 0.055};

// atom types that can donate or accept hydrogen bonds are marked as 1 all others as 0
const int dreiding_hDA[] = {0, 0, 0,
			    0, 0,
			    0, 0, 0, 0,
			    1, 1, 1, 1,
			    1, 1, 1, 1, 1,
			    0, 0, 0, 0, 0,
			    0, 0, 0, 0, 0,
			    0, 0, 0, 0, 0,
			    0, 0, 0, 0};

