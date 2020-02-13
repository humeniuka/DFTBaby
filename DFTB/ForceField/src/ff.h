
// forward declarations
typedef struct FORCE_FIELD FORCE_FIELD;
typedef struct ATOM ATOM;
typedef struct CRYSTAL CRYSTAL;
typedef struct CHROMOPHORE CHROMOPHORE;

FORCE_FIELD *build_force_field(int nat, ATOM *atoms, CRYSTAL *crystal,
			       int nchromo, CHROMOPHORE *chromophores,
			       int *conmat, 
			       int verbose);

void set_current_state(FORCE_FIELD *ff, int state);
void evaluate_force_field(FORCE_FIELD *ff);
double collect_force_field(FORCE_FIELD *ff);
void evaluate_transition_dipoles(FORCE_FIELD *ff);
void delete_force_field(FORCE_FIELD *ff);

void local_chromophore_axes(FORCE_FIELD *ff);
  
// access functions
void set_force_field_coordinates(FORCE_FIELD *ff, int n, double x, double y, double z);
void get_force_field_gradient(FORCE_FIELD *ff, int, double *gradx, double *grady, double *gradz);
int get_number_of_atoms(FORCE_FIELD *ff);
int count_internal_coordinates(FORCE_FIELD *ff);
int get_number_of_excitons(FORCE_FIELD *ff);
// eigenenergies and coefficients of exciton wavefunctions
void get_exciton_energy(FORCE_FIELD *ff, double *energy, int state);
void get_exciton_coefficient(FORCE_FIELD *ff, double *coef, int i, int state);
// electric and magnetic transition dipoles
void get_transition_dipoles(FORCE_FIELD *ff,
			    double *tx, double *ty, double *tz,
			    double *mx, double *my, double *mz, int b);

// 
void get_coordinate_definitions(FORCE_FIELD *ff,
				char **types_ptr,
				int **atoms_ptr,
				int *length);
void collect_internal_coordinates(FORCE_FIELD *ff,
				  double **internal_ptr,
				  double **bmatrix_ptr,
				  int *dimensions);

