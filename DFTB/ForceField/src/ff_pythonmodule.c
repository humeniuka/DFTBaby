#include <Python.h>
#include "numpy/arrayobject.h"
#include "structmember.h"
#include <stdio.h>

#include "input.h"
#include "ff.h"

typedef struct {
  PyObject_HEAD
  // object data
  FORCE_FIELD *ff;
} ForceField;

static void
ForceField_dealloc(ForceField *self)
{
  self->ob_type->tp_free((PyObject*)self);
  if (self->ff != NULL) {
    delete_force_field(self->ff);
  }
}

static PyObject *
ForceField_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  ForceField *self;

  self = (ForceField *)type->tp_alloc(type, 0);
  
  return (PyObject *)self;
}

static int
ForceField_init_from_file(ForceField *self, PyObject *args, PyObject *kwds)
{
  char *ff_file_name;
  ATOM *atoms;
  int n, nat;
  CRYSTAL crystal;

  self->ff = NULL;
  
  if (!PyArg_ParseTuple(args, "s", &ff_file_name)) {
    return -1;
  }
  // read force field from file
  printf("reading force field definition from file %s\n", ff_file_name);
  FILE *fh = fopen(ff_file_name, "r");
  if (fh == NULL) {
    PyErr_SetString(PyExc_IOError, "File not found!");
    return -1;
  }
  atoms = read_force_field(fh, &nat, &crystal);
  fclose(fh);
  if (atoms == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "Syntax error in force field file!");
    return -1;
  }
  assign_parameters(atoms, nat, 1);
  
  //
  printf("Atom in Unit Cell:\n");
  for (n=0; n<nat; n++) {
    printf("%d  %+f %+f %+f\n", atoms[n].atom_type, atoms[n].x, atoms[n].y, atoms[n].z);
  }
  printf("Lattice vectors:\n");
  printf("a  %f %f %f\n", crystal.ax, crystal.ay, crystal.az);
  printf("b  %f %f %f\n", crystal.bx, crystal.by, crystal.bz);
  printf("c  %f %f %f\n", crystal.cx, crystal.cy, crystal.cz);
  //
  
  self->ff = build_force_field(nat, atoms, &crystal, 0, NULL, NULL, 1); 
  
  return 0;
}

static int
ForceField_init_from_lists(ForceField *self, PyObject *args, PyObject *kwds)
{
  PyObject *py_atomlist;
  PyObject *py_typelist;
  PyObject *py_chargelist;
  PyObject *py_lattice_vectors;
  PyObject *py_chromophorelist;
  PyObject *py_atomtuple, *py_poslist, *py_posseq;
  PyObject *py_chromotuple, *py_indexlist, *py_energylist, *py_statelist;
  PyArrayObject *py_conmatarray = NULL;
  int verbose = 0;
  
  ATOM *atoms;
  ATOM *atom;
  int n, nat;
  CRYSTAL crystal;
  int nchromo;
  CHROMOPHORE *chromophores;
  CHROMOPHORE *chromo;
  int i,j;
  int *conmat;
  
  static char *kwlist[] = {"atomlist", "atomtypes", "partial_charges",
			   "lattice_vectors", "chromophores", "connectivity_matrix", "verbose", NULL};
  
  self->ff = NULL;
  
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!O!O!O!|O!i", kwlist,
			&PyList_Type, &py_atomlist, 
			&PyList_Type, &py_typelist,
			&PyList_Type, &py_chargelist,
			&PyList_Type, &py_lattice_vectors,
		        &PyList_Type, &py_chromophorelist,
			&PyArray_Type, &py_conmatarray, 
			&verbose	   )) {
    return -1;
  }
  nat = (int) PyList_Size(py_atomlist);
  if (nat != PyList_Size(py_typelist)) {
    PyErr_SetString(PyExc_ValueError, "There should be one atomtype for each atom!");
    return -1;
  }
  if (nat != PyList_Size(py_chargelist)) {
    PyErr_SetString(PyExc_ValueError, "There should be one partial charge for each atom!");
    return -1;
  }
  if (3 != PyList_Size(py_lattice_vectors)) {
    PyErr_SetString(PyExc_ValueError, "3 lattice vectors are required!");
    return -1;
  }
  // extract all atoms
  atoms = (ATOM *) malloc(sizeof(ATOM)*nat);
  for(n=0;n<nat;n++) {
    py_atomtuple = PyList_GetItem(py_atomlist, (Py_ssize_t) n);
    py_poslist = PyTuple_GetItem(py_atomtuple, (Py_ssize_t) 1);
    py_posseq = PySequence_Fast(py_poslist, "Atom positions should be a list [x,y,z] or tuple (x,y,z)!");

    atom = &(atoms[n]);

    atom->x = PyFloat_AsDouble(PySequence_GetItem(py_posseq, (Py_ssize_t) 0));
    atom->y = PyFloat_AsDouble(PySequence_GetItem(py_posseq, (Py_ssize_t) 1));
    atom->z = PyFloat_AsDouble(PySequence_GetItem(py_posseq, (Py_ssize_t) 2));
    
    Py_DECREF(py_posseq);
    // atom type
    atom->atom_type = (int) PyInt_AsLong(PyList_GetItem(py_typelist, (Py_ssize_t) n));
    if ((atom->atom_type < 0) || (atom->atom_type >= dreiding_nr_atom_types)) {
      PyErr_SetString(PyExc_ValueError, "Invalid atom type!");
      free(atoms);
      return -1;
    }
    // partial charge
    atom->charge = (float) PyFloat_AsDouble(PyList_GetItem(py_chargelist, (Py_ssize_t) n));
  }
  // read lattice vectors
  py_poslist = PyList_GetItem(py_lattice_vectors,(Py_ssize_t) 0);
  crystal.ax = PyFloat_AsDouble(PyList_GetItem(py_poslist, (Py_ssize_t) 0));
  crystal.ay = PyFloat_AsDouble(PyList_GetItem(py_poslist, (Py_ssize_t) 1));
  crystal.az = PyFloat_AsDouble(PyList_GetItem(py_poslist, (Py_ssize_t) 2));
  py_poslist = PyList_GetItem(py_lattice_vectors,(Py_ssize_t) 1);
  crystal.bx = PyFloat_AsDouble(PyList_GetItem(py_poslist, (Py_ssize_t) 0));
  crystal.by = PyFloat_AsDouble(PyList_GetItem(py_poslist, (Py_ssize_t) 1));
  crystal.bz = PyFloat_AsDouble(PyList_GetItem(py_poslist, (Py_ssize_t) 2));
  py_poslist = PyList_GetItem(py_lattice_vectors,(Py_ssize_t) 2);
  crystal.cx = PyFloat_AsDouble(PyList_GetItem(py_poslist, (Py_ssize_t) 0));
  crystal.cy = PyFloat_AsDouble(PyList_GetItem(py_poslist, (Py_ssize_t) 1));
  crystal.cz = PyFloat_AsDouble(PyList_GetItem(py_poslist, (Py_ssize_t) 2));
  // assign atom types
  assign_parameters(atoms, nat, verbose);

  // If the option 'connectivity_matrix' is provided, bonds
  // are created between atoms I and J if connectivity_matrix[I,J] = 1, 
  // otherwise the connectivity is determined from the covalent
  // atom radii inside `build_force_field()`.
  if (py_conmatarray != NULL) {
    printf("overriding default connectivity matrix\n");
    // check connectivity matrix has the right dimensions
    if (py_conmatarray->descr->type_num != NPY_INT64
	   || py_conmatarray->nd != 2
	   || py_conmatarray->dimensions[0] != nat
	|| py_conmatarray->dimensions[1] != nat) {
      PyErr_SetString(PyExc_ValueError, "Connectivity matrix should be an integer array of shape (nat,nat) !");
      free(atoms);
      return -1;
    }
    // copy data from numpy array to C array
    conmat = (int *) malloc(sizeof(int)*nat*nat);
    for(i=0; i<nat; i++) {
      for(j=0; j<nat; j++) {
	conmat[i*nat+j] = *(int *) PyArray_GETPTR2(py_conmatarray, i,j);
      }
    }
  } else {
    conmat = NULL;
  }
    
  //
  if (verbose > 0) {
    printf("Atom in Unit Cell:\n");
    for (n=0; n<nat; n++) {
      printf("%d  %+f %+f %+f   Q=%+f\n", atoms[n].atom_type, atoms[n].x, atoms[n].y, atoms[n].z, atoms[n].charge);
    }
    printf("Lattice vectors:\n");
    printf("a  %f %f %f\n", crystal.ax, crystal.ay, crystal.az);
    printf("b  %f %f %f\n", crystal.bx, crystal.by, crystal.bz);
    printf("c  %f %f %f\n", crystal.cx, crystal.cy, crystal.cz);
  }

  // extract chromophore data (if available)
  nchromo = (int) PyList_Size(py_chromophorelist);
  if (nchromo > 0) {
    chromophores = (CHROMOPHORE *) malloc(sizeof(CHROMOPHORE)*nchromo);
    for(n=0; n<nchromo; n++) {
      py_chromotuple = PyList_GetItem(py_chromophorelist, (Py_ssize_t) n);
      py_indexlist  = PyTuple_GetItem(py_chromotuple, (Py_ssize_t) 0);
      py_energylist = PyTuple_GetItem(py_chromotuple, (Py_ssize_t) 1);
      py_chargelist = PyTuple_GetItem(py_chromotuple, (Py_ssize_t) 2);
      
      chromo = &(chromophores[n]);
      chromo->nat = PyList_Size(py_indexlist);  // number of atoms belonging to this chromophore
      chromo->nst = PyList_Size(py_energylist); // number of excited states this chromophore can be in
      // allocate memory
      chromo->indeces = (int *) malloc(sizeof(int)*chromo->nat);
      chromo->excitation_energies = (double *) malloc(sizeof(double)*chromo->nst);
      chromo->transition_charges = (double *) malloc(sizeof(double)*4*chromo->nat*chromo->nst);
      // copy chromophore data
      for(j=0; j<chromo->nst; j++) {
	chromo->excitation_energies[j] = PyFloat_AsDouble(PyList_GetItem(py_energylist, (Py_ssize_t) j));
      }
      for(i=0; i<chromo->nat; i++) { // iterate over atoms
	chromo->indeces[i] = (int) PyInt_AsLong(PyList_GetItem(py_indexlist, (Py_ssize_t) i));
	py_statelist = PyList_GetItem(py_chargelist, (Py_ssize_t) i);
	for(j=0; j<4*chromo->nst; j++) { // iterate over excited states,
	                                 // with four numbers per state (q mx my mz)
	  chromo->transition_charges[i*4*chromo->nst+j] = PyFloat_AsDouble(PyList_GetItem(py_statelist, (Py_ssize_t) j));
	}
      }
    }
  // 
  } else {
    chromophores = NULL;
  }


  // create lists of bonds, angles, dihedrals, improper torsions and vdW interactions
  self->ff = build_force_field(nat, atoms, &crystal,
			       nchromo, chromophores,
			       conmat, 
			       verbose);
  // determine local coordinate system around each chromophore
  local_chromophore_axes(self->ff);

  return 0;
}

static PyObject *
ForceField_getEnergyAndGradient(ForceField *self, PyObject *args)
{
  PyArrayObject *coords_arr, *grads_arr;
  int state;
  int dimensions[1];
  int ndim, n, nat;
  double *x_ptr, *y_ptr, *z_ptr;
  double x,y,z;
  double energy;
  
  if (!PyArg_ParseTuple(args, "O!i", &PyArray_Type, &coords_arr, &state)) {
    return NULL;
  }

  if (coords_arr->nd != 1 || coords_arr->descr->type_num != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "Coordinates vector must be one-dimensional and of type double!");
    return NULL;
  }
  ndim = coords_arr->dimensions[0];
  nat = get_number_of_atoms(self->ff);
  if (ndim != 3*nat) {
    PyErr_SetString(PyExc_ValueError, "Wrong number of coordinates!");
    return NULL;
  }
  // update coordinates of atoms in the force field
  for(n=0; n<nat; n++) {
    x = *(double *) (coords_arr->data + (3*n+0)*coords_arr->strides[0]);
    y = *(double *) (coords_arr->data + (3*n+1)*coords_arr->strides[0]);
    z = *(double *) (coords_arr->data + (3*n+2)*coords_arr->strides[0]);
    // set coordinates of n-th atom
    set_force_field_coordinates(self->ff, n, x,y,z);
  }

  // set current electronic state: 0 - ground state
  //                               1 - first exciton state
  set_current_state(self->ff, state);
  evaluate_force_field(self->ff);
  energy = collect_force_field(self->ff);

  // copy gradient back
  dimensions[0] = ndim;
  grads_arr = (PyArrayObject *) PyArray_FromDims(1, dimensions, PyArray_DOUBLE);
  if (grads_arr == NULL) {
    return NULL;
  }
  
  for(n=0; n<nat; n++) {
    x_ptr = (double *) (grads_arr->data + (3*n+0)*grads_arr->strides[0]);
    y_ptr = (double *) (grads_arr->data + (3*n+1)*grads_arr->strides[0]);
    z_ptr = (double *) (grads_arr->data + (3*n+2)*grads_arr->strides[0]);

    // get gradient on n-th atom
    get_force_field_gradient(self->ff, n, &x, &y, &z);
    *x_ptr = x;
    *y_ptr = y;
    *z_ptr = z;
  }

  return Py_BuildValue("(fN)", (double) energy, (PyObject *) grads_arr);
}

static PyObject *
ForceField_getRedundantInternalCoordinates(ForceField *self, PyObject *args)
{
  PyArrayObject *coords_arr, *internal_arr, *bmatrix_arr;
  int state;
  int dimensions[2];
  npy_intp npy_dim[2];
  int ndim, n, nat;
  double x,y,z;
  double *internal, *bmatrix;
  
  if (!PyArg_ParseTuple(args, "O!i", &PyArray_Type, &coords_arr, &state)) {
    return NULL;
  }

  if (coords_arr->nd != 1 || coords_arr->descr->type_num != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "Coordinates vector must be one-dimensional and of type double!");
    return NULL;
  }
  ndim = coords_arr->dimensions[0];
  nat = get_number_of_atoms(self->ff);
  if (ndim != 3*nat) {
    PyErr_SetString(PyExc_ValueError, "Wrong number of coordinates!");
    return NULL;
  }
  // update coordinates of atoms in the force field
  for(n=0; n<nat; n++) {
    x = *(double *) (coords_arr->data + (3*n+0)*coords_arr->strides[0]);
    y = *(double *) (coords_arr->data + (3*n+1)*coords_arr->strides[0]);
    z = *(double *) (coords_arr->data + (3*n+2)*coords_arr->strides[0]);
    // set coordinates of n-th atom
    set_force_field_coordinates(self->ff, n, x,y,z);
  }

  // set current electronic state: 0 - ground state
  //                               1 - first exciton state
  set_current_state(self->ff, state);
  evaluate_force_field(self->ff);
  // compute redundant internal coordinates and Wilson B-matrix
  collect_internal_coordinates(self->ff, &internal, &bmatrix, dimensions);
  npy_dim[0] = (npy_intp) dimensions[0];
  npy_dim[1] = (npy_intp) dimensions[1];
  
  // wrap memory allocated inside this function into
  // numpy array
  internal_arr = (PyArrayObject *) PyArray_SimpleNewFromData(1, npy_dim, NPY_DOUBLE, (void *) internal);
  if (internal_arr == NULL) {
    return NULL;
  }
  internal_arr->flags |= NPY_ARRAY_OWNDATA | NPY_ARRAY_C_CONTIGUOUS;
  
  bmatrix_arr = (PyArrayObject *) PyArray_SimpleNewFromData(2, npy_dim, NPY_DOUBLE, (void *) bmatrix);
  if (bmatrix_arr == NULL) {
    return NULL;
  }
  bmatrix_arr->flags |= NPY_ARRAY_OWNDATA | NPY_ARRAY_C_CONTIGUOUS;

  return Py_BuildValue("(NN)", (PyObject *) internal_arr, (PyObject *) bmatrix_arr);
}

static PyObject *
ForceField_getInternalCoordinateDefinitions(ForceField *self, PyObject *args)
{
  PyArrayObject *types_arr, *atoms_arr;
  npy_intp npy_dim[1];
  int ndim;

  char *types;
  int *atoms;

  // types is a list with character labels for each internal coordinates ('B' for bond, 'A' for angle, 'D' for dihedral)
  // and atoms contains the atom indices defining each bond, angle, etc. in blocks of four
  get_coordinate_definitions(self->ff, &types, &atoms, &ndim);

  // wrap memory allocated inside this function into
  // numpy array
  npy_dim[0] = (npy_intp) ndim;
  types_arr = (PyArrayObject *) PyArray_SimpleNewFromData(1, npy_dim, NPY_CHAR, (void *) types);
  if (types_arr == NULL) {
    return NULL;
  }
  types_arr->flags |= NPY_ARRAY_OWNDATA | NPY_ARRAY_C_CONTIGUOUS;

  // atom indices, in blocks of four
  npy_dim[0] = (npy_intp) (4*ndim);
  atoms_arr = (PyArrayObject *) PyArray_SimpleNewFromData(1, npy_dim, NPY_INT, (void *) atoms);
  if (atoms_arr == NULL) {
    return NULL;
  }
  atoms_arr->flags |= NPY_ARRAY_OWNDATA | NPY_ARRAY_C_CONTIGUOUS;
  
  return Py_BuildValue("(NN)", (PyObject *) types_arr, (PyObject *) atoms_arr);
}

static PyObject *
ForceField_getExcitonStates(ForceField *self, PyObject *args)
{
  PyArrayObject *energies_arr, *coeffs_arr;
  int nh, dim1[1], dim2[2];
  int i, j;
  double *en_ptr, *coef_ptr;
  
  // copy excitation energies
  nh = get_number_of_excitons(self->ff);  // number of exciton states
  dim1[0] = nh;
  energies_arr = (PyArrayObject *) PyArray_FromDims(1, dim1, PyArray_DOUBLE);
  if (energies_arr == NULL) {
    return NULL;
  }
  for(j=0; j<nh; j++) {
    en_ptr = (double *) (energies_arr->data + j*energies_arr->strides[0]);
    get_exciton_energy(self->ff, en_ptr, j);
  }
  // copy coefficients of exciton wavefunctions
  dim2[0] = nh;
  dim2[1] = nh;
  coeffs_arr = (PyArrayObject *) PyArray_FromDims(2, dim2, PyArray_DOUBLE);
  for(i=0; i<nh; i++) {
    for(j=0; j<nh; j++) {
      coef_ptr = (double *) (coeffs_arr->data + i*coeffs_arr->strides[0] + j*coeffs_arr->strides[1]);
      get_exciton_coefficient(self->ff, coef_ptr, i, j);
    }
  }
  
  return Py_BuildValue("(NN)", (PyObject *) energies_arr, (PyObject *) coeffs_arr);
}

static PyObject *
ForceField_getBasicTransitionDipoles(ForceField *self, PyObject *args)
{
  PyArrayObject *elec_tdip_arr, *magn_tdip_arr;
  int nh, dim2[2];
  int b; // enumerat exciton basis states
  double *tx_ptr, *ty_ptr, *tz_ptr;
  double *mx_ptr, *my_ptr, *mz_ptr;

  evaluate_transition_dipoles(self->ff);
  
  nh = get_number_of_excitons(self->ff);  // number of exciton states
  // copy transition dipole moments for S0 -> exciton basis vector
  dim2[0] = 3;
  dim2[1] = nh;
  elec_tdip_arr = (PyArrayObject *) PyArray_FromDims(2, dim2, PyArray_DOUBLE);
  magn_tdip_arr = (PyArrayObject *) PyArray_FromDims(2, dim2, PyArray_DOUBLE);
  for(b=0; b<nh; b++) {
    tx_ptr = (double *) (elec_tdip_arr->data + 0*elec_tdip_arr->strides[0] + b*elec_tdip_arr->strides[1]);
    ty_ptr = (double *) (elec_tdip_arr->data + 1*elec_tdip_arr->strides[0] + b*elec_tdip_arr->strides[1]);
    tz_ptr = (double *) (elec_tdip_arr->data + 2*elec_tdip_arr->strides[0] + b*elec_tdip_arr->strides[1]);

    mx_ptr = (double *) (magn_tdip_arr->data + 0*magn_tdip_arr->strides[0] + b*magn_tdip_arr->strides[1]);
    my_ptr = (double *) (magn_tdip_arr->data + 1*magn_tdip_arr->strides[0] + b*magn_tdip_arr->strides[1]);
    mz_ptr = (double *) (magn_tdip_arr->data + 2*magn_tdip_arr->strides[0] + b*magn_tdip_arr->strides[1]);

    get_transition_dipoles(self->ff,
			   tx_ptr, ty_ptr, tz_ptr,
			   mx_ptr, my_ptr, mz_ptr, b);
  }
  
  return Py_BuildValue("(NN)", (PyObject *) elec_tdip_arr, (PyObject *) magn_tdip_arr);
}


static PyMethodDef ForceField_methods[] = {
  {"getInternalCoordinateDefinitions", (PyCFunction)ForceField_getInternalCoordinateDefinitions, METH_VARARGS,
"definitions of internal coordinates in terms of the atom indices\n\
\n\
Returns:\n\
========\n\
types   : numpy array with character labels designating the coordinate type\n\
atoms   : atom indices defining each internal coordinate in blocks of four\n\
\n\
types[k] contains the character label for the k-th internal coordinate\n\
     'B' - bond\n\
     'A' - valence angle\n\
     'D' - dihedral angle\n\
\n\
Each internal coordinate is defined by at most 4 atom indices (starting from 0)\n\
    I J      -  bond\n\
    I J K    -  valence angle\n\
    I J K L  -  dihedral angle\n\
\n\
I = atoms[k*4], J = atoms[k*4+1], K = atoms[k*4+2], L = atoms[k*4+3]\n\
are the atom indices defining the k-th internal coordinate. For bonds\n\
and angles the indices K and L should not be used.\n\
\n\
"},
  {"getEnergyAndGradient", (PyCFunction)ForceField_getEnergyAndGradient, METH_VARARGS,
"Computes the energy and gradient on the atoms in the central unit cell.\n\
The number and order of atoms should remain the same as used when defining the force field.\n\n\
Parameters:\n\
===========\n\
coords: numpy array with shape (3*Nat), coords[3*i:3*(i+1)] are the x,y,z coordinates of the i-th atom in bohr\n\
\n\
Returns:\n\
========\n\
energy: total force field energy in Hartree\n\
grads: numpy array with shape (3*Nat), grads[3*i:3*(i+1)] is the gradient on the i-th atom in Hartree/bohr\n\
    "},
  {"getRedundantInternalCoordinates", (PyCFunction)ForceField_getRedundantInternalCoordinates, METH_VARARGS,
"Convert the cartesian coordinates into redundant internal coordinates and compute the\n\
Wilson B-matrix.\n\
\n\
Parameters:\n\
===========\n\
coords  : numpy array with shape (3*Nat), coords[3*i:3*(i+1)] are the x,y,z coordinates\n\
          of the i-th atom in bohr\n\
\n\
Returns:\n\
========\n\
internal  : numpy array with shape (M) with internal coordinates, first all bonds (in bohr), then\n\
            angles and torsions (in radians).\n\
bmatrix   : numpy array with shape (M,3*Nat) with Wilson B-matrix,\n\
            bmatrix[k,3*i:3*(i+1)] is the gradient of the k-th internal coordinate on\n\
            atom i.\n\
"},
  {"getExcitonStates", (PyCFunction)ForceField_getExcitonStates, METH_NOARGS,
"Fetches the eigenenergies and coefficients of the exciton wavefunctions.\n\
Calls to this functions must be preceeded by a call to .getEnergyAndGradient(...).\n\
\n\
Returns:\n\
========\n\
en     : numpy array with shape (Nst), eigen energies of exciton states (in Hartree)\n\
coeffs : numpy array with shape (Nst,Nst), coefficients of the exciton wavefunction\n\
              coeffs[:,i] are the eigenvectors belonging to the eigenenergy en[i]\n\
    "},
  {"getBasicTransitionDipoles", (PyCFunction)ForceField_getBasicTransitionDipoles, METH_NOARGS,
"Fetches the electric and magnetic transition dipole moments between the ground state (S0)\n\
and the excitonic basis states.\n\
\n\
Returns:\n\
========\n\
elec_tdip : numpy array with shape (3,Nst), electric transition dipoles\n\
magn_tdip : numpy array with shape (3,Nst), magnetic transition dipoles\n\
    "},
  {NULL, NULL, 0, NULL}  /* Sentinel */
};

static PyTypeObject ForceFieldType = {
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "ff.ForceField",             /*tp_name*/
  sizeof(ForceField),             /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)ForceField_dealloc, /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  "Builds a forces field with periodic boundary conditions.\n\
valence terms:\n\
  bond stretching, angle bending, torsion and spectroscopic inversion\n\
non-bonding terms:\n\
  van der Waals (Lennard Jones potential)\n\
The Coulomb interaction is not included, so all atoms should be approximately neutral.\n\
The parameters are taken from the DREIDING force field. \n\
\n\
So far only two atom-types exist:\n\
   0: hydrogen\n\
   1: aromatic carbon\n\
\n\
Parameters:\n\
===========\n\
atomlist: geometry which defined the bonding topology, \n\
  list of tuples (Zat,[x,y,z]) for each atom           \n\
atomtypes: list with indeces of atom types             \n\
lattice_vectors: 3 lists with translation vectors of lattice\n\
  [[ax,ay,az],[bx,by,bz],[cx,cy,cz]]\n\
\n\
Optional:\n\
=========\n\
connectivity_matrix: integer numpy array of shape (Nat,Nat) with connectivity \n\
  matrix that is used to construct bonds, angles and dihedrals.               \n\
  If atoms I and J are connected, connectivity_matrix[I,J] = 1, otherwise 0.  \n\
  If absent, atoms are assumed to be connected if their distance is smaller than \n\
  1.3 times the sum of the covalent radii.\n\
verbose: print information about the constructed force field\n\
",      /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  0,               /* tp_richcompare */
  0,               /* tp_weaklistoffset */
  0,               /* tp_iter */
  0,               /* tp_iternext */
  ForceField_methods,             /* tp_methods */
  0,             /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)ForceField_init_from_lists,      /* tp_init */
  0,                         /* tp_alloc */
  ForceField_new,                 /* tp_new */
};

static PyMethodDef module_methods[] = {
  {NULL, NULL, 0, NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initff(void)
{
  PyObject* m;

  if (PyType_Ready(&ForceFieldType) < 0)
    return;

  m = Py_InitModule3("ff", module_methods,
		     "Force Field for periodic QM/MM calculations of molecular crystals.");
  import_array();
  if (m == NULL)
    return;

  Py_INCREF(&ForceFieldType);
  PyModule_AddObject(m, "ForceField", (PyObject *)&ForceFieldType);
}

