#include "Python.h"
#include "arrayobject.h"
#include <cstdio>
#include <math.h>
#include <xc.h>

typedef struct {
    PyObject_HEAD
    xc_func_type func;
} libXCFunctional;

static int
libXCFunctional_init(libXCFunctional *self, PyObject *args, PyObject *kwds)
{
  char *xc_name;
  int func_id;
  if (!PyArg_ParseTuple(args, "i", &func_id)) {
    PyErr_Clear();
    if (!PyArg_ParseTuple(args, "s", &xc_name)) {
      PyErr_SetString(PyExc_ValueError, "Expected the name or ID of an exchange correlation functional.");
      return -1;
    } else {
      func_id = XC(functional_get_number)(xc_name);
    }
  }

  if (xc_func_init(&(self->func), func_id, XC_UNPOLARIZED) != 0) {
    xc_func_init(&(self->func), 1, XC_UNPOLARIZED);
    PyErr_SetString(PyExc_ValueError, "Unknown exchange correlation potential. To get a list of all possibilities run pylibxc/pylibxc.py");
    return -1;
  }
  return 0;
}

static void
libXCFunctional_destroy(libXCFunctional* self)
{
  xc_func_end(&(self->func));
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *libXCFunctional_exc(PyObject *self, PyObject *args)
{
  PyArrayObject *PyRho, *PySigma;
  PyObject *PyExc;
  double *rho, *sigma, *exc;
  if (!PyArg_ParseTuple(args, "O!O!",
		       &PyArray_Type, &PyRho,
		       &PyArray_Type, &PySigma)) return NULL;
  if (PyRho == NULL) return NULL;
  if (PySigma == NULL) return NULL;
  // need contiguous memory
  PyRho = (PyArrayObject *) PyArray_GETCONTIGUOUS(PyRho);
  PySigma = (PyArrayObject *) PyArray_GETCONTIGUOUS(PySigma);
  int npoints = PyRho->dimensions[0];

  rho = (double *) PyRho->data;
  sigma = (double *) PySigma->data;
  exc = (double *) malloc(sizeof(double)*npoints);

  //
  xc_func_type func = ((libXCFunctional *) self)->func;
  switch(func.info->family) {
  case XC_FAMILY_LDA:
    xc_lda_exc(&func, npoints, rho, exc);
    break;
  case XC_FAMILY_GGA:
  case XC_FAMILY_HYB_GGA:
    xc_gga_exc(&func, npoints, rho, sigma, exc);
    break;
  }
  //
  
  npy_intp pydims[1];
  pydims[0] = (npy_intp) npoints;

  PyExc = PyArray_SimpleNewFromData(1, pydims, NPY_DOUBLE, exc);
  ((PyArrayObject *) PyExc)->flags |= NPY_OWNDATA;

  return PyExc;
}

static PyObject *libXCFunctional_vxc(PyObject *self, PyObject *args)
{
  PyArrayObject *PyRho, *PySigma;
  PyObject *PyVxc;
  double *rho, *sigma, *vxc, *vsigma;
  if (!PyArg_ParseTuple(args, "O!O!",
		       &PyArray_Type, &PyRho,
		       &PyArray_Type, &PySigma)) return NULL;
  if (PyRho == NULL) return NULL;
  if (PySigma == NULL) return NULL;
  // need contiguous memory
  PyRho = (PyArrayObject *) PyArray_GETCONTIGUOUS(PyRho);
  PySigma = (PyArrayObject *) PyArray_GETCONTIGUOUS(PySigma);
  int npoints = PyRho->dimensions[0];

  rho = (double *) PyRho->data;
  sigma = (double *) PySigma->data;
  vxc = (double *) malloc(sizeof(double)*npoints);
  vsigma = (double *) malloc(sizeof(double)*npoints);

  //
  xc_func_type func = ((libXCFunctional *) self)->func;
  switch(func.info->family) {
  case XC_FAMILY_LDA:
    xc_lda_vxc(&func, npoints, rho, vxc);
    break;
  case XC_FAMILY_GGA:
  case XC_FAMILY_HYB_GGA:
    xc_gga_vxc(&func, npoints, rho, sigma, vxc, vsigma);
    break;
  }
  //
  
  npy_intp pydims[1];
  pydims[0] = (npy_intp) npoints;

  PyVxc = PyArray_SimpleNewFromData(1, pydims, NPY_DOUBLE, vxc);
  ((PyArrayObject *) PyVxc)->flags |= NPY_OWNDATA;

  return PyVxc;
}

static PyObject *libXCFunctional_fxc(PyObject *self, PyObject *args)
{
  PyArrayObject *PyRho, *PySigma;
  PyObject *PyFxc;
  double *rho, *sigma, *v2rho2, *v2rhosigma, *v2sigma2;
  if (!PyArg_ParseTuple(args, "O!O!",
		       &PyArray_Type, &PyRho,
		       &PyArray_Type, &PySigma)) return NULL;
  if (PyRho == NULL) return NULL;
  if (PySigma == NULL) return NULL;
  // need contiguous memory
  PyRho = (PyArrayObject *) PyArray_GETCONTIGUOUS(PyRho);
  PySigma = (PyArrayObject *) PyArray_GETCONTIGUOUS(PySigma);
  int npoints = PyRho->dimensions[0];

  rho = (double *) PyRho->data;
  sigma = (double *) PySigma->data;
  v2rho2 = (double *) malloc(sizeof(double)*npoints);
  v2rhosigma = (double *) malloc(sizeof(double)*npoints);
  v2sigma2 = (double *) malloc(sizeof(double)*npoints);

  //
  xc_func_type func = ((libXCFunctional *) self)->func;
  switch(func.info->family) {
  case XC_FAMILY_LDA:
    xc_lda_fxc(&func, npoints, rho, v2rho2);
    break;
  case XC_FAMILY_GGA:
  case XC_FAMILY_HYB_GGA:
    xc_gga_fxc(&func, npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);
    break;
  }
  //
  
  npy_intp pydims[1];
  pydims[0] = (npy_intp) npoints;

  PyFxc = PyArray_SimpleNewFromData(1, pydims, NPY_DOUBLE, v2rho2);
  ((PyArrayObject *) PyFxc)->flags |= NPY_OWNDATA;

  return PyFxc;
}

/* copied from xc-info.c */
const char *get_kind(const xc_func_type *func) {
  switch(func->info->kind) {
  case(XC_EXCHANGE):
    return "XC_EXCHANGE";

  case(XC_CORRELATION):
    return "XC_CORRELATION";

  case(XC_EXCHANGE_CORRELATION):
   return "XC_EXCHANGE_CORRELATION";

  case(XC_KINETIC):
   return "XC_KINETIC";

  default:
    printf("Internal error in get_kind.\n");
    return "";
  }
}

const char *get_family(const xc_func_type *func) {
  switch(func->info->family) {
  case(XC_FAMILY_UNKNOWN):
    return "XC_FAMILY_UNKNOWN";

  case(XC_FAMILY_LDA):
    return "XC_FAMILY_LDA";

  case(XC_FAMILY_GGA):
    return "XC_FAMILY_GGA";

  case(XC_FAMILY_MGGA):
    return "XC_FAMILY_MGGA";

  case(XC_FAMILY_LCA):
    return "XC_FAMILY_LCA";

  case(XC_FAMILY_OEP):
    return "XC_FAMILY_OEP";

  case(XC_FAMILY_HYB_GGA):
    return "XC_FAMILY_HYB_GGA";

  case(XC_FAMILY_HYB_MGGA):
    return "XC_FAMILY_HYB_MGGA";

   default:
    printf("Internal error in get_family.\n");
    return "";
  }
}
/**/

static PyObject *libXCFunctional_description(PyObject *self)
{
  PyObject *PyS;
  char *fname;
  char s[5000];
  int n;

  xc_func_type func = ((libXCFunctional *) self)->func;
  int func_id = func.info->number;
  /* Get functional name */
  fname=XC(functional_get_name)(func_id);

  /* Print out info */
  n=  sprintf(s, "%10s: %-20i\t%10s: %-25s\n","func_id",func_id,"name",fname);
  n+= sprintf(s+n,"%10s: %-20s\t%10s: %-25s\n","family",get_family(&func),"kind",get_kind(&func));
  n+= sprintf(s+n,"%10s: %s\n","comment",func.info->name);

  n+= sprintf(s+n,"\nReference(s):\n");
  n+= sprintf(s+n,"%s\n",func.info->refs);

  free(fname);

  PyS = PyString_FromString(s);
  return PyS;
}

static PyMethodDef libXCFunctional_methods[] = {
  {"exc", (PyCFunction) libXCFunctional_exc, METH_VARARGS,
   "evaluate exchange correlation energy at given density"
  },
  {"vxc", (PyCFunction) libXCFunctional_vxc, METH_VARARGS,
   "evaluate exchange correlation potential dexc/dn at given density"
  },
  {"fxc", (PyCFunction) libXCFunctional_fxc, METH_VARARGS,
   "evaluate exchange correlation kernel d^2(exc)/d(n)^2 at given density"
  },
  {"description", (PyCFunction) libXCFunctional_description, METH_NOARGS,
   "Short description of functional with references",
  },
  {NULL}  /* Sentinel */
};

static PyTypeObject libXCFunctional_Type = {
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "_pylibxc.libXCFunctional",             /*tp_name*/
  sizeof(libXCFunctional), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor) libXCFunctional_destroy,  /*tp_dealloc*/
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
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
  "exchange correlation functional",           /* tp_doc */
  0,               /* tp_traverse */
  0,               /* tp_clear */
  0,               /* tp_richcompare */
  0,               /* tp_weaklistoffset */
  0,               /* tp_iter */
  0,               /* tp_iternext */
  libXCFunctional_methods,             /* tp_methods */
  0,               /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc) libXCFunctional_init,      /* tp_init */
  0,                         /* tp_alloc */
  0,          /* tp_new */
};

//static PyObject *list_functionals(PyObject *self, PyObject *args)
//{
//  
//}


/* #### Globals #################################### */

/* ==== Set up the methods table ====================== */
static PyMethodDef pylibxc_methods[] = {
  //	{"list_functionals", list_functionals, METH_VARARGS},
	{NULL, NULL}     /* Sentinel - marks the end of this structure */
};

/* ==== Initialize the C_test functions ====================== */
// Module name must be _C_arraytest in compile and linked 
PyMODINIT_FUNC init_pylibxc()  {

  PyObject *m;

  libXCFunctional_Type.tp_new = PyType_GenericNew;
  if (PyType_Ready(&libXCFunctional_Type) < 0)
    return;
  
  m = Py_InitModule3("_pylibxc", pylibxc_methods, "python wrapper around libxc");
  if (m == NULL)
    return;

  Py_INCREF(&libXCFunctional_Type);
  PyModule_AddObject(m, "libXCFunctional", (PyObject *) &libXCFunctional_Type);

  import_array();  // Must be present for NumPy.  Called first after above line.
}

