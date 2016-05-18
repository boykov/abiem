#include <Python.h>
#include <numpy/arrayobject.h>
#include <complex.h>

#include "ioAHMED.h"
#include "timeAHMED.h"
#include "vec3d.h"
#include "cmplx.h"

#include "bemcluster.h"
#include "bemblcluster.h"
#include "matgen_omp.h"


#include "H.h" // HCholesky
#include "solvers.h" // PCG
#include "bemmatrix.h"

using namespace AHMED;

class MATGENKERNEL{

public:
  //! Constructor for initialization
  MATGENKERNEL(const vec3d* vectors, const unsigned* op_perm, dcomp (*kernel)(unsigned , unsigned ))
    : vectors_(vectors), op_perm_(op_perm), kernel_(kernel) { }

  void cmpbl(const unsigned b1, const unsigned n1, const unsigned b2, const unsigned n2,
	     dcomp* pt) const {

#pragma omp critical (matgen_prog)
    {
      for(unsigned j=b2; j<b2+n2; ++j)
	for(unsigned i=b1; i<b1+n1; ++i)
	  *pt++ = (*kernel_)(op_perm_[i]+1, op_perm_[j]+1);
    }
  }

  //! returns expected size of the matrix entries
  double scale (const unsigned /*b1*/, const unsigned /*n1*/,
		const unsigned /*b2*/, const unsigned /*n2*/) const {
    return 1.;
  }

private:
  dcomp (*kernel_)(unsigned , unsigned );
  const vec3d* vectors_;
  const unsigned* op_perm_;
};

double eps_matgen;  //!< relative accuracy
double eps_aggl;    //!< relative accuracy
double eps_gmres;   //!< relative accuracy
unsigned steps_gmres;
double eta;   //!< control parameter for admissibility condition
unsigned bmin; //!< minimal size of a cluster
unsigned rankmax;  //!< maximal blockwise rank

unsigned nVectors;
vec3d* vectors;
dcomp *solution;

static PyObject *py_kernel = NULL;
dcomp (*kernel)(unsigned i, unsigned j);


static PyObject *
set_kernel(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  PyObject *temp;

  if (PyArg_ParseTuple(args, "O:set_kernel", &temp)) {
    if (!PyCallable_Check(temp)) {
      PyErr_SetString(PyExc_TypeError, "parameter must be callable");
      return NULL;
    }
    Py_XINCREF(temp);         /* Add a reference to new callback */
    Py_XDECREF(py_kernel);  /* Dispose of previous callback */
    py_kernel = temp;       /* Remember new callback */
    /* Boilerplate to return "None" */
    Py_INCREF(Py_None);
    result = Py_None;
  }
  return result;
};

dcomp
callback_kernel(unsigned i, unsigned j)
{
  double C_result_imag,C_result_real;
  PyObject *arglist;
  PyObject *result;

  arglist = Py_BuildValue("(ii)", i, j);
  result = PyObject_CallObject(py_kernel, arglist);
  Py_DECREF(arglist);
  C_result_real = PyComplex_RealAsDouble(result);
  C_result_imag = PyComplex_ImagAsDouble(result);
  return dcomp(C_result_real, C_result_imag);
};


static PyObject *py_vectorb = NULL;
dcomp (*vectorb)(unsigned i);


static PyObject *
set_vectorb(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  PyObject *temp;

  if (PyArg_ParseTuple(args, "O:set_vectorb", &temp)) {
    if (!PyCallable_Check(temp)) {
      PyErr_SetString(PyExc_TypeError, "parameter must be callable");
      return NULL;
    }
    Py_XINCREF(temp);         /* Add a reference to new callback */
    Py_XDECREF(py_vectorb);  /* Dispose of previous callback */
    py_vectorb = temp;       /* Remember new callback */
    /* Boilerplate to return "None" */
    Py_INCREF(Py_None);
    result = Py_None;
  }
  return result;
};

dcomp
callback_vectorb(unsigned i)
{
  double C_result_imag,C_result_real;
  PyObject *arglist;
  PyObject *result;

  arglist = Py_BuildValue("(i)", i);
  result = PyObject_CallObject(py_vectorb, arglist);
  Py_DECREF(arglist);
  C_result_real = PyComplex_RealAsDouble(result);
  C_result_imag = PyComplex_ImagAsDouble(result);
  return dcomp(C_result_real, C_result_imag);
};


static PyObject *
set_Hmatrix(PyObject *dummy, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "dddidii", &eps_matgen, &eps_aggl, &eps_gmres, &steps_gmres, &eta, &bmin, &rankmax))
    return NULL;
  return Py_None;
}

static PyObject *
set_points(PyObject *dummy, PyObject *args)
{
  PyArrayObject *points;
  NpyIter *points_iter;
  NpyIter_IterNextFunc *points_iternext;
  double ** points_dataptr;
  unsigned i;

  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &points))
    return NULL;
  points_iter = NpyIter_New(points, NPY_ITER_READONLY, NPY_KEEPORDER,
			NPY_NO_CASTING, NULL);
  if (points_iter == NULL)
    goto fail;
  points_iternext = NpyIter_GetIterNext(points_iter, NULL);
  if (points_iternext == NULL) {
    NpyIter_Deallocate(points_iter);
    goto fail;
  }

  points_dataptr = (double **) NpyIter_GetDataPtrArray(points_iter);

  nVectors = points->dimensions[0];
  vectors = new vec3d[nVectors];

  i = 0;
  do {
    vectors[i][0] = **points_dataptr;
    points_iternext(points_iter);
    vectors[i][1] = **points_dataptr;
    points_iternext(points_iter);
    vectors[i][2] = **points_dataptr;
    i++;
  } while(points_iternext(points_iter));


  return Py_None;

 fail:
  return NULL;
}

static PyObject *
get_q(PyObject *dummy, PyObject *args)
{
  PyArrayObject *q;
  NpyIter *q_iter;
  NpyIter_IterNextFunc *q_iternext;
  dcomp ** q_dataptr;
  unsigned i;

  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &q))
    return NULL;
  q_iter = NpyIter_New(q, NPY_ITER_READONLY, NPY_KEEPORDER,
			NPY_NO_CASTING, NULL);
  if (q_iter == NULL)
    goto fail;
  q_iternext = NpyIter_GetIterNext(q_iter, NULL);
  if (q_iternext == NULL) {
    NpyIter_Deallocate(q_iter);
    goto fail;
  }

  q_dataptr = (dcomp **) NpyIter_GetDataPtrArray(q_iter);

  nVectors = q->dimensions[0];
  vectors = new vec3d[nVectors];

  i = 0;
  do {
    **q_dataptr = solution[i];
    i++;
  } while(q_iternext(q_iter));

  return Py_None;

 fail:
  return NULL;
}

static PyObject *
solve_slae(PyObject *dummy, PyObject *args)
{
  Perm vecPerm(nVectors);

  bemcluster<vec3d>* clTreeVec =
    new bemcluster<vec3d>(vectors, vecPerm.op_perm, 0, nVectors);

  clTreeVec->createClusterTree(bmin, vecPerm.op_perm, vecPerm.po_perm);
  const unsigned nClustersPan = clTreeVec->getncl();

  std::cout << "done, " << nClustersPan << " clusters -- ";
  std::cout << inMB(nClustersPan * sizeof(bemcluster<vec3d>)) << " MB." << std::endl;

  bemblcluster<vec3d,vec3d>* blclTreeVec =
    new bemblcluster<vec3d,vec3d>(clTreeVec, clTreeVec);

  unsigned nblcksVec = 0;
  blclTreeVec->subdivide(clTreeVec,clTreeVec, eta * eta, nblcksVec);
  std::cout << "done, " << nblcksVec << " blocks -- ";
  std::cout << inMB(blclTreeVec->size()) << " MB." << std::endl;


  RealTimer timer;
  kernel = &callback_kernel;
  MATGENKERNEL MatGen(vectors, vecPerm.op_perm, kernel);
  BEMMatrix<vec3d,vec3d> A(nVectors, blclTreeVec);

  matgenGeH_omp(MatGen, nblcksVec, A.blclTree, eps_matgen, rankmax, A.blcks);
  {
    const double allmem = sizeH(A.blclTree, A.blcks);
    io::displayInfo(allmem, nVectors, timer.current(), sizeof(dcomp));
  }

  // std::cout << "Agglomerating matrix ... " << std::flush;
  // timer.restart();

  // agglH(A.blclTree, A.blcks, eps_aggl, rankmax);
  // std::cout << "done." << std::endl;

  // {
  //   const double allmem = sizeH(A.blclTree, A.blcks);
  //   io::displayInfo(allmem, nVectors, timer.current(), sizeof(dcomp));
  // }

  dcomp *b = new dcomp[nVectors];

  unsigned i = 0;
  do {
    b[i] = callback_vectorb(vecPerm.op_perm[i] + 1);
    i++;
  } while( i < nVectors );

  dcomp *x = new dcomp[nVectors];
  std::fill_n(x, nVectors, dcomp(0.,0.));

  double acc = eps_gmres;
  unsigned steps = steps_gmres;

  if (GMRes(A, b, x, acc, 100, steps))
    std::cout << "GMRes: iteration did not converge.";
  else
    std::cout << "GMRes converged to " << acc << " in "
  	      << steps << " steps.";

  std::cout << " Solution took " << timer << "." << std::endl;

  solution = new dcomp[nVectors];
  i = 0;
  do {
    solution[i] = x[vecPerm.po_perm[i]];
    i++;
  } while( i < nVectors );

  delete clTreeVec;

  return Py_None;
}

static PyObject *
internal(PyObject *dummy, PyObject *args)
{
  PyArrayObject *in_array;
  NpyIter *in_iter;
  NpyIter_IterNextFunc *in_iternext;
  double ** in_dataptr;

  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &in_array))
    return NULL;
  in_iter = NpyIter_New(in_array, NPY_ITER_READONLY, NPY_KEEPORDER,
  			NPY_NO_CASTING, NULL);
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  in_dataptr = (double **) NpyIter_GetDataPtrArray(in_iter);

  **in_dataptr = eps_matgen;
  in_iternext(in_iter);
  **in_dataptr = eps_aggl;
  in_iternext(in_iter);
  **in_dataptr = eps_gmres;
  in_iternext(in_iter);
  **in_dataptr = steps_gmres;
  in_iternext(in_iter);
  **in_dataptr = eta;
  in_iternext(in_iter);
  **in_dataptr = bmin;
  in_iternext(in_iter);
  **in_dataptr = rankmax;
  in_iternext(in_iter);
  NpyIter_Deallocate(in_iter);

  return Py_None;
}

char slaeahmed_set_Hmatrix_doc[] = "...";
char slaeahmed_set_points_doc[] = "...";
char slaeahmed_get_q_doc[] = "...";
char slaeahmed_internal_doc[] = "...";
char slaeahmed_solve_slae_doc[] = "...";
char slaeahmed_set_kernel_doc[] = "...";
char slaeahmed_set_vectorb_doc[] = "...";

static PyMethodDef slaeahmed_methods[] = {
  {"set_kernel",
   set_kernel,
   METH_VARARGS,
   slaeahmed_set_kernel_doc},
  {"set_vectorb",
   set_vectorb,
   METH_VARARGS,
   slaeahmed_set_vectorb_doc},
  {"set_points",
   set_points,
   METH_VARARGS,
   slaeahmed_set_points_doc},
  {"get_q",
   get_q,
   METH_VARARGS,
   slaeahmed_get_q_doc},
  {"set_Hmatrix",
   set_Hmatrix,
   METH_VARARGS,
   slaeahmed_set_Hmatrix_doc},
  {"internal",
   internal,
   METH_VARARGS,
   slaeahmed_internal_doc},
  {"solve_slae",
   solve_slae,
   METH_VARARGS,
   slaeahmed_solve_slae_doc},
  {NULL, NULL}
};

char slaeahmed_doc[] = "...";

PyMODINIT_FUNC initslaeahmed()
{
  Py_InitModule3("slaeahmed", slaeahmed_methods, slaeahmed_doc);
  import_array();
}
