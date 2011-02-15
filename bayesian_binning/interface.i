%module interface
%{
#include "interface.h"

PyObject * matrixToPyList(Matrix *m) {
    PyObject* mat = PyList_New(m->rows);
    int i, j;
    for(i = 0; i < m->rows; i++) {
       PyObject* item = PyList_New(m->columns);
       PyList_SetItem(mat, i, item);
       for(j = 0; j < m->columns; j++) {
          PyList_SetItem(item, j, PyFloat_FromDouble(m->mat[i][j]));
       }
    }
    return mat;
}

PyObject * vectorToPyList(Vector *v) {
    PyObject* vec = PyList_New(v->size);
    int i;
    for(i = 0; i < v->size; i++) {
       PyList_SetItem(vec, i, PyFloat_FromDouble(v->vec[i]));
    }
    return vec;
}
%}


%typemap(out) BinningResult * {
    PyObject * dict    = PyDict_New();

    if ($1->moments != NULL) {
            PyDict_SetItemString(dict, "moments", matrixToPyList($1->moments));
    }
    else {
            PyDict_SetItemString(dict, "moments", PyList_New(0));
    }
    if ($1->marginals != NULL) {
            PyDict_SetItemString(dict, "marginals", matrixToPyList($1->marginals));
    }
    else {
            PyDict_SetItemString(dict, "marginals", PyList_New(0));
    }
    PyDict_SetItemString(dict, "bprob",   vectorToPyList($1->bprob));
    PyDict_SetItemString(dict, "mpost",   vectorToPyList($1->mpost));
    PyDict_SetItemString(dict, "differential_gain", vectorToPyList($1->differential_gain));
    PyDict_SetItemString(dict, "effective_counts", vectorToPyList($1->effective_counts));
    PyDict_SetItemString(dict, "multibin_entropy", vectorToPyList($1->multibin_entropy));

    $result = dict;
}

%typemap(freearg) BinningResult * {
    freeMatrix($1->moments);
    freeMatrix($1->marginals);
    freeVector($1->bprob);
    freeVector($1->mpost);
    freeVector($1->differential_gain);
    freeVector($1->effective_counts);
    freeVector($1->multibin_entropy);
    free($1);
}

%typemap(in) Options * {
    $1 = (Options *)malloc(sizeof(Options));
    PyObject *verbose        = PyDict_GetItemString($input, "verbose");
    PyObject *prombsTest     = PyDict_GetItemString($input, "prombsTest");
    PyObject *epsilon        = PyDict_GetItemString($input, "epsilon");
    PyObject *gmp            = PyDict_GetItemString($input, "gmp");
    PyObject *bprob          = PyDict_GetItemString($input, "bprob");
    PyObject *which          = PyDict_GetItemString($input, "which");
    PyObject *marginal       = PyDict_GetItemString($input, "marginal");
    PyObject *marginal_step  = PyDict_GetItemString($input, "marginal_step");
    PyObject *marginal_range = PyDict_GetItemString($input, "marginal_range");
    PyObject *n_moments      = PyDict_GetItemString($input, "n_moments");
    PyObject *differential_gain = PyDict_GetItemString($input, "differential_gain");
    PyObject *effective_counts  = PyDict_GetItemString($input, "effective_counts");
    PyObject *multibin_entropy  = PyDict_GetItemString($input, "multibin_entropy");
    PyObject *model_posterior   = PyDict_GetItemString($input, "model_posterior");
    $1->verbose    = (verbose    == Py_True ? 1 : 0);
    $1->prombsTest = (prombsTest == Py_True ? 1 : 0);
    $1->gmp        = (gmp        == Py_True ? 1 : 0);
    $1->bprob      = (bprob      == Py_True ? 1 : 0);
    $1->differential_gain = (differential_gain  == Py_True ? 1 : 0);
    $1->effective_counts  = (effective_counts   == Py_True ? 1 : 0);
    $1->multibin_entropy  = (multibin_entropy   == Py_True ? 1 : 0);
    $1->model_posterior   = (model_posterior    == Py_True ? 1 : 0);
    $1->which         = PyInt_AsLong(which);
    $1->marginal      = PyInt_AsLong(marginal);
    $1->marginal_step = PyFloat_AsDouble(marginal_step);
    $1->marginal_range.from = PyFloat_AsDouble(PyTuple_GetItem(marginal_range, 0));
    $1->marginal_range.to   = PyFloat_AsDouble(PyTuple_GetItem(marginal_range, 1));
    $1->n_marginals   = (int)floor(1.0/$1->marginal_step) + 1;
    $1->n_moments     = PyInt_AsLong(n_moments);
    $1->epsilon       = PyFloat_AsDouble(epsilon);
}

%typemap(freearg) Options * {
    free($1);
}

%typemap(in) Vector * {
    if (!PyList_Check($input)) {
        PyErr_SetString(PyExc_TypeError, "Not a list.");
        return NULL;
    }
    int size = PyList_Size($input);
    $1 = allocVector(size);
    int i;
    for (i = 0; i < size; i++) {
        PyObject *item = PyList_GetItem($input, i);
        if (PyFloat_Check(item)) {
	   $1->vec[i] = PyFloat_AsDouble(item);
        }
        else
	if (PyInt_Check(item)) {
	   $1->vec[i] = (double)PyInt_AsLong(item);
	}
	else {
            PyErr_SetString(PyExc_TypeError, "List item is not a float or integer.");
	    freeVector($1);
            return NULL;
	}
    }
}

%typemap(out) Vector * {
    $result = vectorToPyList($1);
}

%typemap(freearg) Vector * {
    freeVector((Vector *)$1);
}

%typemap(in) Matrix * {
    if (!PyList_Check($input)) {
        PyErr_SetString(PyExc_TypeError, "Not a matrix.");
        return NULL;
    }
    int rows    = PyList_Size($input);
    int columns = PyList_Size(PyList_GetItem($input, 0));
    int i, j;
    for (i = 0; i < rows; i++) {
        PyObject *item = PyList_GetItem($input, i);
	if (!PyList_Check(item) || PyList_Size(item) != columns) {
            PyErr_SetString(PyExc_TypeError, "Not a matrix.");
       	    return NULL;
    	}
	for (j = 0; j < columns; j++) {
            PyObject *cell = PyList_GetItem(item, j);
            if (!(PyFloat_Check(cell) || PyInt_Check(cell))) {
               PyErr_SetString(PyExc_TypeError, "Matrix item is not a float or integer.");
               return NULL;
            }
	}
    }
    $1 = allocMatrix(rows, columns);
    for (i = 0; i < rows; i++) {
        PyObject *item = PyList_GetItem($input, i);
	for (j = 0; j < columns; j++) {
            PyObject *cell = PyList_GetItem(item, j);
            if (PyFloat_Check(cell)) {
               $1->mat[i][j] = PyFloat_AsDouble(cell);
	    }
	    else {
               $1->mat[i][j] = (double)PyInt_AsLong(cell);
	    }
	}
    }
}

%typemap(out) Matrix * {
    $result = matrixToPyList($1);
}

%typemap(freearg) Matrix * {
    freeMatrix((Matrix *)$1);
}

%include interface.c
