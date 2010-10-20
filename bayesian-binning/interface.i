%module interface
%{
#include "interface.h"
%}

%typemap(in) Options * {
    $1 = (Options *)malloc(sizeof(Options));
    PyObject *verbose    = PyDict_GetItemString($input, "verbose");
    PyObject *gmp        = PyDict_GetItemString($input, "gmp");
    PyObject *bprob      = PyDict_GetItemString($input, "bprob");
    PyObject *likelihood = PyDict_GetItemString($input, "likelihood");
    PyObject *which      = PyDict_GetItemString($input, "which");
    $1->verbose    = (verbose == Py_True ? 1 : 0);
    $1->gmp        = (gmp     == Py_True ? 1 : 0);
    $1->bprob      = (bprob   == Py_True ? 1 : 0);
    $1->likelihood = PyInt_AsLong(likelihood);
    $1->which      = PyInt_AsLong(which);
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
    PyObject* list = PyList_New($1->size);
    int i;
    for(i = 0; i < $1->size; i++) {
       PyList_SetItem(list, i, PyFloat_FromDouble($1->vec[i]));
    }
    $result = list;
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
    PyObject* mat = PyList_New($1->rows);
    int i, j;
    for(i = 0; i < $1->rows; i++) {
       PyObject* item = PyList_New($1->columns);
       PyList_SetItem(mat, i, item);
       for(j = 0; j < $1->columns; j++) {
          PyList_SetItem(item, j, PyFloat_FromDouble($1->mat[i][j]));
       }
    }
    $result = mat;
}

%typemap(freearg) Matrix * {
    freeMatrix((Matrix *)$1);
}

%include interface.c
