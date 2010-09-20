%module example
%{
#include "example.h"
%}

%typemap(in) Vector * {
    if (!PyList_Check($input)) {
        PyErr_SetString(PyExc_TypeError, "Not a list.");
        return NULL;
    }
    int size = PyList_Size($input);
    $1 = (Vector *)malloc(sizeof(Vector));
    double *vec = (double *)malloc(sizeof(double) * size);
    $1->vec  = vec;
    $1->size = size;
    int i;
    for (i = 0; i < size; i++) {
        PyObject *item = PyList_GetItem($input, i);
        if (PyFloat_Check(item)) {
	   vec[i] = PyFloat_AsDouble(item);
        }
        else
	if (PyInt_Check(item)) {
	   vec[i] = (double)PyInt_AsLong(item);
	}
	else {
            PyErr_SetString(PyExc_TypeError, "List item is not a float or integer.");
            free(vec);
            free($1);
            return NULL;
	}
    }
}

%typemap(in) Matrix * {
    if (!PyList_Check($input)) {
        PyErr_SetString(PyExc_TypeError, "Not a matrix.");
        return NULL;
    }
    int rows    = PyList_Size($input);
    int columns = PyList_Size(PyList_GetItem($input, 0));
    $1 = (Matrix *)malloc(sizeof(Matrix));
    double **mat = (double **)malloc(sizeof(double *) * rows);
    $1->mat      = mat;
    $1->rows     = rows;
    $1->columns  = columns;
    int i, j;
    for (i = 0; i < rows; i++) {
        PyObject *item = PyList_GetItem($input, i);
	if (!PyList_Check(item) || PyList_Size(item) != columns) {
            PyErr_SetString(PyExc_TypeError, "Not a matrix.");
            free(mat);
            free($1);
       	    return NULL;
    	}
	for (j = 0; j < columns; j++) {
            PyObject *cell = PyList_GetItem(item, j);
            if (!(PyFloat_Check(cell) || PyInt_Check(cell))) {
               PyErr_SetString(PyExc_TypeError, "Matrix item is not a float or integer.");
               free(mat);
               free($1);
               return NULL;
            }
	}
    }
    for (i = 0; i < rows; i++) {
        PyObject *item = PyList_GetItem($input, i);
	mat[i] = (double *)malloc(sizeof(double) * columns);
	for (j = 0; j < columns; j++) {
            PyObject *cell = PyList_GetItem(item, j);
            if (PyFloat_Check(cell)) {
               mat[i][j] = PyFloat_AsDouble(cell);
	    }
	    else {
               mat[i][j] = (double)PyInt_AsLong(cell);
	    }
	}
    }
}

%typemap(freearg) Vector * {
    free(((Vector *)$1)->vec);
    free((Vector *)$1);
}

%typemap(freearg) Matrix * {
    int i;
    for (i = 0; i < ((Matrix *)$1)->rows; i++) {
        free(((Matrix *)$1)->mat[i]);
    }
    free(((Matrix *)$1)->mat);
    free((Matrix *)$1);
}

%include example.c
