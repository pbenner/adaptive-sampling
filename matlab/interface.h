/* Copyright (C) 2012 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <mex.h>

#include <adaptive-sampling/probtype.h>
#include <adaptive-sampling/linalg.h>
#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/interface.h>

void invalidOptions(const char* name);
double getScalar(const mxArray* array, const char* name);
mxArray* copyMatrixToMatlab(matrix_t* in);
mxArray* copyVectorToMatlab(vector_t* in);
mxArray* copyArrayToMatlab(prob_t* in, size_t size);
options_t* getOptions(const mxArray *array);
void copyMatrix(matrix_t* to, const mxArray* from);
void copy3DMatrix(matrix_t* to, const mxArray* from, size_t which);
void copyVector(vector_t* to, const mxArray* from);
void copyArray(prob_t* to, const mxArray* from, size_t size);
matrix_t** getCounts(const mxArray* array, size_t K, size_t L);
void freeCounts(matrix_t** counts);
matrix_t** getAlpha(const mxArray* array, size_t K, size_t L);
void freeAlpha(matrix_t** alpha);
vector_t* getBeta(const mxArray* array, size_t L);
matrix_t* getGamma(const mxArray* array, size_t L);
