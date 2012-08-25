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

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <adaptive-sampling/prombs.h>

static
void copyMatrixToC(prob_t** to, SEXP from, size_t L) {
        double *r_counts = REAL(from);
        size_t i, j;

        for (i = 0; i < L; i++) {
                for (j = 0; j < L; j++) {
                        to[i][j] = r_counts[j*L+i];
                }
        }
}

static
void copyVectorToC(prob_t* to, SEXP from, size_t L) {
        double *r_counts = REAL(from);
        size_t i;

        for (i = 0; i < L; i++) {
                to[i] = r_counts[i];
        }
}

static
SEXP copyMatrixToR(prob_t** from, size_t L) {
        SEXP r_matrix;
        PROTECT(r_matrix = allocMatrix(REALSXP, L, L));
        double *rp_matrix = REAL(r_matrix);
        size_t i, j;

        for (i = 0; i < L; i++) {
                for (j = 0; j < L; j++) {
                        rp_matrix[j*L+i] = from[i][j];
                }
        }
        UNPROTECT(1);
        return r_matrix;
}

static
SEXP copyVectorToR(prob_t* from, size_t L) {
        SEXP r_vector;
        PROTECT(r_vector = allocVector(REALSXP, L));
        double *rp_vector = REAL(r_vector);
        size_t i;

        for (i = 0; i < L; i++) {
                rp_vector[i] = from[i];
        }
        UNPROTECT(1);
        return r_vector;
}

SEXP check_input(
        SEXP r_g,
        SEXP r_f,
        SEXP r_h,
        SEXP r_epsilon,
        SEXP r_m)
{
        /* check g */
        SEXP dim = getAttrib(r_g, R_DimSymbol);
        if(length(dim) != 1) {
                printf("dim: %d\n", (int)length(dim));
                error("g is not a vector");
        }
        size_t L = INTEGER(dim)[0];
        /* check f */
        dim = getAttrib(r_f, R_DimSymbol);
        if(length(dim) != 2 || INTEGER(dim)[0] != L || INTEGER(dim)[1] != L) {
                error("f has invalid dimension");
        }
        /* check h */
        if (r_h) {
                dim = getAttrib(r_h, R_DimSymbol);
                if(length(dim) != 2 || INTEGER(dim)[0] != L || INTEGER(dim)[1] != L) {
                        error("h has invalid dimension");
                }
        }
        /* check m */
        dim = getAttrib(r_m, R_DimSymbol);
        if(length(dim) != 1 || INTEGER(dim)[0] != 1) {
                error("m is not a scalar");
        }
        if (INTEGER(r_m)[0] == 0 || INTEGER(r_m)[0] > L) {
                error("m has invalid value");
        }
        /* check epsilon */
        if (r_epsilon) {
                dim = getAttrib(r_epsilon, R_DimSymbol);
                if(length(dim) != 1 || INTEGER(dim)[0] != 1) {
                        error("epsilon is not a scalar");
                }
        }

        return R_NilValue;
}


/******************************************************************************
 * tools
 *****************************************************************************/

static __inline__
prob_t** alloc_prombs_matrix(size_t L) {
        size_t i;
        prob_t** m = (prob_t **)malloc(L*sizeof(prob_t*));
        for (i = 0; i < L; i++) {
                m[i] = (prob_t *)malloc(L*sizeof(prob_t));
        }
        return m;
}

static __inline__
void free_prombs_matrix(prob_t** m, size_t L) {
        size_t i;
        for (i = 0; i < L; i++) {
                free(m[i]);
        }
        free(m);
}

/******************************************************************************
 * prombs interface
 *****************************************************************************/

SEXP call_prombs(
        SEXP r_g,
        SEXP r_f,
        SEXP r_m)
{
        SEXP r_result;

        /* check whether the arguments are ok */
        check_input(r_g, r_f, NULL, NULL, r_m);

        /* dimension */
        SEXP dim = getAttrib(r_g, R_DimSymbol);
        size_t L = INTEGER(dim)[0];
        size_t m = INTEGER(r_m)[0] - 1;

        /* copy arguments */
        prob_t* result = (prob_t *)malloc(L*sizeof(prob_t));
        prob_t*  g     = (prob_t *)malloc(L*sizeof(prob_t));
        prob_t** f     = alloc_prombs_matrix(L);
        copyVectorToC(g, r_g, L);
        copyMatrixToC(f, r_f, L);

        prombs(result, f, g, NULL, L, m, NULL);

        PROTECT(r_result = copyVectorToR(result, L));

        free(result);
        free(g);
        free_prombs_matrix(f, L);

        UNPROTECT(1);
        return r_result;
}

/******************************************************************************
 * extended prombs interface
 *****************************************************************************/

typedef struct _data_t_ {
        const double* f;
        const double* h;
        const size_t L;
} data_t;

prob_t prombs_f(int i, int j, void* _data)
{
        data_t* data = (data_t *)_data;

        return data->f[j*data->L+i];
}

prob_t prombs_h(int i, int j, void* _data)
{
        data_t* data = (data_t *)_data;

        return data->h[j*data->L+i];
}

SEXP call_prombs_extended(
        SEXP r_g,
        SEXP r_f,
        SEXP r_h,
        SEXP r_epsilon,
        SEXP r_m)
{
        SEXP r_result;

        /* check whether the arguments are ok */
        check_input(r_g, r_f, r_h, r_epsilon, r_m);

        /* set precision */
        __init_prombs__(REAL(r_epsilon)[0]);

        /* dimension */
        SEXP dim = getAttrib(r_g, R_DimSymbol);
        size_t L = INTEGER(dim)[0];
        size_t m = INTEGER(r_m)[0] - 1;

        /* copy arguments */
        data_t  data   = { REAL(r_f), REAL(r_h), L };
        prob_t* result = (prob_t *)malloc(L*sizeof(prob_t));
        prob_t*  g     = (prob_t *)malloc(L*sizeof(prob_t));
        prob_t** ak    = alloc_prombs_matrix(L);
        copyVectorToC(g, r_g, L);

        prombsExt(result, ak, g, prombs_f, prombs_h, L, m, &data);

        PROTECT(r_result = copyVectorToR(result, L));

        free(result);
        free(g);
        free_prombs_matrix(ak, L);

        UNPROTECT(1);
        return r_result;
}
