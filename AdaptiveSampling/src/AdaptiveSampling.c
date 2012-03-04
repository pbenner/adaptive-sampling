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

#include <adaptive-sampling/linalg.h>
#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/interface.h>

static
double getreal(SEXP env, const char* name, size_t n)
{
        SEXP ans;
        double result;

        if(!isEnvironment(env)) {
                error("should be an environment");
        }
        PROTECT(ans = findVar(install(name), env));
        if(isSymbol(ans) && isNull(TAG(ans))) {
                error("options entry for %s is null", name);
        }
        result = REAL(ans)[n];
        UNPROTECT(1);

        return result;
}

static
void copyMatrixToC(matrix_t* to, SEXP from, size_t rows, size_t columns, size_t offset) {
        double *r_counts = REAL(from);
        size_t i, j;

        for (i = 0; i < rows; i++) {
                for (j = 0; j < columns; j++) {
                        to->content[i][j] = r_counts[j*columns+i+offset];
                }
        }
}

static
void copyVectorToC(vector_t* to, SEXP from, size_t length) {
        double *r_counts = REAL(from);
        size_t i;

        for (i = 0; i < length; i++) {
                to->content[i] = r_counts[i];
        }
}

static
SEXP copyMatrixToR(matrix_t* from) {
        SEXP r_matrix;
        PROTECT(r_matrix = allocMatrix(REALSXP, from->rows, from->columns));
        double *rp_matrix = REAL(r_matrix);
        size_t i, j;

        for (i = 0; i < from->rows; i++) {
                for (j = 0; j < from->columns; j++) {
                        rp_matrix[j*from->rows+i] = from->content[i][j];
                }
        }
        UNPROTECT(1);
        return r_matrix;
}

static
SEXP copyVectorToR(vector_t* from) {
        SEXP r_vector;
        PROTECT(r_vector = allocVector(REALSXP, from->size));
        double *rp_vector = REAL(r_vector);
        size_t i;

        for (i = 0; i < from->size; i++) {
                rp_vector[i] = from->content[i];
        }
        UNPROTECT(1);
        return r_vector;
}

static
matrix_t** getCounts(size_t L, size_t K, SEXP r_counts) {
        matrix_t** counts = (matrix_t**)malloc((K+1)*sizeof(matrix_t*));
        size_t k;

        for (k = 0; k < K; k++) {
                counts[k] = alloc_matrix(L, L);
                copyMatrixToC(counts[k], r_counts, L, L, k*L*L);
        }
        counts[k] = NULL;

        return counts;
}

static
void freeCounts(matrix_t** counts) {
        size_t i;

        for (i = 0; counts[i]; i++) {
                free_matrix(counts[i]);
        }
        free(counts);
}

static
vector_t* getBeta(SEXP r_beta, size_t L) {
        vector_t* beta = alloc_vector(L);

        copyVectorToC(beta, r_beta, L);

        return beta;
}

static
matrix_t* getGamma(SEXP r_gamma, size_t L) {
        matrix_t* gamma = alloc_matrix(L, L);

        copyMatrixToC(gamma, r_gamma, L, L, 0);

        return gamma;
}

static
Options* getOptions(SEXP r_options)
{
        Options* options = (Options*)malloc(sizeof(Options));

        options->model_posterior     = getreal(r_options, "model_posterior", 0);
        options->utility             = getreal(r_options, "utility", 0);
        options->kl_component        = getreal(r_options, "kl_component", 0);
        options->kl_multibin         = getreal(r_options, "kl_multibin",  0);
        options->effective_counts    = getreal(r_options, "effective_counts", 0);
        options->bprob               = getreal(r_options, "bprob", 0);
        options->marginal            = getreal(r_options, "marginal", 0);
        options->marginal_step       = getreal(r_options, "marginal_step", 0);
        options->marginal_range.from = getreal(r_options, "marginal_range", 0);
        options->marginal_range.to   = getreal(r_options, "marginal_range", 1);
        options->n_moments           = getreal(r_options, "n_moments", 0);
        options->n_marginals         = floor(1.0/options->marginal_step) + 1;
        options->epsilon             = getreal(r_options, "epsilon", 0);
        options->verbose             = 0;
        options->prombsTest          = 0;
        options->threads             = getreal(r_options, "threads", 0);
        options->stacksize           = getreal(r_options, "stacksize", 0);
        options->algorithm           = getreal(r_options, "algorithm", 0);
        options->which               = getreal(r_options, "which", 0);
        options->samples[0]          = getreal(r_options, "samples", 0);
        options->samples[1]          = getreal(r_options, "samples", 1);

        return options;
}

static
SEXP copyResult(BinningResult* result) {
        SEXP r_matrix;
        SEXP r_vector;
        SEXP r_result;

        PROTECT(r_result = allocSExp(ENVSXP));

        if (result->moments) {
                PROTECT(r_matrix = copyMatrixToR(result->moments));
                defineVar(install("moments"), r_matrix, r_result);
                UNPROTECT(1);
        }
        if (result->marginals) {
                PROTECT(r_matrix = copyMatrixToR(result->marginals));
                defineVar(install("marginals"), r_matrix, r_result);
                UNPROTECT(1);
        }
        if (result->bprob) {
                PROTECT(r_vector = copyVectorToR(result->bprob));
                defineVar(install("bprob"), r_vector, r_result);
                UNPROTECT(1);
        }
        if (result->mpost) {
                PROTECT(r_vector = copyVectorToR(result->mpost));
                defineVar(install("mpost"), r_vector, r_result);
                UNPROTECT(1);
        }
        if (result->utility) {
                PROTECT(r_vector = copyVectorToR(result->utility));
                defineVar(install("utility"), r_vector, r_result);
                UNPROTECT(1);
        }
        UNPROTECT(1);

        return r_result;
}

static
void freeResult(BinningResult * result) {
        if (result->moments) {
                free_matrix(result->moments);
        }
        if (result->marginals) {
                free_matrix(result->marginals);
        }
        if (result->bprob) {
                free_vector(result->bprob);
        }
        if (result->mpost) {
                free_vector(result->mpost);
        }
        if (result->utility) {
                free_vector(result->utility);
        }
        free(result);
}

static
BinningResult * callBinning(
        SEXP r_counts,
        SEXP r_alpha,
        SEXP r_beta,
        SEXP r_gamma,
        SEXP r_options) {

        SEXP dim = getAttrib(r_counts, R_DimSymbol);
        size_t L = INTEGER(dim)[0];
        size_t K = INTEGER(dim)[2];
        matrix_t** counts;
        matrix_t** alpha;
        vector_t* beta;
        matrix_t* gamma;
        Options* options;
        BinningResult * result;

        counts  = getCounts(L, K, r_counts);
        alpha   = getCounts(L, K, r_alpha);
        beta    = getBeta(r_beta, L);
        gamma   = getGamma(r_gamma, L);
        options = getOptions(r_options);

        result = binning(K, counts, alpha, beta, gamma, options);

        freeCounts(counts);
        freeCounts(alpha);
        free_vector(beta);
        free_matrix(gamma);
        free(options);

        return result;
}

SEXP check_input(
        SEXP r_counts,
        SEXP r_alpha,
        SEXP r_beta,
        SEXP r_gamma,
        SEXP r_options)
{
        /* check counts */
        SEXP dim = getAttrib(r_counts, R_DimSymbol);
        if(length(dim) != 3) {
                error("counts has invalid dimension");
        }
        size_t L = INTEGER(dim)[0];
        size_t K = INTEGER(dim)[2];
        /* check alpha */
        dim = getAttrib(r_alpha, R_DimSymbol);
        if(length(dim) != 3 || INTEGER(dim)[0] != L || INTEGER(dim)[2] != K) {
                error("alpha has invalid dimension");
        }
        /* check beta */
        if(length(r_beta) != L) {
                error("beta has invalid dimension");
        }
        /* check gamma */
        dim = getAttrib(r_gamma, R_DimSymbol);
        if(length(dim) != 2 || INTEGER(dim)[0] != L || INTEGER(dim)[1] != L) {
                error("gamma has invalid dimension");
        }
        /* check r_options */
        if(!isEnvironment(r_options)) {
                error("options should be an environment");
        }

        return R_NilValue;
}

SEXP adaptive_sampling(
        SEXP r_counts,
        SEXP r_alpha,
        SEXP r_beta,
        SEXP r_gamma,
        SEXP r_options)
{
        BinningResult *result;
        SEXP r_result;

        check_input(r_counts, r_alpha, r_beta, r_gamma, r_options);

        result = callBinning(r_counts, r_alpha, r_beta, r_gamma, r_options);
        PROTECT(r_result = copyResult(result));
        freeResult(result);
        UNPROTECT(1);

        return r_result;
}
