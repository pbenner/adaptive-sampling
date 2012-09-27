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
int getbool(SEXP env, const char* name, size_t n)
{
        SEXP ans;
        int result;

        if(!isEnvironment(env)) {
                error("should be an environment");
        }
        PROTECT(ans = findVar(install(name), env));
        if(isSymbol(ans) && isNull(TAG(ans))) {
                error("options entry for %s is null", name);
        }
        result = LOGICAL(ans)[n];
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
options_t* getOptions(SEXP r_options)
{
        options_t* options = (options_t*)malloc(sizeof(options_t));

        options->model_posterior            = getbool(r_options, "model.posterior", 0);
        options->kl_psi                     = getbool(r_options, "kl.psi", 0);
        options->kl_multibin                = getbool(r_options, "kl.multibin",  0);
        options->effective_counts           = getbool(r_options, "effective.counts", 0);
        options->effective_posterior_counts = getbool(r_options, "effective.posterior.counts", 0);
        options->bprob                      = getbool(r_options, "bprob", 0);
        options->density                    = getbool(r_options, "density", 0);
        options->density_step               = getreal(r_options, "density.step", 0);
        options->density_range.from         = getreal(r_options, "density.range", 0);
        options->density_range.to           = getreal(r_options, "density.range", 1);
        options->n_moments                  = getreal(r_options, "n.moments", 0);
        options->n_density                  = floor(1.0/options->density_step) + 1;
        options->epsilon                    = getreal(r_options, "epsilon", 0);
        options->verbose                    = 0;
        options->prombsTest                 = 0;
        options->threads                    = getreal(r_options, "threads", 0);
        options->stacksize                  = getreal(r_options, "stacksize", 0);
        options->algorithm                  = getreal(r_options, "algorithm", 0);
        options->which                      = getreal(r_options, "which", 0);
        options->samples[0]                 = getreal(r_options, "samples", 0);
        options->samples[1]                 = getreal(r_options, "samples", 1);
        options->hmm                        = getbool(r_options, "hmm", 0);
        options->rho                        = getreal(r_options, "rho", 0);

        return options;
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

/******************************************************************************
 * posterior
 *****************************************************************************/

static
SEXP copyPosterior(marginal_t* result) {
        SEXP r_matrix;
        SEXP r_vector;
        SEXP r_result;

        PROTECT(r_result = allocSExp(ENVSXP));

        if (result->moments) {
                PROTECT(r_matrix = copyMatrixToR(result->moments));
                defineVar(install("moments"), r_matrix, r_result);
                UNPROTECT(1);
        }
        if (result->density) {
                PROTECT(r_matrix = copyMatrixToR(result->density));
                defineVar(install("density"), r_matrix, r_result);
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
        UNPROTECT(1);

        return r_result;
}

static
void freePosterior(marginal_t * result) {
        if (result->moments) {
                free_matrix(result->moments);
        }
        if (result->density) {
                free_matrix(result->density);
        }
        if (result->bprob) {
                free_vector(result->bprob);
        }
        if (result->mpost) {
                free_vector(result->mpost);
        }
        free(result);
}

static
marginal_t * callPosterior(
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
        options_t* options;
        marginal_t * result;

        counts  = getCounts(L, K, r_counts);
        alpha   = getCounts(L, K, r_alpha);
        beta    = getBeta(r_beta, L);
        gamma   = getGamma(r_gamma, L);
        options = getOptions(r_options);

        result  = posterior(K, counts, alpha, beta, gamma, options);

        freeCounts(counts);
        freeCounts(alpha);
        free_vector(beta);
        free_matrix(gamma);
        free(options);

        return result;
}

SEXP call_posterior(
        SEXP r_counts,
        SEXP r_alpha,
        SEXP r_beta,
        SEXP r_gamma,
        SEXP r_options)
{
        marginal_t *result;
        SEXP r_result;

        check_input(r_counts, r_alpha, r_beta, r_gamma, r_options);

        result = callPosterior(r_counts, r_alpha, r_beta, r_gamma, r_options);
        PROTECT(r_result = copyPosterior(result));
        freePosterior(result);
        UNPROTECT(1);

        return r_result;
}

/******************************************************************************
 * utility
 *****************************************************************************/

static
utility_t* callUtility(
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
        options_t* options;
        utility_t * result;

        counts  = getCounts(L, K, r_counts);
        alpha   = getCounts(L, K, r_alpha);
        beta    = getBeta(r_beta, L);
        gamma   = getGamma(r_gamma, L);
        options = getOptions(r_options);

        result  = utility(K, counts, alpha, beta, gamma, options);

        freeCounts(counts);
        freeCounts(alpha);
        free_vector(beta);
        free_matrix(gamma);
        free(options);

        return result;
}

static
SEXP copyUtility(utility_t* result) {
        SEXP r_matrix;
        SEXP r_vector;
        SEXP r_result;

        PROTECT(r_result = allocSExp(ENVSXP));

        if (result->expectation) {
                PROTECT(r_matrix = copyMatrixToR(result->expectation));
                defineVar(install("expectation"), r_matrix, r_result);
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
void freeUtility(utility_t * result) {
        if (result->expectation) {
                free_matrix(result->expectation);
        }
        if (result->utility) {
                free_vector(result->utility);
        }
}

SEXP call_utility(
        SEXP r_counts,
        SEXP r_alpha,
        SEXP r_beta,
        SEXP r_gamma,
        SEXP r_options)
{
        utility_t* result;
        SEXP r_result;
        SEXP r_matrix;

        check_input(r_counts, r_alpha, r_beta, r_gamma, r_options);

        result = callUtility(r_counts, r_alpha, r_beta, r_gamma, r_options);
        PROTECT(r_result = copyUtility(result));
        freeUtility(result);
        UNPROTECT(1);

        return r_result;
}
