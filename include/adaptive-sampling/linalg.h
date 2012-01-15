/* Copyright (C) 2010 Philipp Benner
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

#ifndef LINALG_H
#define LINALG_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stddef.h>

#include <gsl/gsl_matrix.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct {
        int size;
        double *content;
} vector_t;

typedef struct {
        int rows;
        int columns;
        double **content;
} matrix_t;

static __inline__
vector_t * alloc_vector(int size) {
        vector_t *v = (vector_t *)malloc(sizeof(vector_t));
        v->content  = (double   *)calloc(size, sizeof(double));
        v->size     = size;
        return v;
}

static __inline__
matrix_t * alloc_matrix(int rows, int columns) {
        matrix_t *m = (matrix_t *)malloc(sizeof(matrix_t));
        m->content  = (double  **)calloc(rows, sizeof(double *));
        m->rows     = rows;
        m->columns  = columns;
        int i;
        for (i = 0; i < rows; i++) {
                m->content[i] = (double *)calloc(columns, sizeof(double));
        }
        return m;
}

static __inline__
void free_vector(vector_t *v) {
        free(v->content);
        free(v);
}

static __inline__
void free_matrix(matrix_t *m) {
        int i;
        for (i = 0; i < m->rows; i++) {
                free(m->content[i]);
        }
        free(m->content);
        free(m);
}

static __inline__
gsl_vector * to_gsl_vector(const vector_t *vector)
{
        gsl_vector *v = gsl_vector_alloc(vector->size);
        int i;

        for(i = 0; i < vector->size; i++) {
                gsl_vector_set(v, i, vector->content[i]);
        }
        return v;
}

static __inline__
vector_t * from_gsl_vector(const gsl_vector * vector)
{
        int i;
        int size = vector->size;
        vector_t *v = alloc_vector(size);

        for(i = 0; i < size; i++) {
                v->content[i] = gsl_vector_get(vector, i);
        }
        return v;
}

static __inline__
gsl_matrix * to_gsl_matrix(const matrix_t *matrix)
{
        gsl_matrix *m = gsl_matrix_alloc(matrix->rows, matrix->columns);
        int i, j;

        for(i = 0; i < matrix->rows; i++) {
                for(j = 0; j < matrix->columns; j++) {
                        gsl_matrix_set(m, i, j, matrix->content[i][j]);
                }
        }
        return m;
}

static __inline__
matrix_t * from_gsl_matrix(const gsl_matrix * matrix)
{
        int i, j;
        int rows    = matrix->size1;
        int columns = matrix->size2;
        matrix_t *m    = alloc_matrix(rows, columns);

        for(i = 0; i < rows; i++) {
                for(j = 0; j < columns; j++) {
                        m->content[i][j] = gsl_matrix_get(matrix, i, j);
                }
        }
        return m;
}

static __inline__
void print_vector(vector_t* v) {
        int i;
        for (i = 0; i < v->size; i++) {
                printf("%f ", v->content[i]);
        }
        printf("\n");
}

static __inline__
void print_matrix(matrix_t* m) {
        int i, j;
        for (i = 0; i < m->rows; i++) {
                for (j = 0; j < m->columns; j++) {
                        printf("%f ", m->content[i][j]);
                }
                printf("\n");
        }
}

__END_DECLS

#endif /* LINALG_H */
