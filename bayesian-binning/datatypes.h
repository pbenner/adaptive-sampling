
typedef struct vector {
        int size;
        double *vec;
} Vector;

typedef struct matrix {
        int rows;
        int columns;
        double **mat;
} Matrix;

typedef struct options {
        long double epsilon;
        int verbose;
        int gmp;
        int likelihood;
        int bprob;
        int which;
} Options;
