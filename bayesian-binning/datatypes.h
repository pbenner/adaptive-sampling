
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
        int verbose;
        int likelihood;
        unsigned int sigma;
        unsigned int gamma;
} Options;
