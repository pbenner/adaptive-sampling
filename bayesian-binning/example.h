
struct vector {
        int size;
        double *vec;
};

struct matrix {
        int rows;
        int columns;
        double **mat;
};

typedef struct vector Vector;
typedef struct matrix Matrix;

extern int cube( int n );
extern int test( Matrix *m);
extern double product(Vector *l);
