#ifndef HAVE_HEADER
#define HAVE_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <assert.h>

typedef struct {int size; double *data;} vector;
typedef struct {int size1, size2; double *data;} matrix;

vector* vector_alloc                (int n);
void    vector_free                 (vector* v);
void    vector_set                  (vector* v, int i, double x);
double  vector_get                  (vector* v, int i);
void    vector_set_random           (vector* v);
void    vector_set_zero             (vector* v);
void    vector_set_basis            (vector* v, int i);
void    vector_fprintf              (FILE* fp, vector* v);

void vector_memcpy (vector* original, vector* copy);
void matrix_memcpy (matrix* original, matrix* copy);
void matrix_column_view_vector (matrix* A, int j, vector* v);
void matrix_row_view_vector (matrix* A, int j, vector* v);
void vector_view_matrix_column (matrix* A, int j, vector* v);

void matrix_vector_product (matrix* A, vector* b, vector* Ab);
void matrix_matrix_product (matrix* A, matrix* B, matrix* AB);
void matrix_transpose_memcpy (matrix* src, matrix* dest);
void matrix_inverse_upper_triangular (matrix* R, matrix* R_inv);

matrix* matrix_alloc                (int n, int m);
void    matrix_free                 (matrix* A);
void    matrix_set                  (matrix* A, int i, int j, double x);
double  matrix_get                  (matrix* A, int i, int j);
void    matrix_set_zero             (matrix* A);
void    matrix_set_identity         (matrix* A);
void    matrix_scale                (matrix* A, double s);
void    matrix_set_random           (matrix* A);
void    matrix_set_random_symmetric (matrix* A);
void    matrix_fprintf              (FILE*   fp, matrix* A);
double  matrix_column_norm          (matrix* A,  int     j);

void    GS_decomp           (matrix* A, matrix* Q,  matrix* R);
void    GS_solve            (matrix* Q, matrix* R,  vector* b, vector* x);
void    GS_inverse          (matrix* Q, matrix* R, matrix* QRI);
void    backsub             (matrix* A,  vector* b);

int real_equal (double a, double b, double tau, double epsilon);
int vector_vector_test_equal (vector* v, vector* w, double epsabs, double epsrel);
int matrix_matrix_test_equal (matrix* A, matrix* B, double epsabs, double epsrel);
void fprintf_result_vector_vector_test_equal (FILE* fp, vector* v, vector* w, double epsabs, double epsrel);
void fprintf_result_matrix_matrix_test_equal (FILE* fp, matrix* A, matrix* B, double epsabs, double epsrel);




void timesJ (matrix* A, int p, int q, double theta);
void Jtimes (matrix* A, int p, int q, double theta);
void sweep (matrix* A, matrix* V);
int jacobi_diag (matrix* A, vector* eigenvals, matrix* eigenvecs);


#endif

