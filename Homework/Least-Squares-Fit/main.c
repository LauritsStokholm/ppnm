#include "main.h"
#include "library.c"
#include "utilities.c"

// Define the system of functions
int k = 2; // Number of function switches.

// Set of functions
double function_switch (int i, double x)
{
  switch(i){
    case 0: return 1; break;
    case 1: return x; break;
    default: return NAN;
  }
}


// Ordinary Least-Squares Fit by QR decomp (using GS)
void
lsfit (vector* x, vector* y, vector* dy, funlst* F)
{
  int n = x->size; int m = F->size;
  matrix* A = matrix_alloc (n, m); matrix_set_zero (A);
  vector* b = vector_alloc (n);
  double (*f)(int, double) = F->function;

  for (int i=0; i<n; i++)
  {
    double xi  = vector_get (x, i);
    double yi  = vector_get (y, i);
    double dyi = vector_get (dy, i);
    vector_set (b, i, yi/dyi);

    for (int j=0; j<m; j++)
    {
      double val = f(j, xi) / dyi;
      matrix_set (A, i, j, val);
    }
  }

  // Solve system by QR-decomp using GS method
  matrix* Q = matrix_alloc (n, m); matrix_set_zero (Q);
  matrix* R = matrix_alloc (m, m); matrix_set_zero (R);
  matrix* QR_inverse = matrix_alloc (m, n);

  GS_decomp (A, Q, R);
  GS_inverse (Q, R, QR_inverse);

  vector* c = vector_alloc (m);
  matrix_vector_product (QR_inverse, b, c);

  // Covariance
  matrix* R_inverse = matrix_alloc (m, m);
  matrix* R_inverse_transposed = matrix_alloc (m, m);
  matrix* Covariance = matrix_alloc (m, m);

  matrix_inverse_upper_triangular (R, R_inverse);
  matrix_transpose_memcpy(R_inverse, R_inverse_transposed);
  matrix_matrix_product(R_inverse, R_inverse_transposed, Covariance);

  F->coefficients = c;
  F->covariance = Covariance;

  vector_free (b);
  matrix_free (A);
  matrix_free (Q);
  matrix_free (R);
  matrix_free (QR_inverse);
  matrix_free (R_inverse);
  matrix_free (R_inverse_transposed);
}


int
main (int argc, char* argv[])
{

  funlst F;
  F.size = k;
  F.function = &function_switch;


  // DATA (quick and dirty - otherwise; read from Rutherford.data)
  double x_list[] = {1,   2,   3,  4,  6,  9,    10,   13,   15};
  double y_list[] = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};

  vector* x  = vector_alloc (9);
  vector* y  = vector_alloc (9);
  vector* dy = vector_alloc (9);

  for (int i=0; i<9; i++)
  {
    vector_set (x, i, x_list[i]);
    vector_set (y, i, log(y_list[i]));
    vector_set (dy, i, 1./20.);
  }

  // Calculate coefficients
  lsfit (x, y, dy, &F);

  FILE* fp = fopen("myplot.data", "w");
  fprintf (fp, "x\tfmin\tfval\tfmax\n");


  double ci, ci_err;
  for (double dx=x_list[0]; dx<x_list[8]; dx+=1e-3)
  {
    double fmin = 0, fval = 0, fmax = 0;
    for (int i=0; i<F.size; i++)
    {
      ci = vector_get (F.coefficients, i);
      ci_err = sqrt(matrix_get (F.covariance, i, i));
      fmin += (ci - ci_err) * F.function(i, dx);
      fval +=  ci           * F.function(i, dx);
      fmax += (ci + ci_err) * F.function(i, dx);
    }
  fprintf(fp, "%lg\t%lg\t%lg\t%lg\n", dx, exp(fmin), exp(fval), exp(fmax));
  }

  // We have the decay constant, lambda, in the expression exp(-lambda*t),
  // which the half-life time is related by t_1/2 = log(2) / lambda
  // Since, we fitted exp(lambda*t), our parameter differs by sign.

  double c0 = vector_get (F.coefficients, 0);
  double c1 = vector_get (F.coefficients, 1);
  double c0_err = sqrt(matrix_get (F.covariance, 0, 0) );
  double c1_err = sqrt(matrix_get (F.covariance, 1, 1) );

  printf("\n\n\n Fitted coefficients \n");
  printf("c0 \t \t = \t \t %lg +/- %lg\n", c0, c0_err);
  printf("c1 \t \t = \t \t %lg +/- %lg\n", c1, c1_err);
  printf("\n\n\n");
  printf ("The half-life-time of ThX is (by fit):    \t\t t_1/2 = %lg +/- %lg \t[days]\n", log(2)/(-c1), c1_err*log(2)/(c1*c1));
  printf ("The half-life-time of 224Ra is (by table):\t\t t_1/2 = %lg +/- %lg \t[days]\n", 3.6319, 0.0023);
  printf("\n\n\n");


  fclose(fp);
  vector_free (x);
  vector_free (y);
  vector_free (dy);

  return 0;
}
