#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <assert.h>


#include <gsl/gsl_linalg.h> // SVG algorithm
#include <gsl/gsl_vector.h> // Vector structure
#include <gsl/gsl_matrix.h> // Matrix structure
//#include <gsl/gsl_eigen.h>
//#include <gsl/gsl_blas.h>
/* ............................................................
 * This flie implements the GSL defined SVG method. This is 
 * only done for a comparison for my own implementation. A convergence
 * test can be found in Makefile (the plot).
 * ............................................................ */

// Stochastic matrix
void
matrix_set_random (gsl_matrix* A)
{
  // Test for square matrix
  if (A->size1 != A->size2)
  {fprintf(stderr, "Symmetry err; dimensions are not equal.");}

  // Set elements stochastically
  srand(time(NULL)); // Initialisation of RNG
  double val;
  for (int i=0; i<A->size1; i++)
  {
    for (int j=0; j<A->size2; j++)
    {
    val = rand() / (double) RAND_MAX;// in [0, 1]
    gsl_matrix_set (A, i, j, val);
    }
  }
}

void matrix_fprintf (FILE* fp, gsl_matrix* A)
{
  for (int i=0; i<A->size1; i++)
  {
    int j = 0;
    while (j < A->size2)
    {
      fprintf(fp, "%.4lg\t", gsl_matrix_get(A, i, j));
      j++;
      if (j == A->size2){fprintf(fp, "\n");}
    }
  }
}

void vector_fprintf (FILE* fp, gsl_vector* v)
{
  for (int i=0; i<v->size; i++)
  {
    fprintf(fp, "%.4lg\t", gsl_vector_get(v, i));
    if (i+1 == v->size){fprintf(fp, "\n");}
  }
}

int
main (int argc, char* argv[])
  {
   /* ............................................................
   * Parameters:
   * (n, n) dimensions of matrix in (row, column) type.
   *
   * epsabs and epsrel is for testing equality of objects within
   * absolut and relative epsilon
   *
   * fp is a filepointer, used as alternative to pipelines
   * ./main > file
   * ............................................................ */

  // Setting parameters
  int n = 3;
  FILE* fp = fopen ("gsl.txt", "w");

  /* Check for specified parameters */
  while (1)
  {
    int opt = getopt (argc, argv, "n:");
    if (opt == -1) break;
    switch (opt)
    {
      case 'n': n = atoi (optarg); break;
      default:
      fprintf (stderr, "Usage: %s [-n dimensions of matrix]\n", argv[0]);
      exit (EXIT_FAILURE);
    }
  }


  gsl_matrix* A     = gsl_matrix_alloc (n,  n); // Randomised matrix

  gsl_matrix* U     = gsl_matrix_alloc (n,  n); // Decomposistion
  gsl_matrix* D     = gsl_matrix_alloc (n,  n); // Decomposistion
  gsl_matrix* V     = gsl_matrix_alloc (n,  n); // Decomposistion
  gsl_matrix* DVT   = gsl_matrix_alloc (n,  n); // Product D*VT
  gsl_matrix* UDVT  = gsl_matrix_alloc (n,  n); // Product U*D*VT
  gsl_vector* S     = gsl_vector_alloc (n);     // Singular values

  matrix_set_random (A);
  gsl_matrix_memcpy (U, A);
  gsl_matrix_set_identity (V);
  gsl_matrix_set_zero (D);
  gsl_vector_set_zero (S);

  fprintf (fp, "Randomised matrix\n");
  matrix_fprintf (fp, A);

  // Diagonalisation of A
  gsl_linalg_SV_decomp_jacobi (U, V, S); // input: U=A; output U is corrected

  // Set diagonal of singular values
  for (int i=0; i<n; i++)
  {
    gsl_matrix_set (D, i, i, gsl_vector_get (S, i));
  }

  // Testing:
  fprintf (fp, "\nAfter SVG algorithm\n");
  gsl_blas_dgemm (CblasNoTrans, CblasTrans,   1, D, V,   0, DVT);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, U, DVT, 0, UDVT);
  fprintf (fp, "\nTest of diagonalisation\n");
  matrix_fprintf (fp, UDVT);


  fclose(fp);
  gsl_matrix_free (A);
  gsl_matrix_free (U);
  gsl_matrix_free (D);
  gsl_matrix_free (V);
  gsl_matrix_free (DVT);
  gsl_matrix_free (UDVT);
  gsl_vector_free (S);

  return 0;
  }

