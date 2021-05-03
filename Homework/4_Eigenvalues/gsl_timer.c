#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

void
matrix_set_random_symmetric (gsl_matrix* A)
{
  if (A->size1 != A->size2)
  {fprintf(stderr, "Symmetry err; dimensions are not equal.");}

  srand(time(NULL)); // Initialisation of RNG
  double val;
  for (int i=0; i<A->size1; i++)
  {
    // Set diagonal
    val = rand() / (double) RAND_MAX;// in [0, 1]
    gsl_matrix_set (A, i, i, val);

    for (int j=i+1; j<A->size1; j++)
    {
      // New random value
      val = rand() / (double) RAND_MAX;
      // Set off-diagonal
      gsl_matrix_set (A, i, j, val);
      gsl_matrix_set (A, j, i, val);
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
  gsl_matrix* A_cpy = gsl_matrix_alloc (n,  n); // Randomised matrix
  gsl_matrix* V     = gsl_matrix_alloc (n,  n); // Eigenvectors
  gsl_matrix* AV    = gsl_matrix_alloc (n,  n); // Eigenvectors
  gsl_matrix* D     = gsl_matrix_alloc (n,  n); // Randomised matrix
  gsl_vector* e     = gsl_vector_alloc (n);     // Eigenvalues

  matrix_set_random_symmetric (A);
  gsl_matrix_memcpy (A_cpy, A);
  gsl_matrix_set_identity (V);
  gsl_vector_set_zero (e);

  fprintf (fp, "Randomised matrix\n");
  matrix_fprintf (fp, A);
  matrix_fprintf (fp, A_cpy);

  // Diagonalisation of A
  gsl_eigen_symmv_workspace* W = gsl_eigen_symmv_alloc (n);
  gsl_eigen_symmv (A, e, V, W);

  fprintf (fp, "\nAfter diagonalisation\n");
  matrix_fprintf (fp, A);
  fprintf (fp, "\nWith eigenspace\n");
  matrix_fprintf (fp, V);
  fprintf (fp, "\nWith eigenvalues\n");
  vector_fprintf (fp, e);


  // Test of diagonalisation:
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, A_cpy, V, 0, AV);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1, V, AV, 0, D);
  fprintf (fp, "\nTest of diagonalisation\n");
  matrix_fprintf (fp, D);




  fclose(fp);
  gsl_matrix_free (A);
  gsl_matrix_free (V);
  gsl_matrix_free (AV);
  gsl_matrix_free (D);
  gsl_vector_free (e);
  gsl_eigen_symmv_free (W);

  return 0;
  }

