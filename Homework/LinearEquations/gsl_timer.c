#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#include<getopt.h>
#include <time.h>

void matrix_set_random (gsl_matrix* A)
{
  srand(time(NULL)); // Initialisation of RNG
  for (int i=0; i<A->size1; i++)
  {
    for (int j=0; j<A->size2; j++)
    {
      double val = rand() / (double) RAND_MAX;// in [0, 1]
      gsl_matrix_set (A, i, j, val);
    }
  }
}

int
main (int argc, char* argv[])
  {
   /* ............................................................
   * Parameters:
   * (n, m) dimensions of matrix in (row, column) type.
   *
   * epsabs and epsrel is for testing equality of objects within
   * absolut and relative epsilon
   *
   * fp is a filepointer, used as alternative to pipelines
   * ./main > file
   * ............................................................ */

  // Setting parameters
  int n = 4, m = 4;

  /* Check for specified parameters */
  while (1)
  {
    int opt = getopt (argc, argv, "n:m:");
    if (opt == -1) break;
    switch (opt)
    {
      case 'n': n = atoi (optarg); break;
      case 'm': m = atoi (optarg); break;
      default:
      fprintf (stderr, "Usage: %s [-n -m dimensions of matrix]\n", argv[0]);
      exit (EXIT_FAILURE);
    }
  }

    // Allocate memory
    gsl_matrix* A = gsl_matrix_alloc (n, m);
    gsl_vector* v = gsl_vector_alloc (m);

    matrix_set_random (A);
    gsl_linalg_QR_decomp (A, v);
    return 0;
  }

