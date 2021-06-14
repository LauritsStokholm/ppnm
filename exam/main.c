// C Standard libraries for API
#include <stdio.h>     // Standard i/o functions,  i.e. printf
#include <stdlib.h>    // Standard lib. functions, i.e. RANDOM
#include <math.h>      // Mathematical library
#include <assert.h>    // assertions used for error-corrections
#include <time.h>      // OS time; used for random seed
#include <getopt.h>    // Get variable options from command line

// Self-implemented libraries
#include "main.h"      // header-file for functions, structures
#include "library.h"   // Library for linear algebra implementation
#include "utilities.h" // The functions used in SVD
#include "SVD.c"       // The one-sided Jacobi (SVD) algorithm
#include "QR_decomp.c" // The QR-decomposition algorithm (for n>m)
#include "GS_decomp.c" // The Gram-Schmitt algorithm for QR_decomp

/* ............................................................
 * The main.c file: calling all seperate functions; here only
 * the one sided jacobi algorithm for Singular Value Decomposition
 * (SVD) is called.
 * ............................................................ */

int
main (int argc, char* argv[])
{
  /* Default setting of parameters */
  // By default; look at small, tall (thin) matrix
  int n = 3;   // Number of rows
  int m = n-1; // number of columns
  // These can of course be changed by ./main -n X -m Y

  // File pointer for output-stream
  FILE* fp = fopen ("results.txt", "w");

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

  /* -- Begin calculations -- */
  matrix* A = matrix_alloc (n, m);
  matrix* U     = matrix_alloc (n, m); // Orthogonal tall
  matrix* D     = matrix_alloc (m, m); // Orthg, square, diag. singular values
  matrix* V     = matrix_alloc (m, m); // Orthogonal square
  matrix* VT    = matrix_alloc (m, m); // V transposed
  matrix* UD    = matrix_alloc (n, m); // Product matrix U*D
  matrix* UDVT  = matrix_alloc (n, m); // Product matrix U*D*VT


  // The generated matrix; subject to algorithm
  fprintf(fp, "Setting stochastic matrix\n");
  matrix_set_random (A);
  matrix_fprintf (fp, A);
  fprintf (fp, "\n\n");


  // Algorithm assumes tall/thin or square matrix
  if (n < m) { fprintf(stderr, "Error: Algorithm assumes tall/thin or square matrix\n"); }

  // If square; do SVD algorithm
  if (n == m) {one_sided_jacobi_SVD (A, U, D, V); }

  // If tall/thin; do QR first then SVD algorithm
  else
  {
    matrix* Uprime = matrix_alloc (m, m);
    matrix* Q = matrix_alloc (n, m);
    matrix* R = matrix_alloc (m, m);
    QR_decomp (A, Q, R); // Q (orthogonal), R (upper-triangular)
    one_sided_jacobi_SVD (R, Uprime, D, V);
    matrix_matrix_product (Q, Uprime, U);
  }

  fprintf (fp, "Converged: Orthg. equality and other tests can be found in respective .txt files.\n");
  fprintf (fp, "Here only equality is checked:\n");
  matrix_transpose_memcpy (V, VT);
  matrix_matrix_product (U, D, UD);
  matrix_matrix_product (UD, VT, UDVT);

  fprintf (fp, "Is A = U*D*VT?\n");
  fprintf_result_matrix_matrix_test_equal (fp, A, UDVT, 1e-10, 1e-10);


  matrix_free (A);
  matrix_free (U); matrix_free (D); matrix_free (V);
  matrix_free (VT); matrix_free (UD); matrix_free (UDVT);
  return 0;
}
