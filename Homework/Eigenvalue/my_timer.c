#include "main.h"
#include "library.h"
#include "utilities.h"

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
  int n = 5;
  FILE* fp = fopen ("task_c.txt", "w");

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


  // Do:

  /*-- Generate a random symmetric (n, n) matrix A --*/
  matrix* A     = matrix_alloc (n, n);
  matrix* A_cpy = matrix_alloc (n, n);
  matrix* AV    = matrix_alloc (n, n);
  matrix* V     = matrix_alloc (n, n);
  matrix* VT    = matrix_alloc (n, n);
  matrix* VTAV  = matrix_alloc (n, n);
  matrix* VTV   = matrix_alloc (n, n);
  matrix* VVT   = matrix_alloc (n, n);
  matrix* ID    = matrix_alloc (n, n);

  vector* eigenvals = vector_alloc (n);

  matrix_set_identity (ID);
  matrix_set_identity (V);
  vector_set_zero (eigenvals);

  matrix_set_random_symmetric (A);
  fprintf(fp, "Setting random symmetric matrix\n");
  matrix_fprintf (fp, A);
  fprintf (fp, "\n\n");

  matrix_memcpy (A, A_cpy); // Store A in A_cpy

  //sweep (A_cpy, V);
  // Changes A_cpy and V by Jacobi operations
  int sweep_counter = jacobi_diag (A_cpy, eigenvals, V);

  fprintf(fp, "after %i sweeps; eigenvalues calculated:\n", sweep_counter);
  vector_fprintf(fp, eigenvals);
  fprintf(fp, "\n\n");


  matrix_transpose_memcpy (V, VT);
  matrix_matrix_product (A, V, AV);
  matrix_matrix_product (VT, AV, VTAV);

  fprintf(fp, "Calculated value of VT*A*V\n");
  matrix_fprintf(fp, VTAV);

  matrix_matrix_product (V, VT, VVT);
  matrix_matrix_product (VT, V, VTV);

  fprintf(fp, "\n\n\n TESTING IMPLEMENTATION \n");
  fprintf (fp, "Is V*VT == IDENTITY ?\n");
  matrix_fprintf (fp, VVT);
  fprintf_result_matrix_matrix_test_equal (fp, VVT, ID, 1e-15, 1e-15);

  fprintf (fp, "Is VT*V == IDENTITY ?\n");
  matrix_fprintf (fp, VTV);
  fprintf_result_matrix_matrix_test_equal (fp, VTV, ID, 1e-15, 1e-15);


  matrix_free (A);
  matrix_free (A_cpy);
  matrix_free (AV);
  matrix_free (V);
  matrix_free (VT);
  matrix_free (VTAV);
  matrix_free (VTV);
  matrix_free (VVT);
  matrix_free (ID);
  vector_free (eigenvals);
  fclose(fp);

  return 0;
  }

