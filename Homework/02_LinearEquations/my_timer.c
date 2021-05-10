#include "main.h"
#include "library.c"
#include "utilities.c"

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
  double epsabs = 1e-10;
  double epsrel = 1e-10;
  char* name = "output.txt";

  /* Check for specified parameters */
  while (1)
  {
    int opt = getopt (argc, argv, "n:m:a:r:f:");
    if (opt == -1) break;
    switch (opt)
    {
      case 'n': n = atoi (optarg); break;
      case 'm': m = atoi (optarg); break;
      case 'a': epsabs = atof (optarg); break;
      case 'r': epsrel = atof (optarg); break;
      case 'f': name = optarg; break;
      default:
      fprintf (stderr, "Usage: %s [-n -m dimensions of matrix] [-a -r absolute/relative tolerance for testing] [-f filepointer\n", argv[0]);
      exit (EXIT_FAILURE);
    }
  }
  FILE* fp = fopen (name, "w");

    // Allocate memory
    matrix* A   = matrix_alloc (n, m);
    matrix* Q   = matrix_alloc (n, m);
    matrix* R   = matrix_alloc (n, m);
    matrix* QT  = matrix_alloc (m, n);
    matrix* QR  = matrix_alloc (n, m);
    matrix* QTQ = matrix_alloc (n, m);
    matrix* QQT = matrix_alloc (n, m);
    matrix* IDENTITY = matrix_alloc (n, n);

    fprintf(fp, "Set A randomly (x: 0<= x <=1)\n");
    matrix_set_random (A);
    matrix_set_zero   (Q);
    matrix_set_zero   (R);
    matrix_set_identity(IDENTITY);

    matrix_fprintf (fp, A);

    // Task 1: QR- decomposition
    fprintf (fp, "\n\n\nTesting QR-decomposition by GS-method\n");
    GS_decomp (A, Q, R);

    fprintf (fp, "Q\n");
    matrix_fprintf(fp, Q);

    fprintf(fp, "\nR\n");
    matrix_fprintf(fp, R);

    fprintf(fp, "\nQT\n");
    matrix_transpose_memcpy (Q, QT);
    matrix_fprintf (fp, QT);

    fprintf (fp, "\n\n\nTesting the properties of the decomposition\n");
    matrix_matrix_product (Q, R, QR);
    matrix_matrix_product (QT, Q, QTQ);
    matrix_matrix_product (Q, QT, QQT);

    fprintf(fp, "Q*QT\n");
    matrix_fprintf(fp, QQT);

    fprintf(fp, "QT*Q\n");
    matrix_fprintf(fp, QTQ);

    fprintf(fp, "QR\n");
    matrix_fprintf(fp, QR);

    fprintf(fp, "\n\n\n");
    fprintf(fp, "TEST: A==QR with epsabs=%lg, epsrel=%lg?\n", epsabs, epsrel);
    fprintf_result_matrix_matrix_test_equal (fp, A, QR, epsabs, epsrel);

    fprintf(fp, "TEST: Q*QT==IDENTITY with epsabs=%lg, epsrel=%lg?\n", epsabs, epsrel);
    fprintf_result_matrix_matrix_test_equal (fp, QQT, IDENTITY, epsabs, epsrel);

    fprintf(fp, "TEST: QT*Q==IDENTITY with epsabs=%lg, epsrel=%lg?\n", epsabs, epsrel);
    fprintf_result_matrix_matrix_test_equal (fp, QTQ, IDENTITY, epsabs, epsrel);
    matrix_free (A);
    matrix_free (R);
    matrix_free (Q);
    matrix_free (QR);
    matrix_free (QT);
    matrix_free (QQT);
    matrix_free (QTQ);
    matrix_free (IDENTITY);
    return 0;
  }

