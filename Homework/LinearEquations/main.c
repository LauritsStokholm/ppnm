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
    matrix* R   = matrix_alloc (m, m);
    matrix* QT  = matrix_alloc (m, n);
    matrix* QR  = matrix_alloc (n, m);
    matrix* QTQ = matrix_alloc (m, m);
    matrix* ID_nn = matrix_alloc (n, n);
    matrix* ID_mm = matrix_alloc (m, m);

    fprintf(fp, "Set A randomly (x: 0<= x <=1)\n");
    matrix_set_random (A);
    matrix_set_zero   (Q);
    matrix_set_zero   (R);
    matrix_set_identity(ID_nn);
    matrix_set_identity(ID_mm);

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

    fprintf(fp, "QT*Q\n");
    matrix_fprintf(fp, QTQ);

    fprintf(fp, "QR\n");
    matrix_fprintf(fp, QR);

    fprintf(fp, "\n\n\n");
    fprintf(fp, "TEST: A==QR with epsabs=%lg, epsrel=%lg?\n", epsabs, epsrel);
    fprintf_result_matrix_matrix_test_equal (fp, A, QR, epsabs, epsrel);

    fprintf(fp, "TEST: QT*Q==IDENTITY with epsabs=%lg, epsrel=%lg?\n", epsabs, epsrel);
    fprintf_result_matrix_matrix_test_equal (fp, QTQ, ID_nn, epsabs, epsrel);


    // Task 2: Test GS-solve
    fprintf(fp, "\n\n\nTask 2: Test GS-solve\n");
    vector* x = vector_alloc (m); // Solution A*x = (QR)*x = b
    vector* b = vector_alloc (n); // RHS (n, m)(m) = (n)
    vector* QRx = vector_alloc (n);

    fprintf (fp, "Set b randomly\n");
    vector_set_random (b);
    vector_fprintf(fp, b);

    fprintf (fp, "Solve the system A*x = Q*R*x = b\n");
    GS_solve (Q, R, b, x);

    fprintf (fp, "Is the result QRx, as calculated by GS_solve, equal to b?\n");
    matrix_vector_product (QR, x, QRx);
    fprintf (fp, "QRx\n");
    vector_fprintf(fp, QRx);

    fprintf (fp, "\nTEST: QRx==b with epsabs=%lg, epsrel=%lg?\n", epsabs, epsrel);
    fprintf_result_vector_vector_test_equal (fp, QRx, b, epsabs, epsrel);


    // Test GS-inverse
    matrix* inv_QR = matrix_alloc (m, n);
    matrix* id_nn = matrix_alloc (n, n);
    matrix* id_mm = matrix_alloc (m, m);

    fprintf (fp, "\n\n\nTask3: Test GS-inverse\n");
    GS_inverse(Q, R, inv_QR);

    matrix_matrix_product(QR, inv_QR, id_nn);
    matrix_matrix_product(inv_QR, QR, id_mm);

    fprintf (fp, "is inv_QR*QR = identity?\n");
    fprintf (fp, "\nTEST: inv_QR*QR==IDENTITY with epsabs=%lg, epsrel=%lg?\n", epsabs, epsrel);
    matrix_fprintf(fp, id_nn);
    fprintf_result_matrix_matrix_test_equal(fp, id_nn, ID_nn, epsabs, epsrel);

    fprintf (fp, "\nTEST: QR*inv_QR==IDENTITY with epsabs=%lg, epsrel=%lg?\n", epsabs, epsrel);
    matrix_fprintf(fp, id_mm);
    fprintf_result_matrix_matrix_test_equal(fp, id_mm, ID_mm, epsabs, epsrel);


    // Not necessary?
    //matrix* QQT = matrix_alloc (n, n);

    //matrix_matrix_product (Q, QT, QQT);
    //fprintf(fp, "Q*QT\n");
    //matrix_fprintf(fp, QQT);

    //fprintf(fp, "TEST: Q*QT==IDENTITY with epsabs=%lg, epsrel=%lg?\n", epsabs, epsrel);
    //fprintf_result_matrix_matrix_test_equal (fp, QQT, IDENTITY1, epsabs, epsrel);

    //matrix_free (QQT);


    matrix_free (A);
    matrix_free (R);
    matrix_free (Q);
    matrix_free (QR);
    matrix_free (QT);
    matrix_free (QTQ);
    matrix_free (inv_QR);
    matrix_free (ID_nn);
    matrix_free (ID_mm);
    matrix_free (id_nn);
    matrix_free (id_mm);

    vector_free (x);
    vector_free (b);
    vector_free (QRx);
    return 0;
  }

