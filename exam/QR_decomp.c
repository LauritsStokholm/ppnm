


int
QR_decomp (matrix* A, matrix* Q, matrix* R)
{
  // Setting parameters
  int n = A->size1;
  int m = A->size2;
  double epsabs = 1e-10, epsrel = 1e-10;
  FILE* fp = fopen ("QR_output.txt", "w");

  // Allocate memory
  matrix* QT  = matrix_alloc (m, n);
  matrix* QR  = matrix_alloc (n, m);
  matrix* QTQ = matrix_alloc (m, m);
  matrix* ID_n = matrix_alloc (n, n);
  matrix* ID_m = matrix_alloc (m, m);

  fprintf(fp, "Set A randomly (x: 0<= x <=1)\n");
  matrix_set_zero   (Q);
  matrix_set_zero   (R);
  matrix_set_identity(ID_n);
  matrix_set_identity(ID_m);

  matrix_fprintf (fp, A);

  // QR- decomposition
  fprintf (fp, "\n\n\nTesting QR-decomposition by GS-method\n");
  GS_decomp (A, Q, R);

  // Results
  fprintf (fp, "Q\n");
  matrix_fprintf(fp, Q);

  fprintf(fp, "\nR\n");
  matrix_fprintf(fp, R);

  fprintf(fp, "\nQT\n");
  matrix_transpose_memcpy (Q, QT);
  matrix_fprintf (fp, QT);

  // Testing conditions
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
  fprintf_result_matrix_matrix_test_equal (fp, QTQ, ID_m, epsabs, epsrel);

  //matrix_free (R);    matrix_free (Q);
  matrix_free (QR);   matrix_free (QT); matrix_free (QTQ);
  matrix_free (ID_n); matrix_free (ID_m);

  return 0;
}
