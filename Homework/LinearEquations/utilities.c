// Gram-Schmidt decomposition
/* ............................................................
 * Given a matrix A, return the QR decomposition using the Gram-Schmidt
 * orthogonalisation method (GS-method).
 * Upon return, A stores Q and R stores the computed matrix R s.t
 * A = QR, QTQ = QQT = 1 and R is upper-triangular.
 * ............................................................ */
void GS_decomp (matrix* A, matrix* Q, matrix* R)
{
  matrix_memcpy (A, Q);

  //double aij, aik, dot;
  double qij, aik, dot;

  // Iterate over column-vectors of A;
  for (int j=0; j<A->size2; j++)
  {
    // Calculate norm of column-vector
    double norm = matrix_column_norm (Q, j);

    // Set diagonal of R
    matrix_set (R, j, j, norm);

    // Normalise column vector of A (qi = ai/norm(ai))
    for (int i=0; i<A->size1; i++)
    {
      qij = matrix_get (Q, i, j) / norm;
      matrix_set (Q, i, j, qij);
    }

    // Orthogonalise the column-vector from rest of set
    // and calculate dot product for fixed qj, and variable ak (row of R)
      for (int k=j+1; k<A->size2+1; k++)// +1
      {
          // Calculate dot-product
          dot = 0;
          // elementwise multiplication
          for (int l=0; l<A->size1; l++)
          {
            dot += matrix_get (Q, l, j) * matrix_get (Q, l, k);
          }
          matrix_set (R, j, k, dot);

        // Elementwise
        for (int i=0; i<A->size1; i++)
        {
          aik = matrix_get (Q, i, k);
          qij = matrix_get (Q, i, j);

          // Elementwise GS
          matrix_set(Q, i, k, aik - dot*qij);
        }
      }
  }
}


/* ............................................................
 * Gram-Schmidt solver
 * Given matrices Q, R from GS_decomp, and vectors b and x,
 * solve the system
 * Q*R*x = b
 * by applying QT to the vector b, saving the result in vector x,
 * QT(Q*R*x) = QT*b, as QTQ = 1
 * R*x = QT*b
 * and then performing backward/in-place substitution on x
 * (as R is upper triangular).
 *
 * returns x
 * ............................................................ */
void GS_solve (matrix* Q, matrix* R, vector* b, vector* x)
{
  matrix* QT = matrix_alloc (Q->size2, Q->size1);
  matrix_transpose_memcpy (Q, QT);

  matrix_vector_product (QT, b, x);
  backsub (R, x);
  matrix_free (QT);
  // Stores result in x
}

/* ............................................................
 * Given matrices Q and R from GS_decomp calculates inverse of
 * A = Q*R
 * inv(A) = inv(Q*R) = inv(R)*inv(Q) = inv(R)*QT
 * and stores this in QRI
 * ............................................................ */
void GS_inverse (matrix* Q, matrix* R, matrix* inv_QR)
{
  matrix* inv_Q = matrix_alloc (Q->size2, Q->size1);
  matrix* inv_R = matrix_alloc (R->size2, R->size1);

  // Calculate inverses
  matrix_transpose_memcpy (Q, inv_Q);
  matrix_inverse_upper_triangular (R, inv_R);

  // inv_QR = inv_R * inv_Q
  matrix_matrix_product (inv_R, inv_Q, inv_QR);

  free (inv_Q);
  free (inv_R);
}
