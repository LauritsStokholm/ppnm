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
