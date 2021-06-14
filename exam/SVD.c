// C Standard libraries for API
#include <stdio.h>     // Standard i/o functions,  i.e. printf
#include <stdlib.h>    // Standard lib. functions, i.e. RANDOM
#include <math.h>      // Mathematical library
#include <assert.h>    // assertions used for error-corrections

// Self-implemented libraries
#include "main.h"      // header-file for functions, structures
#include "library.h"   // Library for linear algebra implementation
#include "utilities.h" // The functions used in SVD

void
one_sided_jacobi_SVD (matrix* A, matrix* U, matrix* D, matrix* V)
{
  FILE* fp = fopen ("one-sided-jacobi.txt", "w");
  // Setting general parameters
  int n = A->size1, m = A->size2;

  /*-- Allocation of matrices used --*/
  matrix* A_cpy = matrix_alloc (n, m); // Copy for algorithm (destroyed)

  // Decomposition A = U*D*VT
  //matrix* U     = matrix_alloc (n, m); // Orthogonal tall
  //matrix* D     = matrix_alloc (m, m); // Orthogonal, square diagonal (singular values)
  //matrix* V     = matrix_alloc (m, m); // Orthogonal square
  matrix* UT    = matrix_alloc (m, n); // U transposed
  matrix* VT    = matrix_alloc (m, m); // V transposed

  // Product matrices
  matrix* DVT   = matrix_alloc (m, m); // D*V
  matrix* UDVT  = matrix_alloc (n, m); // U*D*VT

  // Testing for Identities (orthogonal matrices)
  matrix* UTU   = matrix_alloc (m, m); // UT*U
  matrix* UUT   = matrix_alloc (n, n); // U*UT
  matrix* VTV   = matrix_alloc (m, m); // VT*V
  matrix* VVT   = matrix_alloc (m, m); // V*VT
  matrix* ID_n  = matrix_alloc (n, n); // Identity of size n
  matrix* ID_m  = matrix_alloc (m, m); // Identity of size m
  matrix_set_identity (ID_n); matrix_set_identity (ID_m);

  // Storage for diagonal of D
  vector* D_diag = vector_alloc (m);

  // Set values
  matrix_set_identity (V); // For iterative step 0
  matrix_set_identity (U); // non-empty
  matrix_set_identity (D); // non-empty (we set diagonal later)
  vector_set_zero     (D_diag); // non-empty values

  matrix_memcpy (A, A_cpy); // Store A in A_cpy for algorithm

  // Changes A_cpy and V by Jacobi operations
  int sweep_counter = SVD (A_cpy, V);

  for (int i=0; i<m; i++) // i'th column
  {
    double norm_ai = matrix_column_norm (A_cpy, i);
    matrix_set (D, i, i, norm_ai);
    vector_set (D_diag, i, norm_ai);

    for (int j=0; j<n; j++) // j'th row
    {
      double aji = matrix_get (A_cpy, j, i);
      matrix_set (U, j, i, aji / norm_ai);
    }
  }


  fprintf (fp, "converged after %i sweeps:\n", sweep_counter);
  vector_fprintf (fp, D_diag);
  fprintf (fp, "\n\n");

  // Test for SVD
  matrix_transpose_memcpy (V, VT);
  matrix_transpose_memcpy (U, UT);
  matrix_matrix_product (D, VT, DVT);
  matrix_matrix_product (U, DVT, UDVT);

  fprintf (fp, "Calculated value of U*A*VT\n");
  matrix_fprintf (fp, UDVT);
  fprintf (fp, "Is A = U*D*VT ?\n");
  fprintf_result_matrix_matrix_test_equal (fp, UDVT, A, 1e-15, 1e-15);


  // Testing for orthogonalities
  fprintf (fp, "U -- matrix\n");
  matrix_fprintf (fp, U);
  fprintf (fp, "V -- matrix\n");
  matrix_fprintf (fp, V);

  // Calculate products
  matrix_matrix_product (U,  UT, UUT);
  matrix_matrix_product (V,  VT, VVT);
  matrix_matrix_product (UT, U,  UTU);
  matrix_matrix_product (VT, V,  VTV);

  fprintf (fp, "\n\n\n TESTING IMPLEMENTATION FOR ORTHOGONALITY\n");
  fprintf (fp, "Is V*VT == IDENTITY ?\n");
  matrix_fprintf (fp, VVT);
  fprintf_result_matrix_matrix_test_equal (fp, VVT, ID_m, 1e-15, 1e-15);

  fprintf (fp, "Is VT*V == IDENTITY ?\n");
  matrix_fprintf (fp, VTV);
  fprintf_result_matrix_matrix_test_equal (fp, VTV, ID_m, 1e-15, 1e-15);


  fprintf (fp, "Is U*UT == IDENTITY ?\n");
  matrix_fprintf (fp, UUT);
  fprintf_result_matrix_matrix_test_equal (fp, UUT, ID_n, 1e-15, 1e-15);

  fprintf (fp, "Is UT*U == IDENTITY ?\n");
  matrix_fprintf (fp, UTU);
  fprintf_result_matrix_matrix_test_equal (fp, UTU, ID_m, 1e-15, 1e-15);


  // Free memory and close file pointers
  matrix_free (A_cpy);
 // matrix_free (U);   matrix_free (V);     matrix_free (D);
  matrix_free (UT);  matrix_free (VT);
  matrix_free (UUT); matrix_free (UTU);
  matrix_free (VVT); matrix_free (VTV);
  matrix_free (DVT); matrix_free (UDVT);
  vector_free (D_diag);
  matrix_free (ID_n); matrix_free (ID_m); 
  fclose(fp);
}


