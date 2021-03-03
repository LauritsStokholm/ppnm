/* ............................................................................
 * We do not need any guarding for multiple includes, as this is a convention
 * written in all the header files includes. 
 * ..........................................................................*/

#include<assert.h>
// For linear algebra:
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

void fprint_matrix_3x3(FILE* fp, char* s, gsl_matrix* A)
{
  fprintf(fp, "%s\n", s);
  for (int i=0; i<A->size1; i++)
  {
    for (int j=0; j<A->size2; j+=3)
    {
      double A1 = gsl_matrix_get(A, i, j);
      double A2 = gsl_matrix_get(A, i, j+1);
      double A3 = gsl_matrix_get(A, i, j+2);
      fprintf(fp, "[%lg\t%lg\t%lg]\n", A1, A2, A3);
    }
  }
}


/* ............................................................
 * This compares two vectors for equality up to an absolute/relative epsilon
 * return 0 if not equal
 * return 1 if equal
 * ............................................................ */
int
my_gsl_vector_equal(gsl_vector* a, gsl_vector* b, double eps_abs, double eps_rel)
{
  assert(a->size == b->size);
  double x; double y;

  for (int i=0; i<a->size; i++)
  {
    x = fabs ( gsl_vector_get(a, i) - gsl_vector_get(b, i) );
    y = fabs ( gsl_vector_get(a, i) + gsl_vector_get(b, i) );

    if ((x < eps_abs) || (x/y < eps_rel))
    {
      return 1;
    }
  }
  return 0;
}

void task1_linalg(void){
  /* ..........................................................................
   * Consider the following system of linear equations in the matrix form
   *
   * [ 6.13  -2.90   5.86 ][x0]     [6.23]
   * [ 8.08  -6.31  -3.89 ][x1]  =  [5.37]
   * [-4.36   1.00   0.19 ][x2]     [2.29]
   *
   * Solve it using GSL.
   *
   * In this script all methods as listed in the exercise is tried, and
   * matrices are allocated by different methods.
   * ....................................................................... */

  /* We prepare 4 matrices for : 1) HH, 2) LU, 3) QR and 4) to check answers */

  /* Dimensions*/
  size_t n = 3;
  size_t m = 3;

  /* File pointer */
  FILE* fp0 = fopen("input_linalg.txt", "r");
  FILE* fp1 = fopen("output_linalg.txt", "w");

  /* 1st method: Read file for matrix input */
  gsl_matrix* M1 = gsl_matrix_alloc(n, m);
  gsl_matrix_fscanf(fp0, M1);
  fclose(fp0);

  /* 2nd method: Matrix views */
  /* Matrix views are objects, so to obtain the matrix, we point to
   * matrixcomponent */
  double matrix_data[] = { 6.13, -2.90,  5.86,
                           8.08, -6.31, -3.89,
                          -4.36,  1.00,  0.19};

  gsl_matrix_view M2_view = gsl_matrix_view_array(matrix_data, n, m);
  gsl_matrix* M2 = &M2_view.matrix;

  /* 3rd method: Matrix copies */
  gsl_matrix* M3 = gsl_matrix_alloc(n, m);
  gsl_matrix* M4 = gsl_matrix_alloc(n, m);

  gsl_matrix_memcpy(M3, M1);
  gsl_matrix_memcpy(M4, M2);

  /* Now we check our 4 matrices by printing to output file */
// This is the build-in functions to view matrices
//  fprintf(fp1, "Here we check all matrices\n");
//  gsl_matrix_fprintf(fp1, M1, "%lg");
//  gsl_matrix_fprintf(fp1, M2, "%lg");
//  gsl_matrix_fprintf(fp1, M3, "%lg");
//  gsl_matrix_fprintf(fp1, M4, "%lg");

// This is a custommade view of matrices:
  fprintf(fp1, "Custom check\n");
  fprint_matrix_3x3(fp1, "Custom check matrix 1:", M1);
  fprint_matrix_3x3(fp1, "Custom check matrix 2:", M2);
  fprint_matrix_3x3(fp1, "Custom check matrix 3:", M3);
  fprint_matrix_3x3(fp1, "Custom check matrix 4:", M4);

  /* Now we allocate vectors */
  double b_data[] = {6.23, 5.37, 2.29};

  gsl_vector_view b_view = gsl_vector_view_array(b_data, n);
  gsl_vector* b = &b_view.vector;

  gsl_vector* x1 = gsl_vector_alloc(n);
  gsl_vector* x2 = gsl_vector_alloc(n);
  gsl_vector* x3 = gsl_vector_alloc(n);


  /* Methods to solve linear systems */
/*****************************************************************************/
  /* Method 1: Householders method */
  gsl_linalg_HH_solve(M1, b, x1);

/*****************************************************************************/
  /* Method 2: LU method (Lower-Upper triangular) */
  /* We store signum (sign of permutations) and a permutation matrix, which is
  * used for solving the system. */
  int signum;
  gsl_permutation* p = gsl_permutation_alloc(n);
  gsl_linalg_LU_decomp(M2, p, &signum);
  gsl_linalg_LU_solve(M2, p, b, x2);

/*****************************************************************************/
   /* Method 3: QR method */
  gsl_vector* tau = gsl_vector_alloc(n);
  gsl_linalg_QR_decomp(M3, tau);
  gsl_linalg_QR_solve(M3, tau, b, x3);

/*****************************************************************************/
  /* Note that the matrices are ruined, after using any of the methods */
  fprintf(fp1, "Matrices are ruined, as we can see: \n");
  fprint_matrix_3x3(fp1, "Custom check matrix 1:", M1);
  fprint_matrix_3x3(fp1, "Custom check matrix 2:", M2);
  fprint_matrix_3x3(fp1, "Custom check matrix 3:", M3);
  fprint_matrix_3x3(fp1, "Custom check matrix 4:", M4);

  // Solutions: (Solving by CBLAS)
  gsl_vector* y1 = gsl_vector_alloc(n);
  gsl_vector* y2 = gsl_vector_alloc(n);
  gsl_vector* y3 = gsl_vector_alloc(n);

   /* Calculate   */
  gsl_blas_dgemv(CblasNoTrans, 1., M4, x1, 0., y1);
  gsl_blas_dgemv(CblasNoTrans, 1., M4, x2, 0., y2);
  gsl_blas_dgemv(CblasNoTrans, 1., M4, x3, 0., y3);

   /* Testing */
  fprintf(fp1, "Testing results from each method: b is given:\n");
  gsl_vector_fprintf(fp1, b, "%lg");

  fprintf(fp1, "And the solutions for each method is:\n\n\n");

  fprintf(fp1, "Solution from HH method\n");
  gsl_vector_fprintf(fp1, x1, "%lg");
  fprintf(fp1, "and checking the solution for Ax=b\n");
  gsl_vector_fprintf(fp1, y1, "%lg");
  if( my_gsl_vector_equal(b, y1, 1e-3, 1e-3) ){fprintf(fp1, "Successful solution\n\n\n");}
  else{fprintf(fp1, "Not a solution\n\n\n");}

  fprintf(fp1, "Solution from LU method\n");
  gsl_vector_fprintf(fp1, x2, "%lg");
  fprintf(fp1, "and checking the solution for Ax=b\n");
  gsl_vector_fprintf(fp1, y2, "%lg");
  if( my_gsl_vector_equal(b, y2, 1e-3, 1e-3) ){fprintf(fp1, "Successful solution\n\n\n");}
  else{fprintf(fp1, "Not a solution\n\n\n");}

  fprintf(fp1, "Solution from QR method\n");
  gsl_vector_fprintf(fp1, x3, "%lg");
  fprintf(fp1, "and checking the solution for Ax=b\n");
  gsl_vector_fprintf(fp1, y3, "%lg");

  if( my_gsl_vector_equal(b, y3, 1e-3, 1e-3) ){fprintf(fp1, "Successful solution\n\n\n");}
  else{fprintf(fp1, "Not a solution\n\n\n");}


  gsl_permutation_free(p);
  gsl_vector_free(tau);
  gsl_vector_free(x1);
  gsl_vector_free(x2);
  gsl_vector_free(x3);
  gsl_vector_free(y1);
  gsl_vector_free(y2);
  gsl_vector_free(y3);
  gsl_matrix_free(M1);
  gsl_matrix_free(M3);
  gsl_matrix_free(M4);

// These are saved on the stack as they were defined by views, and so they
// shall not be freed.
//gsl_matrix_free(M2);
//gsl_vector_free(b);
}
