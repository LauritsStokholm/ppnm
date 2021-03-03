/* ............................................................................
 *  Guarding for multiple includes is not necessary.
 * ..........................................................................*/
#include<gsl/gsl_eigen.h>

void task2_eigenv(void){
  /* ..........................................................................
   * (OPTIONAL) Compute the eigenvectors of 4-th order Hilbert matrix
   * H_(ij) = 1/(i+j+1)
   * ........................................................................*/

  // We define the dimensions of the system
  size_t n = 4; size_t m = 4;

  // Allocate memory for matrix
  gsl_matrix* H = gsl_matrix_alloc(n, m);

  // Set each element
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      gsl_matrix_set(H, i, j, 1./(i+j+1));
    }
  }

  // File pointer
  FILE* fp = fopen("output_eigenv.txt", "w");

  // Check H for elements
  // gsl_matrix_fprintf(fp, H, "%lg");
  fprintf(fp, "The hilbert matrix looks like:\n");
  for(int i=0; i<n; i++){
  for(int j=0; j<m; j++){
  if(j==0){fprintf(fp, "[ ");}
  fprintf(fp, "%.2lg \t", gsl_matrix_get(H, i, j));
  if(j+1==m){fprintf(fp, "]\n");}
  }}


  /* Prepare eigenvectors and eigenvalues */
  gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(4);
  gsl_vector* eigenvalue = gsl_vector_alloc(4);
  gsl_matrix* eigenvector = gsl_matrix_alloc(4, 4);

  gsl_eigen_symmv(H, eigenvalue, eigenvector, workspace);

  //fprintf(fp, "The eigenvalues are stored in the vector:\n");
  //gsl_vector_fprintf(fp, eigenvalue, "%lg");

  //fprintf(fp, "The eigenvectors are stored in the matrix:\n");
  //gsl_matrix_fprintf(fp, eigenvector, "%lg");

  fprintf(fp, "Now after the process, the matrix is ruined:\n");
  for(int i=0; i<n; i++){
  for(int j=0; j<m; j++){
  if(j==0){fprintf(fp, "[ ");}
  fprintf(fp, "%.2lg \t", gsl_matrix_get(H, i, j));
  if(j+1==m){fprintf(fp, "]\n");}
  }}
  
  fprintf(fp, "The eigenvalues are:\n");
  for(int i=0; i<n; i++)
  {
    fprintf(fp, "x(%i) = %.2lg\n", i, gsl_vector_get(eigenvalue, i));
  }

  fprintf(fp, "The eigenvectors are:\n");
  for(int i=0; i<n; i++){
  for(int j=0; j<m; j++){
  fprintf(fp, "v%i(%i) = %.2lg \t", j, i, gsl_matrix_get(eigenvector, i, j));
  if(j+1 == n){fprintf(fp, "\n");}
  }}

  gsl_matrix_free(eigenvector);
  gsl_matrix_free(H);
  gsl_vector_free(eigenvalue);
  gsl_eigen_symmv_free(workspace);
  fclose(fp);

}
