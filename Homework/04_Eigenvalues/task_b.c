#include "main.h"
#include "library.h"
#include "utilities.h"

void
task_b (void)
{
  FILE* fp = fopen ("energies.txt", "w");
  FILE* plot = fopen ("plot.txt", "w");

  // Build the Hamiltonian matrix
  int n = 20;
  double s = 1.0 / (n+1);
  matrix* H = matrix_alloc (n, n);
  matrix_set_zero (H);

  for (int i=0; i<n-1; i++)
  {
    matrix_set (H, i,   i,  -2);
    matrix_set (H, i,   i+1, 1);
    matrix_set (H, i+1, i,   1);
  }

  matrix_set   (H, n-1, n-1, -2);
  matrix_fprintf (fp, H);
  matrix_scale (H, -1./(s*s));


  // Diagonalise Hamiltonian using Jacobi routine
  matrix* V = matrix_alloc (n, n);
  vector* e = vector_alloc (n);
  matrix_set_identity (V);

  printf(" Jacobi_diag!\n");
  jacobi_diag (H, e, V);
  printf(" Jacobi_diag done!\n");



  fprintf(fp, "nth Energy:\tCalculated\t Exact\n");
  // Check energies
  for (int k=0; k<n/3; k++)
  {
    double exact = M_PI*M_PI*(k+1)*(k+1);
    double calculated = matrix_get (H, k, k);
    fprintf(fp, "%i\t%lg\t%lg\n",k, calculated, exact);
  }

  // Plotting several lowest eigenfunctions (to compare with analytic result)
  // Boundary condition
    fprintf (plot, "0\t0\t0\t0\t0\t0\t0\n");
    for (int i=0; i<n; i++)
    {
      fprintf(plot, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
          (i+1.0)/(n+1),
          matrix_get (V, i, 0), matrix_get (V, i, 1),
          matrix_get (V, i, 2), matrix_get (V, i, 3),
          matrix_get (V, i, 4), matrix_get (V, i, 5));
    }
    // Boundary condition
    fprintf (plot, "1\t0\t0\t0\t0\t0\t0\n");

  fclose (fp);
  fclose (plot);
  matrix_free (H);
  matrix_free (V);
}

