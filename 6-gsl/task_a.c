
/* ............................................................
 * Using the GSL integration routine QAGS calculate numerically
 * the integral: int_{0}^{1} dx (log(x) / sqrt(x))
 * ............................................................ */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>


/* ............................................................
 * Defining integrand
 * ............................................................ */
double
f (double x, void* params)
{
  return log(x) / sqrt(x);
}

/* ............................................................
 * Numerically integration method: QAGS
 * xmin, xmax: integration limits [xmin, xmax)
 * result, error: value of integral and its numerical error
 * epsabs, epsrel: absolute/relative tolerance of integration
 * limit_N: Number of interval partitions
 * ............................................................ */
int
main (int argc, char* argv[])
{
  double xmin = 0, xmax = 1;
  double result, error;
  double epsabs = 1e-6, epsrel = 1e-6;
  size_t limit_N = 1e3;

  FILE* fp = fopen ("task_a.txt", "w");
  fprintf(fp, "Task a: Integration using QAGS:\n");

 /* Allocate workspace and integrate */
  gsl_integration_workspace* W = gsl_integration_workspace_alloc(limit_N);

 /* Function structure for integral */
  double parameters = 0;
  gsl_function F;
  F.function = &f;
  F.params = &parameters;

  gsl_integration_qags (&F, xmin, xmax, epsabs, epsrel, limit_N, W, &result, &error);

 /* Results */
  fprintf(fp, "Result\tError\n");
  fprintf(fp, "%lg\t%lg\n\n\n", result, error);

 /* Free memory after allocation */
  gsl_integration_workspace_free(W);
  fclose(fp);

  return 0;
}
