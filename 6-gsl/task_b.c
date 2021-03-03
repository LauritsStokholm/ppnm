
/* ............................................................
 * Using one of the GSL integration routines implement the error
 * function using its integral representation. Make a plot.
 * ............................................................ */

#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <gsl/gsl_integration.h>


/* ............................................................
 * Defining integrand: Error function integrand
 * ............................................................ */
double
f (double x, void* params)
{
  return 2./sqrt(M_PI) * exp(-x*x);
}

int
main (int argc, char** argv)
{
  // Unpacking parameters
  // Default
  double xmin = 0, xmax = 1;
  double epsabs = 1e-2, epsrel = 1e-2;
  double N = 1e3;

  FILE* fp = fopen ("task_b.txt", "w");

  // Possible setting parameters
  while (1)
  {
    int opt = getopt (argc, argv, "a:b:e:t:N");
    if (opt == -1) break;
    switch (opt)
    {
      case 'a': xmin = atof (optarg); break;
      case 'b': xmax = atof (optarg); break;
      case 'e': epsabs = atof (optarg);
                epsrel = atof (optarg); break;
      case 'N': N = atof (optarg); break;
      default:
      fprintf (stderr,
        "Usage: %s [-a start] [-b end] [e: absolute/relative error limits] [-N: partition number]\n",
        argv[0]);
    exit (EXIT_FAILURE);
    }
  }

  /* Maximal number of double precision intervals (partition of interval)*/
  /* and function evaluation number (for some integration methods only) */
  size_t limit_N = (size_t) N;
  size_t neval   = (size_t) N;

  /* Initiialise results and error (we wont save errors) */
  double res0, res1, res2, res3, err;

  double parameters = 0;
  gsl_function F;
  F.function = &f;
  F.params = &parameters;

  fprintf(fp, "x\tres0\tres1\tres2\tres3\tfunc\n");
  /* Workspaces */
  gsl_integration_workspace* W0         = gsl_integration_workspace_alloc(limit_N);
  gsl_integration_cquad_workspace* W1   = gsl_integration_cquad_workspace_alloc(limit_N);
  gsl_integration_romberg_workspace* W2 = gsl_integration_romberg_alloc(limit_N);
  gsl_integration_glfixed_table* W3     = gsl_integration_glfixed_table_alloc(limit_N);

  double dx = (xmax - xmin) / N;
  for (double x=xmin; x<xmax; x+=dx){
  gsl_integration_qags (&F, xmin, x, epsabs, epsrel, limit_N, W0, &res0, &err);
  gsl_integration_cquad(&F, xmin, x, epsabs, epsrel, W1, &res1, &err, &neval);
  gsl_integration_romberg(&F, xmin, x, epsabs, epsrel, &res2, &neval, W2);
  res3 = gsl_integration_glfixed(&F, xmin, x, W3);

  fprintf(fp, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", x, res0, res1, res2, res3, erf(x));
  }
  gsl_integration_workspace_free(W0);
  gsl_integration_cquad_workspace_free(W1);
  gsl_integration_romberg_free(W2);
  gsl_integration_glfixed_table_free(W3);

  fprintf(stdout, "\n\n\n");
  fclose(fp);


  return 0;
}
