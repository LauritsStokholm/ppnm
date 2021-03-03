
/* ............................................................
 * Using one of the GSL integration routines implement the bessel
 * function using its integral representation. Make a plot.
 * ............................................................ */

#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <gsl/gsl_integration.h>

/* ............................................................
 * Multiple parameter function needs a typedef struct
 * ............................................................ */
typedef struct {int n; double x;} my_params;

/* ............................................................
 * Defining integrand: Error function integrand
 * ............................................................ */
double
my_jn (double t, void* params)
{
  int n = ((my_params*) params)->n;
  double x = ((my_params*) params)->x;

  return cos(n*t - x*sin(t)) / M_PI;
}

int
main (int argc, char** argv)
{
  // Unpacking parameters
  // Default
  double xmin = 0, xmax = M_PI; // Integral limits
  double epsabs = 1e-9, epsrel = 1e-9; // Integral result tolerance
  double N = 1e4; // Interval partition number

  FILE* fp;
  char* filenames[] = {"j0.txt", "j1.txt", "j2.txt"};

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

  my_params parameters;
  gsl_function F;
  F.function = &my_jn;
  F.params = &parameters;

  /* Workspaces */
  gsl_integration_workspace* W0         = gsl_integration_workspace_alloc(limit_N);
  gsl_integration_cquad_workspace* W1   = gsl_integration_cquad_workspace_alloc(limit_N);
  gsl_integration_romberg_workspace* W2 = gsl_integration_romberg_alloc(limit_N);
  gsl_integration_glfixed_table* W3     = gsl_integration_glfixed_table_alloc(limit_N);

  for (int n=0; n<3; n++)
  {
    fp = fopen (filenames[n], "w");
    fprintf(fp, "x\tres0\tres1\tres2\tres3\tfunc\n");
    parameters.n = n;
    for (double x=0; x<20; x+=1e-2)
    {
      parameters.x = x;
  gsl_integration_qags (&F, xmin, xmax, epsabs, epsrel, limit_N, W0, &res0, &err);
  gsl_integration_cquad(&F, xmin, xmax, epsabs, epsrel, W1, &res1, &err, &neval);
  gsl_integration_romberg(&F, xmin, xmax, epsabs, epsrel, &res2, &neval, W2);
  res3 = gsl_integration_glfixed(&F, xmin, xmax, W3);

  fprintf(fp, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", x, res0, res1, res2, res3, jn(n, x));
    }
  }
  gsl_integration_workspace_free(W0);
  gsl_integration_cquad_workspace_free(W1);
  gsl_integration_romberg_free(W2);
  gsl_integration_glfixed_table_free(W3);
  fclose(fp);



  return 0;
}
