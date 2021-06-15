#include "main.h"
#include "adapt.c"
#include "gsl.c"

int calls = 0;
// This has been done
double f (double x) {calls++; return sqrt(x); }
double g (double x) {calls++; return 4*sqrt(1-x*x); }

// This is used for comparison
double h (double x) {calls++; return 1./sqrt(x); }
double i (double x) {calls++; return log(x) / sqrt(x); }

double gaussian (double x) {calls++; return exp(-x*x);}

double gsl_gaussian (double x, void* params)
{
  calls++;
  return exp(-x*x);
}


int
main (void)
{
  printf ("Task A: Implementation of a recursive adaptive integrator\n");
  printf ("Task B: Implementation of clenshaw curtis variable transf\n");
  printf ("Task C: Include integration error, and generalise to infinite limits\n");
  // Set interval for integral limits
  double a = 0, b = 1;
  // Tolerance of integrals
  double abs = 1e-4, rel = 1e-4;

  // Results of integrator
  double result, error;

  // Sqrt(x)
  printf ("integrand\tmethod\t\ta\tb\testimated result\terror\t\t\tactual result\tnumber of calls\n");

  adaptive_integrator (f, a, b, abs, rel, &result, &error);
  printf ("sqrt(x)\t\t(adapt)\t\t%lg\t%lg\t%.15lg\t%.15lg\t%lg\t%d\n" ,a, b, result, error, 2./3, calls);

  calls = 0;
  clenshaw_curtis (f, a, b, abs, rel, &result, &error);
  printf ("sqrt(x)\t\tclenshaw-curtis\t%lg\t%lg\t%.15lg\t%.15lg\t%lg\t%d\n" ,a, b, result, error, 2./3, calls);
  printf ("\n\n");

  // 4*Sqrt(1-x*x)
  calls=0;
  adaptive_integrator (g, a, b, abs, rel, &result, &error);
  printf ("4sqrt(1-xx)\t(adapt)\t\t%lg\t%lg\t%.15lg\t%.15lg\t%lg\t\t%d\n" ,a, b, result, error, M_PI, calls);

  calls=0;
  clenshaw_curtis (g, a, b, abs, rel, &result, &error);
  printf ("4sqrt(1-xx)\tclenshaw-curtis\t%lg\t%lg\t%.15lg\t%.15lg\t%lg\t\t%d\n" ,a, b, result, error, M_PI, calls);
  printf ("\n\n");

  // 1/sqrt(x)
  calls=0;
  adaptive_integrator (h, a, b, abs, rel, &result, &error);
  printf ("1/sqrt(x)\t(adapt)\t\t%lg\t%lg\t%.15lg\t%.15lg\t%d\t\t%d\n", a, b, result, error, 2, calls);

  calls=0;
  clenshaw_curtis (h, a, b, abs, rel, &result, &error);
  printf ("1/sqrt(x)\tclenshaw-curtis\t%lg\t%lg\t%.15lg\t%.15lg\t%d\t\t%d\n", a, b, result, error, 2, calls);
  printf ("\n\n");

  // log(x)/sqrt(x)
  calls=0;
  adaptive_integrator (i, a, b, abs, rel, &result, &error);
  printf ("log(x)/sqrt(x)\t(adapt)\t\t%lg\t%lg\t%.15lg\t%.15lg\t%d\t\t%d\n", a, b, result, error, -4, calls);

  calls=0;
  clenshaw_curtis (i, a, b, abs, rel, &result, &error);
  printf ("log(x)/sqrt(x)\tclenshaw-curtis\t%lg\t%lg\t%.15lg\t%.15lg\t%d\t\t%d\n", a, b, result, error, -4, calls);
  printf ("\n\n");
  printf ("NOTE: Notice how the last two examples differ in calls, since log(x) and 1./sqrt(x) is not defined at x=0\n\n");

  printf ("Calculate the integral of 4*sqrt(1-x*x) over [0, 1], using GSL\n");
  calls=0;
  my_gsl_integration();

  printf ("\n\n\n");
  printf ("Task C: Infinite Limits\n");
  printf ("integrand\tmethod\t\ta\tb\t\testimated result\terror\t\t\tactual result\t\t#calls\n");
  calls=0;
  adapt_inf (gaussian, -INFINITY, INFINITY, abs, rel, &result, &error);
  printf ("gaussian\t(adapt)\t-INFINITY\tINFINITY\t%.15lg\t%.15lg\t%.15lg\t%d\n", result, error, sqrt(M_PI), calls);
  printf ("\n\n");
  calls=0;
  adapt_inf (gaussian, 0, INFINITY, abs, rel, &result, &error);
  printf ("gaussian\t(adapt)\t0\t\tINFINITY\t%.15lg\t%.15lg\t%.15lg\t%d\n", result, error, sqrt(M_PI)/2., calls);
  printf ("\n\n");
  calls=0;
  adapt_inf (gaussian, -INFINITY, 0, abs, rel, &result, &error);
  printf ("gaussian\t(adapt)\t-INFINITY\t0\t\t%.15lg\t%.15lg\t%.15lg\t%d\n", result, error, sqrt(M_PI)/2., calls);
  printf ("\n\n");

  printf ("infinite limits\n");
  my_inf_gsl_integration(-INFINITY, INFINITY, &gsl_gaussian);
  printf ("\n\n");
  my_inf_gsl_integration(0, INFINITY, &gsl_gaussian);
  printf ("\n\n");
  my_inf_gsl_integration(-INFINITY, 0, &gsl_gaussian);
  printf ("\n\n");

  return 0;
}
