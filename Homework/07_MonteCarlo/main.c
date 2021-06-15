#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "plainmc.c"
#include "halton.c"
#include "ssmc.c"
//#include "low-discrepancy.c"

double
f (int dim, double x[])
{
  double product = 1;
  for (int i=0; i<dim; i++)
  {
    product *= cos(x[i]);
  }
  return 1./(1-product)/(M_PI*M_PI*M_PI);
}

double
g (int dim, double x[])
{
  double product = 1;
  for (int i=0; i<dim; i++)
  {
    product *= sin(x[i]);
  }
  return product;///(M_PI*M_PI*M_PI);
}

int
main (int argc, char* argv[])
{
  double dim = 3; // Dimension of integration space
  double a[] = {0, 0, 0}; // Boundary (limit points)
  double b[] = {M_PI, M_PI, M_PI}; // Boundary (limit points)
  int N = 1e6; // Sampling points
  double result, error; // Containers for monte-carlo

  // Task A
  printf ("Calculate the tripple integral from 0 to pi of sin(x)sin(y)sin(z)\n");
  printf ("Actual result is %d\n", 8);

  plainmc (dim, g, a, b, N, &result, &error);
  printf("Using plain monte-carlo, the result is =%.15lg\n with error =%.15lg\n", result, error);
  printf("Using %d sample points\n\n", N);

  // Task B
  quasimc (dim, g, a, b, N, &result, &error);
  printf("Using a quasi-random monte-carlo, the result is =%.15lg\n with error =%.15lg\n", result, error);
  printf("Using %d sample points\n\n", N);

  // Task C
  double abstol = 1e-4, reltol = 1e-4;
  strata (dim, g, a, b, abstol, reltol, N, 0, 0, &result, &error);
  printf("Using a recursive stratified sampling monte-carlo, the result is =%.15lg\n with error =%.15lg\n", result, error);
  printf("Using %lg as absolute and relative tolerance, and %d sample points\n\n", abstol, N);


  printf("\n\n\nA more difficult beast\n\n");
  printf ("Calculate the tripple integral from 0 to pi of (1/pi)^3 1./[1-cos(x)cos(y)cos(z)]\n");
  printf ("Actual result is %.15lg\n", 1.3932039296856768591842462603255);
  plainmc (dim, f, a, b, N, &result, &error);
  printf("Using plain monte-carlo, the result is =%.15lg\n with error =%.15lg\n", result, error);
  printf("Using %d sample points\n\n", N);

  quasimc (dim, f, a, b, N, &result, &error);
  printf("Using a quasi-random monte-carlo, the result is =%.15lg\n with error =%.15lg\n", result, error);
  printf("Using %d sample points\n\n", N);


  abstol = 1e-2, reltol = 1e-2, N=1e3;
  strata (dim, f, a, b, abstol, reltol, N, 0, 0, &result, &error);
  printf("Using a recursive stratified sampling monte-carlo, the result is =%.15lg\n with error =%.15lg\n", result, error);
  printf("Using %lg as absolute and relative tolerance, and %d sample points\n\n", abstol, N);


  return 0;
}
