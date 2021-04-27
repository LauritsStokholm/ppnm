#include<stdio.h>
#include<gsl/gsl_integration.h>
/* ............................................................................
 * Calculate numerically the integral
 *
 * integral_{0}^{1}  ln(x)/sqrt(x) dx
 *
 * (We quickly notice that ln(x) is not defined for x=0 and that is a
 * singularity for 1/sqrt(x) aswell.)
 * ..........................................................................*/
int calls2 = 0;

// We first define our functional (the function to be integrated)
double
integrand(double x, void* params){
  calls2++;
  double f = log(x)/sqrt(x);
  return f;
}

int
my_gsl_integration(void)
{

  // We initialize a workspace for our integration
  // Maximal number of subintervals
  size_t n = 1000;
  gsl_integration_workspace * W = gsl_integration_workspace_alloc(n);

  // We initialize the result, error objects to hold the values after
  // integration
  double result, error;

  // The structure of gsl_function to be integrated
  gsl_function F;
  F.function = &integrand;
  F.params = NULL;

  gsl_integration_qags(&F, 0, 1, 0, 1e-7, n, W, &result, &error);

  printf("The results from the gsl implementation can be found here:\n");
  printf("result = %.18f\n", result);
  printf("estimated error =%.18f\n", error);
  printf("intervals = %.18lu\n", W->size);
  printf ("calls = %d", calls2);


  // When we are done with our allocated memory, we free it.
  gsl_integration_workspace_free(W);
  return 0;
}

int
my_inf_gsl_integration (double a, double b,
    double (*functional)(double x, void* params))
{
  // We initialize a workspace for our integration
  // Maximal number of subintervals
  size_t n = 1000;
  double epsabs = 1e-7, epsrel = 1e-7;
  gsl_integration_workspace * W = gsl_integration_workspace_alloc(n);

  // We initialize the result, error objects to hold the values after
  // integration
  double result, error;

  // The structure of gsl_function to be integrated
  gsl_function F;
  F.function = functional;
  F.params = NULL;

  if ((-1)*isinf(a) && isinf(b))
  {
    gsl_integration_qagi (&F, epsabs, epsrel, n, W, &result, &error);
  }
  else if (isinf(b))
  {
    gsl_integration_qagiu(&F, a, epsabs, epsrel, n, W, &result, &error);
  }
  else if ((-1)*isinf(a))
  {
    gsl_integration_qagil(&F, b, epsabs, epsrel, n, W, &result, &error);
  }
  else
  {
  gsl_integration_qags(&F, a, b, epsabs, epsrel, n, W, &result, &error);
  }

  printf("The results from the gsl implementation can be found here:\n");
  printf("result = %.18f\n", result);
  printf("estimated error =%.18f\n", error);
  printf("intervals = %.18lu\n", W->size);
  //printf ("calls = %d", calls2);


  // When we are done with our allocated memory, we free it.
  gsl_integration_workspace_free(W);
  return 0;
}
