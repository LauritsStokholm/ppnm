#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>

/* ............................................................................
 * This is the method for multiple parameters in function 
 * ..........................................................................*/
struct my_f_params {double a; double b;};
typedef struct my_f_params my_f_params;


/*
int
rosenbrock(const gsl_vector* r, void * params, double * function)
{
  const double x = gsl_vector_get(r, 0);
  const double y = gsl_vector_get(r, 1);
  my_f_params * p = (my_f_params *) params;
  double a = (p->a);
  double b = (p->b);

  *function = pow((a-x), 2) + b*pow(y-pow(x,2),2);
  return GSL_SUCCESS;
}
*/
// I have changed the function for a quick test in gnuplot
double
rosenbrock(double x, double y)
{
  return pow((1-x),2) + 100*pow(y - pow(x,2),2);
}

int
gradient (const gsl_vector * x, void * params, gsl_vector * f)
{
  // This is the vectors values of domain (x0, x1)
  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);

  // Parameters for Rosenbrock function (a = 1, b = 100);
  my_f_params * p = (my_f_params *) params;
  double a = (p->a);
  double b = (p->b);

  // This is the analytically calculated gradient values (df/dx, df/dy)
  const double dfx = (-2)*((a-x0) + 2*b*(x1-pow(x0,2)));
  const double dfy = 2*b*(x1-pow(x0,2));

  gsl_vector_set(f, 0, dfx);
  gsl_vector_set(f, 1, dfy);

  return GSL_SUCCESS;
}

int
task1(void)
{
  /* ............................................................................
   * We will compare the multiple methods defined in the GSL library
   * without derivatives
   *
   * - Hybrids
   * - Hybrid
   * - dnewton
   * - broyden
   * ..........................................................................*/

  // Types of solvers
  int n = 4;
  const gsl_multiroot_fsolver_type * T0 = gsl_multiroot_fsolver_dnewton;
  const gsl_multiroot_fsolver_type * T1 = gsl_multiroot_fsolver_hybrids;
  const gsl_multiroot_fsolver_type * T2 = gsl_multiroot_fsolver_hybrid;
  const gsl_multiroot_fsolver_type * T3 = gsl_multiroot_fsolver_broyden;
  const gsl_multiroot_fsolver_type* T[] = {T0, T1, T2, T3};

  // Dimension of the system
  size_t dimensions = 2;
  my_f_params params = {1, 100};

  // File pointers
  FILE* fp0 = fopen("dnewton.txt", "w");
  FILE* fp1 = fopen("hybrids.txt", "w");
  FILE* fp2 = fopen("hybrid.txt", "w");
  FILE* fp3 = fopen("broyden.txt", "w");
  FILE* fp[] = {fp0, fp1, fp2, fp3};

  // Setting the sytstem
  gsl_multiroot_function F;
  F.f = &gradient;
  F.n = dimensions;
  F.params = &params;

  // My guess:
  gsl_vector * x = gsl_vector_alloc(dimensions);
  gsl_vector_set(x, 0, -2);
  gsl_vector_set(x, 1, 3);


  for (int i=0; i<n; i++){
  // Solver
  gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc(T[i], dimensions);

  // Solve system (0th time)
  gsl_multiroot_fsolver_set(s, &F, x);

  // prepare for print
  int status;
  size_t iter = 0;

  // headers for data file
  fprintf(fp[i], "iter\tx\ty\tf(x,y)\tgradx(x,y)\tgrady(x,y)\n");
  fprintf(fp[i], "%lu\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
      iter, 
      gsl_vector_get (s->x, 0),
      gsl_vector_get (s->x, 1),
      rosenbrock(gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1)),
      gsl_vector_get (s->f, 0),
      gsl_vector_get (s->f, 1));

  // iterate
  do
  {
    ++iter;
    status = gsl_multiroot_fsolver_iterate (s);

    fprintf(fp[i], "%lu\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
        iter,
        gsl_vector_get (s->x, 0),
        gsl_vector_get (s->x, 1),
        rosenbrock(gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1)),
        gsl_vector_get (s->f, 0),
        gsl_vector_get (s->f, 1));

    if (status)
      break;

    status = gsl_multiroot_test_residual (s->f, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < 1000);

//  fprintf(fp[i],"status = %s\n", gsl_strerror (status));
  gsl_multiroot_fsolver_free (s);
  }
  gsl_vector_free(x);

  fclose(fp0); fclose(fp1); fclose(fp2); fclose(fp3);

  return 0;
}

