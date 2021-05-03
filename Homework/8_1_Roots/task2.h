#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>


// This is the s-wave radial Schrodinger equation for the H-atom
int
equation(double r, const double y[], double yprime[], void* params)
{
/* ............................................................................
* This defines the ODE (rank 2), which is written to 2 ODEs (rank 1)
* -(1/2)f'' - (1/r)f = ef
* y = {y0, y1} = (f, f')
* y'= {y0',y1'}= (f', f'') = (y1, -2*[e + (1/r))*y0)
* ...........................................................................*/

  // We only have one parameter (energy, E)
  double e = *(double*) params;
  yprime[0] = y[1];
  yprime[1] = -2*(e + (1/r))*y[0];
  return GSL_SUCCESS;
}

// This is the function f_eps(r).used in the ODE. We will vary epsilon in this,
// to find epsilon_0 s.t. f_e(rmax) = 0
double
f_e(double e, double r)
{
  // Program assumes positive r, if not, abort.
  assert(r >= 0);
  // Limit for r -> 0
  const double rmin = 1e-3;
  if (r < rmin) return r - pow(r, 2);

  // Parameter for the system; driven from t to t1 (will be changed)
  double t = rmin;

  // We have changed the ODE(2) to two ODE(1), so we use y[]
  double y[] = {t - pow(t,2), 1 - 2*t};

  // Setting up the system structure
  gsl_odeiv2_system swave;
  swave.function = equation;
  swave.jacobian = NULL;
  swave.dimension = 2;
  swave.params = (void *) &e;

  // The precision
  double epsabs = 1e-7, epsrel = 1e-7;

  // Starting step
  double hstart = 1e-3;

  // Choose a method
  const gsl_odeiv2_step_type * mystep = gsl_odeiv2_step_rk8pd;

  // Driver
  gsl_odeiv2_driver * mydriver =
    gsl_odeiv2_driver_alloc_y_new(&swave, mystep, hstart, epsabs, epsrel);

  // We test if the driver was succesful
  int status = gsl_odeiv2_driver_apply(mydriver, &t, r, y);
  if (status != GSL_SUCCESS) fprintf(stderr, "fun: status=%i", status);

  return y[0];
}

// This is the function M(eps) = f_e(rmax) used in rootfinding s.t.
// f_e(rmax) has the correct boundary condition (by shooting method)
// s toggle will toggle between the exercises a/b and c
int
M_eps (const gsl_vector * x, void * params, gsl_vector * f)
{
  // Function M(eps) = f_e(rmax) to find root eps_0 s.t. f_eps(rmax)=0
  double epsilon = gsl_vector_get (x, 0);
  // Energy is negative (bound state)
  assert(epsilon < 0);
  double r_max = *(double *) params;
  double M_epsilon = f_e(epsilon, r_max);
  //double M_epsilon = f_e(epsilon, r_max) - r*exp(-sqrt(-*2*epsilon)*r);

  gsl_vector_set (f, 0, M_epsilon);

  return GSL_SUCCESS;
}

// This is the second shooting value used in our method.
int
M_eps_nonzero (const gsl_vector * x, void * params, gsl_vector * f)
{
  // Function M(eps) = f_e(rmax) to find root eps_0 s.t. f_eps(rmax)=0
  double epsilon = gsl_vector_get (x, 0);
  // Energy is negative (bound state)
  assert(epsilon < 0);
  double r_max = *(double *) params;
  double M_epsilon = f_e(epsilon, r_max) - r_max*exp(-sqrt(-2*epsilon)*r_max);

  gsl_vector_set (f, 0, M_epsilon);

  return GSL_SUCCESS;
}
