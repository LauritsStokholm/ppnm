/* ............................................................................
   SHOOTING METHOD:
   Given a boundary value problem of an ODE, the shooting method reduces the
   problem to one finding roots.

   y''(t) = f(t, y(t), y'(t)),    y(t0) = y0,     y(t1) = y1;
   Let y(t; a) denote the solution of the initial value problem
   y''(t) = f(t, y(t), y'(t)),    y(t0) = y0      y(t1) = a
   Then define the function F(a) as the difference between y(t1; a) and the
   specified boundary value y1
   F(a) = y(t1; a) - y1.
   If F has a root, we are done.


 * First we will find the numerical solution f(r) to the differential equation.
 * This is the s-wave radial Schroedinger equation for the Hydrogen atom in
 * units Bohr radius and Hartree.
 *
 * -(1/2)f'' - (1/r)f = e*f
 *
 *  where the bound s-state wave-function satisfies this equation and the two
 *  boundary conditions
 *  f(r->0)   = r - r**2
 *  f(r->inf) = 0
 *
 *  we will replace inf with a large number rmax = 10 (bohr radii), since one
 *  cannot integrate numerically to inf, and as rmax is much larger than the
 *  typical size of the H-atom, but still managable for the numerical
 *  integrator. The second boundary condition then reads:
 *  f(rmax) = 0
 *
 *  Both boundary conditions are only satisfied for certain discrete negeative
 *  values of the energy. We will only use the first boundary condition, to
 *  define our solution to the ODE, and then define an auxilary function
 *  M(eps) = F_eps (rmax)
 *  This function of energy, is then analysed, and when finding the root
 *  M(eps) = 0, we have solved the energies.
 *
 *  a) Find the lowest root epsilon_0 of the equation M(epsilon) = 0 for, say
 *  rmax=8. Plot the resulting function and compare with the exact result,
 *  which is epsilon_0 = -1/2 and f(r)=r*exp(-r)
 *
 * ..........................................................................*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include"task2.h"


int
caller(int toggle, double rmax)
{
  // Basics
  int status, iter = 0, iter_max = 1000;
  size_t dim = 1;

  // Setting system
  gsl_multiroot_function F;
  if(toggle == 0){F.f = &M_eps;}
  if(toggle == 1){F.f = &M_eps_nonzero;}
  F.n = dim;
  F.params = &rmax;

  // My guess of epsilon (negative energy)
  gsl_vector * x = gsl_vector_alloc(dim);
  gsl_vector_set (x, 0, -1);

  // Choose method
  const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrids;

  // Solver
  gsl_multiroot_fsolver * s= gsl_multiroot_fsolver_alloc(T, dim);

  // Solve system (0th time)
  gsl_multiroot_fsolver_set(s, &F, x);

  do
  {
    iter++;
    int flag = gsl_multiroot_fsolver_iterate(s);
    if (flag != 0)
      break;
    status = gsl_multiroot_test_residual(s->f, 1e-3);
  }
  while (status == GSL_CONTINUE && iter < iter_max);

  double e = gsl_vector_get (s->x, 0);
  FILE * fp = fopen("task2.txt", "w");
  fprintf(fp, "rmax\te\n");
  fprintf(fp, "%g\t%g\n", rmax, e);

  fclose(fp);
  FILE * fp2 = fopen("data2.txt", "w");
  fprintf(fp2, "r\tf(e,r)\tExact\n");

  for (double r=0; r<=rmax; r+=0.1)
  {
    fprintf(fp2, "%lg\t%lg\t%lg\n", r, f_e(e,r), r*exp(-r));
  }

  // Free up
  fclose(fp2);
  gsl_vector_free(x);
  gsl_multiroot_fsolver_free(s);

  return iter;
}

// b) test the algorithm for convergence dependency of rmax
int
convergence(int toggle)
{
  assert(toggle<2);
  assert(toggle>-1);
  FILE* fp;
  if (toggle == 0){fp = fopen("convergence.txt", "w");}
  if (toggle == 1){fp = fopen("precision.txt", "w");}
  fprintf(fp, "rmax\titer\n");
  for (double rmax=2; rmax < 8; rmax+=0.1)
  {
    int iter = caller(toggle, rmax);
    fprintf(fp, "%lg\t%i\n", rmax, iter);
  }
  return GSL_SUCCESS;
}

int
task2(void)
{
  // b)
  convergence(0);
  // We also try a more precise boundary condition for bound states
  convergence(1);



  return GSL_SUCCESS;
}
