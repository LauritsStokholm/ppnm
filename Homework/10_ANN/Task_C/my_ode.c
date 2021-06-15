#include "main.h"
#include "library.c"
#include "ode.h"

// Defining ODES
// Example 1: y'' = -y
void f(double t, vector* y, vector* dydt)
{
  /* ............................................................
   * We consider the ODE: u'' = -u
   *
   * Let y = {y0, y1} = {u, u'}
   * the yprime = {y0', y1'} = {u', u''} = {y1, -y0}
   * ............................................................ */
  vector_set (dydt, 0, vector_get (y, 1));
  vector_set (dydt, 1, -vector_get (y, 0));
}

int
main (void)
{
  // System interval
  double a = 0, b = 2*M_PI;
  // h (step-size), abs and eps (absolute, relative precision)
  double h = 0.1, abs = 1e-2, eps = 1e-2;
  int n = 2, max = 100;
  int rkXY_method;

  // Start
  vector* y = vector_alloc (n);
  vector_set (y, 0, 0);
  vector_set (y, 1, 1);

  rkXY_method = 12;
  const char* string_name = "ann_ode.txt";
  my_ode_driver (f, y, a, b, h, abs, eps, max, rkXY_method, string_name);

  vector_set (y, 0, 0);
  vector_set (y, 1, 1);
  rkXY_method = 23;
  my_ode_driver (f, y, a, b, h, abs, eps, max, rkXY_method, string_name);

  vector_free (y);

  return 0;
}

