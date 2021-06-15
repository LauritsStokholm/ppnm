#ifndef HAVE_ODE
#define HAVE_ODE

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>


/* ............................................................
 * Notation:
 * yt = y(t)                  Current evaluation of y at t
 * h                          The step to be taken
 * yh = y(t+h)                Next evaluation
 * dy                         Error in y
 *
 * The system considered it  dy(t)/dt = f(t, y)
 *
 * NOTE: Even though n = 2, y = {y0, y1} will be updated for each
 * iteration as described by the given ODE.
 *
 * NOTE 2: Even though f has a dependency of x, it is mostly for
 * notation; as it reminds us where y is evaluated.
 * ............................................................ */


// rkstep12
/* ............................................................
 * Description:
 * Given domain values x and h; and codomain list of values y(x)
 * calculates dydx at x, to approximate a half-step
 * y(t) = y(x+h/2) = y(x) + y'(x)*(h/2)
 * thus to calculate dydx at half-step t=x+h/2; such to estimate
 * for the full step y(x+h).
 *
 * Error in this embeded method is estimated as difference in y(t)
 * after two half steps compared to a single full step.
 * ............................................................ */
void
rkstep12 (void f(double x, vector* yx, vector* dydx),
          double x, vector* y, double h, vector* yh, vector* dy)
{
  //assert ( (yt->size == yh->size) && (yh->size == dy->size) );
  int i;
  int n = y->size;
  double val;

  vector* yt  = vector_alloc (n); // Intermediary values of y(t)
  vector* k0  = vector_alloc (n); // List for dydx for 1st step
  vector* k12 = vector_alloc (n); // List for dydx for 2nd step 

  // Estimation of y(t) a half step away from x (beginning)
  f(x, y, k0); // Calculates dydx at x (stored in k0)
  for (i=0; i<n; i++)
  {
    val = vector_get(y, i) + vector_get(k0, i)*h/2.; // Estimate y(x+h/2)
    vector_set(yt, i, val);
  }

  // Estimation of y(t) another half steps away in succession from before
  f(x+h/2., yt, k12); // dydx at x+h/2
  for (i=0; i<n; i++)
  {
    val = vector_get(y, i) + vector_get(k12, i)*h; // Estimate y(x+h)
    vector_set(yh, i, val);
  }

   // Error estimation by embeded method
  for (i=0; i<n; i++)
  {
    val = (vector_get (k0, i) - vector_get (k12, i))*h/2.;
    vector_set(dy, i, val);
  }

  vector_free (yt);
  vector_free (k0);
  vector_free (k12);
}


// For description; refer to rkstep12.
void
rkstep23 (void f(double x, vector* y, vector* dydx),
          double x, vector* y, double h, vector* yh, vector* dy)
{
  int i;
  int n = y->size;
  double val, error;

  vector* k1 = vector_alloc (n);
  vector* k2 = vector_alloc (n);
  vector* k3 = vector_alloc (n);
  vector* k4 = vector_alloc (n);
  vector* yt = vector_alloc (n);

  f(x, y, k1); // Calculate dydx for 1st step (stored in k1)
  for (i=0; i<n; i++)
  {
    // Estimate y(x+h/2) = y(x) + (h/2)*y'(x) linearly
    val = vector_get (y, i) + (1./2)*vector_get(k1, i)*h;
    vector_set (yt, i, val);
  }

  f(x + (1./2)*h, yt, k2); // Calculate dydx for 2nd step (stored in k2)
  for (i=0; i<n; i++)
  {
    // Estimate y(x+3/4h) = y(x) + (3h/4)*y'(x)
    val = vector_get (y, i) + (3./4)*vector_get (k2, i)*h;
    vector_set (yt, i, val);
  }

  f(x + (3./4)*h, yt, k3); // Set k3 to dydx at three quarters step
  for (i=0; i<n; i++)
  {
    val = vector_get (y, i)
      + ((2./9)*vector_get(k1, i)
      +  (1./3)*vector_get(k2, i)
      +  (4./9)*vector_get(k3, i))*h;

    vector_set (yh, i, val);
  }

  f(x+h, yh, k4); // Set k4 to dydx at end x (full step)
  for (i=0; i<n; i++)
  {
    val = vector_get (y, i)
          +((7./24)*vector_get (k1, i)
          + (1./4) *vector_get (k2, i)
          + (1./3) *vector_get (k3, i)
          + (1./8) *vector_get (k4, i))*h;
    vector_set (yt, i, val);

    error = vector_get(yh, i) - vector_get (yt, i);
    vector_set (dy, i, error);
    vector_set (y, i, val);

  }

  vector_free (k1);
  vector_free (k2);
  vector_free (k3);
  vector_free (k4);
  vector_free (yt);
}

int
my_ode_driver (void f(double t, vector* y, vector* dydx),
            vector* y, double a, double b, double h,
            double abs, double eps, int max,
            int toggle_switch, const char* string_name)
{
  // rkXY method (12 or 23)
  assert (toggle_switch == 12 || toggle_switch == 23);

  char charname[50];
  if (toggle_switch == 12) {strcpy(charname, "rk12_");}
  else {strcpy(charname, "rk23_");}

  strcat (charname, string_name);
  printf ("%s\n", charname);

  FILE* fp = fopen(charname, "w");

  int i, k=0;
  double error, norm_y, tol, s;

  int n = y->size;
  vector* yh = vector_alloc (n);
  vector* dy = vector_alloc (n);

  // From x[0] = a; to x[end]=b,
  double x = a;
  while (x < b)
  {
    //yval = vector_get (y, k);
    if ((x + h) > b) { h = b-x; }

    if (toggle_switch == 12) { rkstep12(f, x, y, h, yh, dy); }
    if (toggle_switch == 23) { rkstep23(f, x, y, h, yh, dy); }

    // Calculate error
    s = 0;
    for (i=0; i<n; i++)
    {
      s += vector_get(dy, i)* vector_get(dy, i);
    }
    error = sqrt(s);

    // Calculate norm_y
    s = 0;
    for (i=0; i<n; i++)
    {
      s += vector_get (yh, i) * vector_get(yh, i);
    }
    norm_y = sqrt(s);

    // Calculate tolerance
    tol = (norm_y *eps + abs)*sqrt(h/(b-a));
    if (error < tol) /* accept step and continue */
    {
      k++;

      // Max has been reached
      if (k>max-1) { return -k; }

      //vector_set (x, k, x+h);
      x = x+h;
      fprintf (fp, "%lg\t", x);
      for (i=0; i<n; i++)
      {
        vector_set (y, i, vector_get (yh, i));
        fprintf(fp, "%lg\t", vector_get(y, i));
      }
      fprintf(fp, "\n");
    }

    if (error > 0)
    {
      h *= pow(tol/error, 0.25)*0.95;
    }
    else h *= 2;
  } /* end while */

  vector_free (yh);
  vector_free (dy);

  fclose (fp);
  return k+1;
}

#endif
