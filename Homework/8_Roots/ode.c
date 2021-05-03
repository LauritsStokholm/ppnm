#ifndef HAVE_ODE
#define HAVE_ODE
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
 *
 * ............................................................ */
void
rkstep12 (void f(double x, vector* y, vector* yprime, void* params),
          double x, void* params, vector* y, double step_size,
          vector* dy, vector* y2)
{
  int n = y->size;
  double val;

  vector* yt  = vector_alloc (n); // Intermediary values of y(t)
  vector* k0  = vector_alloc (n); // List for dydx for 1st step
  vector* k12 = vector_alloc (n); // List for dydx for 2nd step 

  // Estimation of y(t) a half step away from x (beginning)
  f (x, y, k0, params); // Calculates dydx at x (stored in k0)
  for (int i=0; i<n; i++)
  {
    val = vector_get (y, i) + vector_get (k0, i)*step_size/2.; // Estimate y(x+h/2)
    vector_set(yt, i, val);
  }

  // Estimation of y(t) another half steps away in succession from before
  f(x+step_size/2., yt, k12, params); // dydx at x+h/2
  for (int i=0; i<n; i++)
  {
    val = vector_get(y, i) + vector_get(k12, i)*step_size; // Estimate y(x+h)
    vector_set (y2, i, val);
  }

   // Error estimation by embeded method
  for (int i=0; i<n; i++)
  {
    val = (vector_get (k0, i) - vector_get (k12, i))*step_size/2.;
    vector_set (dy, i, val);
  }

  vector_free (yt);
  vector_free (k0);
  vector_free (k12);
  return; // returns yh (evolution of y(x) at y(x+h) and dy (error)
}


// For description; refer to rkstep12.
void
rkstep23 (void f(double x, vector* y, vector* yprime, void* params),
          double x, void* params, vector* y, double step_size,
          vector* dy, vector* y2)
{
  int n = y->size;
  double val, error;

  vector* k1 = vector_alloc (n);
  vector* k2 = vector_alloc (n);
  vector* k3 = vector_alloc (n);
  vector* k4 = vector_alloc (n);
  vector* yt = vector_alloc (n);

  f(x, y, k1, params); // Calculate dydx for 1st step (stored in k1)
  for (int i=0; i<n; i++)
  {
    // Estimate y(x+h/2) = y(x) + (h/2)*y'(x) linearly
    val = vector_get (y, i) + (1./2)*vector_get(k1, i)*step_size;
    vector_set (yt, i, val);
  }

  f(x + (1./2)*step_size, yt, k2, params); // Calculate dydx for 2nd step (stored in k2)
  for (int i=0; i<n; i++)
  {
    // Estimate y(x+3/4h) = y(x) + (3h/4)*y'(x)
    val = vector_get (y, i) + (3./4)*vector_get (k2, i)*step_size;
    vector_set (yt, i, val);
  }

  f(x + (3./4)*step_size, yt, k3, params); // Set k3 to dydx at three quarters step
  for (int i=0; i<n; i++)
  {
    val = vector_get (y, i)
      + ((2./9)*vector_get(k1, i)
      +  (1./3)*vector_get(k2, i)
      +  (4./9)*vector_get(k3, i))*step_size;

    vector_set (y2, i, val);
  }

  f(x+step_size, y2, k4, params); // Set k4 to dydx at end x (full step)
  for (int i=0; i<n; i++)
  {
    val = vector_get (y, i)
          +((7./24)*vector_get (k1, i)
          + (1./4) *vector_get (k2, i)
          + (1./3) *vector_get (k3, i)
          + (1./8) *vector_get (k4, i))*step_size;
    vector_set (yt, i, val);

    error = vector_get(y2, i) - vector_get (yt, i);
    vector_set (dy, i, error);
    vector_set (y, i, val);

  }

  vector_free (k1);
  vector_free (k2);
  vector_free (k3);
  vector_free (k4);
  vector_free (yt);
  return; // returns yh (evolution of y(x) at y(x+h) and dy (error)
}

int
my_ode_driver
  (
    void f (double x, vector* y, vector* yprime, void* params),
    double x, void* params, vector* y,
    double a, double b, double* step_size, int iter_max, int ODE_type,
    double epsabs, double epsrel
  )
{
  // rkXY method (12 or 23)
  assert (ODE_type == 12 || ODE_type == 23);

  int n = y->size;
  int iter = 0;
  vector* y2 = vector_alloc (n); // y (x+h)
  vector* dy = vector_alloc (n); // error

  while (iter < iter_max)
  {
    // Sets y2 = y(t+h)
    if (ODE_type == 12) { rkstep12 (f, x, params, y, *step_size, dy, y2); }
    if (ODE_type == 23) { rkstep23 (f, x, params, y, *step_size, dy, y2); }

    // Calculate comparison values for current node
    double y_error = vector_norm (dy);
    double norm_y  = vector_norm (y2);
    double tol     = (norm_y*epsrel + epsabs) * sqrt( *step_size / (b-a) );

    // Max has been reached
    if ( iter > iter_max-1 )
    {
      fprintf (stderr, "fun: status=%i\n", iter);
      return -iter;
    }

    if (y_error < tol) /* accept step and continue */
    {
      iter++; // No more than iter_max
      vector_memcpy (y2, y);
      printf ("tol=%lg\n", tol);
      printf ("y=%lg\n", vector_get(y, 0));
      break;
    }

    // Adjust step_size
    if (y_error > 0)
    {
      *step_size *= pow(tol/y_error, 0.25)*0.95;
    }
    else {*step_size *= 2;}
  }

    vector_free (y2);
    vector_free (dy);

    return 0;//iter+1;
}
#endif
