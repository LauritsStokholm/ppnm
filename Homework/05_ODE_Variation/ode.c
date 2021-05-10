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
 * ............................................................ */
void
rkstep12 (void f (double t, vector* y, vector* yprime, void* params),
          double t, my_variable x, double step_size, vector* dy, vector* y2)
{
  //assert ( (yt->size == yh->size) && (yh->size == dy->size) );
  double val;

  // Unpack my_variable
  vector* y       = x.ode;
  void* my_params = x.params;

  int n = y->size;
  vector* yt  = vector_alloc (n); // Intermediary values of y(t)
  vector* k0  = vector_alloc (n); // List for dydx for 1st step
  vector* k12 = vector_alloc (n); // List for dydx for 2nd step 

  // Estimation of y(t) a half step away from x (beginning)
  f (t, y, k0, my_params); // Calculates dydx at x (stored in k0)
  for (int i=0; i<n; i++)
  {
    val = vector_get (y, i) + vector_get (k0, i)*step_size/2.; // Estimate y(x+h/2)
    vector_set (yt, i, val);
  }

  // Estimation of y(t) another half steps away in succession from before
  f (t+step_size/2., yt, k12, NULL); // dydx at x+h/2
  for (int i=0; i<n; i++)
  {
    val = vector_get (y, i) + vector_get (k12, i)*step_size; // Estimate y(x+h)
    vector_set (y2, i, val);
  }

   // Error estimation by embeded method
  for (int i=0; i<n; i++)
  {
    val = (vector_get (k0, i) - vector_get (k12, i))*step_size/2.;
    vector_set(dy, i, val);
  }

  vector_free (yt);
  vector_free (k0);
  vector_free (k12);
}


// For description; refer to rkstep12.
//void
//rkstep23 (void f(double t, vector* y, vector* yprime, void* params),
//          double t, my_variable x, double dt, vector* y2, vector* dy)
//{
//    // Unpack my_variable
//  vector* y       = x.ode;
//  void* my_params = x.params;
//
//  int n = y->size;
//  double val, err;
//
//  vector* k1 = vector_alloc (n);
//  vector* k2 = vector_alloc (n);
//  vector* k3 = vector_alloc (n);
//  vector* k4 = vector_alloc (n);
//  vector* yt = vector_alloc (n);
//
//  f(t, y, k1, NULL); // Calculate dydx for 1st step (stored in k1)
//  for (int i=0; i<n; i++)
//  {
//    // Estimate y(x+h/2) = y(x) + (h/2)*y'(x) linearly
//    val = vector_get (y, i) + (1./2)*vector_get(k1, i)*dt;
//    vector_set (yt, i, val);
//  }
//
//  f(t + (1./2)*dt, yt, k2, NULL); // Calculate dydx for 2nd step (stored in k2)
//  for (int i=0; i<n; i++)
//  {
//    // Estimate y(x+3/4h) = y(x) + (3h/4)*y'(x)
//    val = vector_get (y, i) + (3./4)*vector_get (k2, i)*dt;
//    vector_set (yt, i, val);
//  }
//
//  f(t + (3./4)*dt, yt, k3, NULL); // Set k3 to dydx at three quarters step
//  for (int i=0; i<n; i++)
//  {
//    val = vector_get (y, i)
//      + ((2./9)*vector_get(k1, i)
//      +  (1./3)*vector_get(k2, i)
//      +  (4./9)*vector_get(k3, i))*dt;
//
//    vector_set (y2, i, val);
//  }
//
//  f(t+dt, y2, k4, NULL); // Set k4 to dydx at end x (full step)
//  for (int i=0; i<n; i++)
//  {
//    val = vector_get (y, i)
//          +((7./24)*vector_get (k1, i)
//          + (1./4) *vector_get (k2, i)
//          + (1./3) *vector_get (k3, i)
//          + (1./8) *vector_get (k4, i))*dt;
//
//    vector_set (yt, i, val);
//    err = vector_get(y2, i) - vector_get (yt, i);
//
//    // This was changed from y to yh, so not to do lots of steps
//    vector_set (dy, i, err);
//    vector_set (y2, i, val);// HERE
//  }
//
//  vector_free (k1);
//  vector_free (k2);
//  vector_free (k3);
//  vector_free (k4);
//  vector_free (yt);
//  // set vector y and dy with values y(x+h) and error
//}

int
my_ode_driver
  (
    void f (double t, vector* y, vector* yprime, void* params),
    double t, my_variable x, double step_size, vector* dy, vector* y2, int toggle_switch
  )
{
  // rkXY method (12 or 23)
  assert (toggle_switch == 12 || toggle_switch == 23);

  // Unpack my_variable
  vector* y       = x.ode;
  void* my_params = x.params;

  if (toggle_switch == 12) { rkstep12(f, t, x, step_size, dy, y2); }
  //if (toggle_switch == 23) { rkstep23(f, t, x, step_size, dy, y2); }

  return 0;
}

int
lineback
  (
    void f (double t, vector* y, vector* yprime, void* params),
    double t, my_variable x,
    double tmin, double tmax, double* step_size, int iter_max, int ODE_type,
    double epsabs, double epsrel
  )
{
  // Unpack my_variable
  vector* y = x.ode;
  //void* my_params = x.params;

  int n = y->size;

  // Basic definitions of ODE status and the iteration number
  int status = 0, iter = 0;

  // double values used for calculation of tolerance
  double y_norm, y_error, tol;

  // Vector containers for y(x+h) and errors of y
  vector* y2 = vector_alloc (y->size);
  vector* dy = vector_alloc (y->size);

  // This is to optimise step_size from t to t+h
  while (iter < iter_max )
  {
    // Sets y2 = y(t+h)
    my_ode_driver (f, t, x, *step_size, dy, y2, ODE_type);

    y_error = vector_norm (dy); // Calculate error
    y_norm  = vector_norm (y2); // Calculate norm_y
    tol = (y_norm*epsrel + epsabs)*sqrt(*step_size/(tmax-tmin)); // Calculate tolerance

    // If error is within tolerance, set y = yh
    if (y_error < tol) /* accept step and continue */
    {
      printf ("iter:\t %d\n", iter);
      iter++;

      if (iter > iter_max-1)
      {
      status = -1; return -iter; fprintf (stderr, "max iter has been reached\n");
      }
      vector_memcpy (y2, y);
      break;
    }

    if (y_error > 0)
    {
      *step_size *= pow(tol/y_error, 0.25)*0.95;
    }
    else{ *step_size *= 2;}

  } /* end while */

  vector_free (y2);
  vector_free (dy);
  return 0;//iter+1;
}

#endif
