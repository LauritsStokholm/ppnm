
// Defining ODES
// Example 1: y'' = -y
void f (double t, vector* y, vector* dydt, void* params)
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

void
task_a (void)
{
  FILE* fp = fopen ("task_a.txt", "w");

  // Implementations are found in ode.h
  fprintf (fp, "As an example; we consider u'' = -u\n");

  // System interval
  double t0 = 0, t1 = 2*M_PI;
  // h (step-size), abs and eps (absolute, relative precision)
  double dt = 0.1, epsabs = 1e-2, epsrel = 1e-2;
  int iter_max = 100;
  int ODE_type = 12; // rkstep12

  // ODE vector
  vector* y = vector_alloc (2);
  vector_set (y, 0, 0);
  vector_set (y, 1, 1);

  my_variable x;
  x.ode = y;
  x.params = NULL;

  fprintf (fp, "t\ty[0]\ty[1]\tcos\tsin\titer\n");
  double t = t0;
  while (t < t1)
  {
    if ( t+dt> t1 ) { dt = t1 - t; }
    printf ("t out in while: %lg\n", t);

    int iter = lineback (f, t, x, t0, t1, &dt, iter_max, ODE_type, epsabs, epsrel);
    t += dt;
    fprintf (fp, "%lg\t%lg\t%lg\t%lg\t%lg\t%i\n",
    t,
    vector_get (y, 0),
    vector_get (y, 1),
    cos(t),
    sin(t),
    iter);

  }

  vector_free (y);
  fclose (fp);
}

