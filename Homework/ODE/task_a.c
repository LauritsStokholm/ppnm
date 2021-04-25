
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

// Example 2: SIR model (COVID epidemic in Denmark)
void g(double t, vector* y, vector* dydt)
{
  /* ............................................................
   * dSdt = -I*S / (N*Tc)
   * dIdt = I*S / (N*Tc) - I/Tr
   * dRdt = I/Tr;
   *
   * Let y = {S, I, R}
   * dydt = {dS/dt, dI/dt, dR/dt} = {-y0y1/NTc, 
   * ............................................................ */
  double N = 6e6; // Population size from worldometers.info (2021)
  double Tc = 3;  // Average time between contacts
  double Tr = 10; // Average recovery time

  double y0 = vector_get (y, 0);
  double y1 = vector_get (y, 1);
  //double y2 = vector_get (y, 2);

  //y0*y1/(N*Tc) - y1/Tr );
  vector_set (dydt, 0, -y0*y1/(N*Tc));
  vector_set (dydt, 1, -vector_get(dydt, 0) - y1/Tr );
  vector_set (dydt, 2, -vector_get(dydt, 1) - vector_get(dydt, 0) );
}
  // Tr/Tc is the typical number of new infected per infectious individual


void
task_a ()
{
  FILE* fp = fopen ("task_a.txt", "w");

  // Implementations are found in ode.h
  fprintf (fp, "As an example; we consider u'' = -u\n");

  // System interval
  double a = 0, b = 2*M_PI;
  // h (step-size), abs and eps (absolute, relative precision)
  double h = 0.1, abs = 1e-2, eps = 1e-2;
  int n = 2, max = 100;
  int rkXY_method;

  vector* y = vector_alloc (n);
  vector_set (y, 0, 0);
  vector_set (y, 1, 1);

  // This is simulated by step size h, and the while-loop (from a to b) in
  // driver
  //double delta_x = 0.05;
  //for (double x=0; x<10; x+= delta_x);

  rkXY_method = 12;
  const char* string_name = "harmonic.txt";
  my_ode_driver (f, y, a, b, h, abs, eps, max, rkXY_method, string_name);

  vector_set (y, 0, 0);
  vector_set (y, 1, 1);
  rkXY_method = 23;
  my_ode_driver (f, y, a, b, h, abs, eps, max, rkXY_method, string_name);

  // SIR model
  n = 3;
  b = 100;
  vector* y2 = vector_alloc (n);
  vector_set (y2, 0, 5e6);  // S: suceptible
  vector_set (y2, 1, 1e6);  // I: infectious
  vector_set (y2, 2, 0);    // R: removed (immune or dead)

  rkXY_method = 12;
  const char* string_name2 = "SIR.txt";
  my_ode_driver (g, y2, a, b, h, abs, eps, max, rkXY_method, string_name2);

  vector_set (y2, 0, 5e6); // S
  vector_set (y2, 1, 1e6); // I
  vector_set (y2, 2, 0); // R
  rkXY_method = 23;
  my_ode_driver (g, y2, a, b, h, abs, eps, max, rkXY_method, string_name2);


  vector_free (y);
  vector_free (y2);
  fclose (fp);
}

