 /* find the solution F_e (r) parametrised by energy e,
 * This will in general not satisfy BDC at r=0 and r->infty
 * 
 * Fix rmax = 10 (approx infty), search the roots of varying energy
 * e, of F(e, rmax)
 * ............................................................ */


// Defining ODES
// This is the s-wave radial Schrodinger equation for the H-atom
void swave (double r, vector* y, vector* yprime, void* params)
{
  /* ............................................................
   * We consider the ODE: -1/2u'' - 1/r u = E u
   *
   * Let y = {y0, y1} = {u, u'}
   * then yprime = {y0', y1'} = {u', u''} = {y1, (-2)*(E+1/r)y0}
   * ............................................................ */
  double e = *(double*) params;
  vector_set (yprime, 0, vector_get (y, 1));
  vector_set (yprime, 1, -2*(e + (1./r))*vector_get (y, 0));
}

// This is the energy-parametrised function F_e(r) which solves the ODE
double
F_e (double e, double r)
{
  assert (r >= 0);
  // Limit for r->0
  double rmin = 1e-3;
  double rmax = 10;
  if (r < rmin){ return r - r*r;}

  double t = rmin; // Parameter for the system; to be driven from t to t1
  double step_size = 1e-3;// Starting step
  double epsabs = 1e-3, epsrel=1e-3;  // Precision / tolerance
  int iter_max = 10000;

  void* params = (void*) &e;

  // Solve ODE
  vector* y = vector_alloc (2);
  vector_set (y, 0, t - t*t);
  vector_set (y, 1, 1 - 2*t);

  int ODE_method = 12; //rkstep12 = 12 or rkstep23 = 23

  // Solve ODE and set the answer in y
  my_ode_driver (swave, r, params, y, rmin, rmax, &step_size, iter_max,
      ODE_method, epsabs, epsrel);
  double val = vector_get (y, 0);
  vector_free (y);
  return val;
}

// This is the auxiliary function used in rootfinding at rmax to find the
// correct boundary condition (by shooting method)

// x contains the parameter, epsilon and the rmax
void
M_e (vector* x, void* params, vector* fx)
{
  double epsilon = vector_get (x, 0);
  //assert ( epsilon < 0 ); // Energy is negative (bound state)
  double rmax = *(double *) params;
  double M_e = F_e (epsilon, rmax);
  vector_set (fx, 0, M_e);
  return;
}

void
task_b (void)
{
  FILE* fp = fopen ("test.txt", "w");
  fprintf (fp, "Task b\n");

  // Initialisation of coordinate and function vectors
  vector* epsilon  = vector_alloc (1);
  vector* fx = vector_alloc (1); // M_e = F_e (eps, rmax);

  // Function parameters and initial guess of root coordinates
  double rmax = 10; // rmax
  double epsilon_guess = -1;
  vector_set (epsilon, 0, epsilon_guess);

  // Acceptance tolerance
  double acc = 1e-4;

  // Change epsilon, as to set M_e (rmax) = 0
  newton (M_e, epsilon, (void*) &rmax, acc);

  fprintf (fp, "Checking the given root\n");
  M_e (epsilon, (void*) &rmax, fx);
  vector_fprintf (fp, fx);

  double e = vector_get (epsilon, 0);
  vector_fprintf (fp, epsilon);
  fprintf (fp,"epsilon = %lg\n", e);

  FILE * fp2 = fopen("data2.txt", "w");
  fprintf(fp2, "r\tf(e,r)\tExact\n");

  for (double ri=0; ri<=rmax; ri+=0.001){
    fprintf(fp2, "%lg\t%lg\t%lg\n", ri, F_e(e,ri), ri*exp(-ri));
  }

  vector_free (epsilon);
  vector_free (fx);
  fclose (fp2);

}
