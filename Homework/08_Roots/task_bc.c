 /* find the solution F_e (r) parametrised by energy e,
 * This will in general not satisfy BDC at r=0 and r->infty
 * 
 * Fix rmax = 10 (approx infty), search the roots of varying energy
 * e, of F(e, rmax)
 * ............................................................ */


int converged_iteration = 0;

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
  double rmin = 1e-3; // rmax= 10, but integration is up to variable r
  if (r < rmin){ return r - r*r;}

  double t = rmin; // Parameter for the system; to be driven from t=rmin: t=r
  double step_size = 1e-1;// Starting step
  double epsabs = 1e-3, epsrel=1e-3;  // Precision / tolerance

  void* params = (void*) &e;

  // Set y(a)
  vector* y = vector_alloc (2);
  vector_set (y, 0, t - t*t);
  vector_set (y, 1, 1 - 2*t);

  int ODE_method = 12; //rkstep12 = 12 or rkstep23 = 23

  // Integration (solve from y(a) to y(b))
  int iter = my_ode_driver (swave, params, y, rmin, r, step_size, ODE_method,
                 epsabs, epsrel);

  // For convergence test
  converged_iteration = iter;

  // To free vector and still return value.
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
  assert ( epsilon < 0 ); // Energy is negative (bound state)
  double rmax = *(double *) params;
  double M_e = F_e (epsilon, rmax);
  vector_set (fx, 0, M_e);
}

void
M_e_nonzero (vector* x, void* params, vector* fx)
{
  double epsilon = vector_get (x, 0);
  assert (epsilon < 0);
  double rmax = *(double *) params;
  double M_e = F_e (epsilon, rmax) - rmax*exp(-sqrt(-2*epsilon)*rmax);
  vector_set (fx, 0, M_e);
}

/* ............................................................
 * toggle: task b or task c (M_e or M_e_nonzero)
 * rmax: variable
 * ............................................................ */
int
caller (int toggle, double rmax)
{
  FILE* fp;
  assert (toggle == 0 || toggle == 1);
  if (toggle == 0) {fp = fopen("task_b.dat", "w");}
  if (toggle == 1) {fp = fopen("task_c.dat", "w");}


  // Initialisation of coordinate and function vectors
  vector* epsilon  = vector_alloc (1);
  vector* fx = vector_alloc (1); // M_e = F_e (eps, rmax);

  // Function parameters and initial guess of root coordinates
  double epsilon_guess = -1;
  vector_set (epsilon, 0, epsilon_guess);

  // Acceptance tolerance
  double acc = 1e-7;


  void (*f) (vector* x, void* p, vector* z);
  // Change epsilon, as to set M_e (rmax) = 0
  if (toggle == 0){f = &M_e;}
  if (toggle == 1){f = &M_e_nonzero;}

  newton (f, epsilon, (void*) &rmax, acc);


  fprintf (fp, "Checking the given root\n");
  f (epsilon, (void*) &rmax, fx);
  vector_fprintf (fp, fx);

  double e = vector_get (epsilon, 0);
  FILE* fp2 = fopen ("task2.txt", "w");
  fprintf (fp2, "rmax\tepsilon\n");
  fprintf (fp2, "%lg\t%lg\n", rmax, e);

  fprintf(fp, "r\tf(e,r)\tExact\n");
  for (double ri=0; ri<=rmax; ri+=0.001){
    fprintf(fp, "%lg\t%lg\t%lg\n", ri, F_e(e, ri), ri*exp(-ri));
  }

  vector_free (epsilon);
  vector_free (fx);
  fclose (fp);
  fclose (fp2);

  return converged_iteration;
}

// Convergence test
int
convergence (int toggle)
{

  assert (toggle == 0 || toggle == 1);
  FILE* fp;
  if (toggle == 0) {fp = fopen("convergence.txt", "w");}
  else { fp = fopen ("precision.txt", "w");}

  fprintf (fp, "rmax\titer\n");
  for (double rmax=2; rmax < 10; rmax+=0.1)
  {
    int iter = caller (toggle, rmax);
    fprintf (fp, "%lg\t%i\n", rmax, iter);
  }
  return 0;
}


void
task_b (void)
{
  convergence (0);
}


void
task_c(void)
{
  convergence (1);
}

