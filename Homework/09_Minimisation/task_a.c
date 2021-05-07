struct my_f_params {double a; double b;};
typedef struct my_f_params my_f_params;

/* Rosenbrock's valley function */
double
Rosenbrock (vector* z, void* params)
{
  // Unpack variables
  double x = vector_get (z, 0);
  double y = vector_get (z, 1);

  // Unpack parameters
  my_f_params* p = (my_f_params*) params;
  double a = p->a;
  double b = p->b;

  //return (a-x)*(a-x) + b*(y-x*x)*(y-x*x);
  return pow(a-x, 2) + b*pow(y - x*x, 2);
}

/* Himmelblau's function */
double
Himmelblau (vector* z, void* params)
{
  double x = vector_get (z, 0);
  double y = vector_get (z, 1);

  return (x*x + y - 11)*(x*x + y - 11) + (x + y*y - 7)*(x + y*y - 7);
}

void
fprintf_results (FILE* fp, int nsteps, vector* x, void* params,
    double f(vector* x, void* params))
{
  fprintf (fp, "nsteps = %i\n", nsteps);
  fprintf (fp, "The minimum found is at:\n");
  vector_fprintf (fp, x);
  fprintf (fp, "Checking f(x,y)=%lg\n", f(x, params));
}


void
task_a (void)
{
  FILE* fp = fopen ("output.txt", "w");
  int nsteps;

  // Find a minimum of the Rosenbrock's valley function
  // Parameters (a=1, b=100)
  my_f_params p = {1, 100};

  // Initial guess-point
  vector* x = vector_alloc (2);
  vector_set (x, 0, 1); vector_set (x, 1, 3);
  double eps = 1e-12;
  double simplex_step_goal = eps;

  fprintf (fp, "\nROSENBROCK'S VALLEY FUNCTION:\n");
  fprintf (fp, "Using qnewton\n");
  nsteps = qnewton (Rosenbrock, x, (void*) &p, eps);
  fprintf_results  (fp, nsteps, x, (void*) &p, Rosenbrock);

  fprintf (fp, "\nUsing nmsimplex\n");
  vector* start = vector_alloc (2);
  vector* step  = vector_alloc (2);

  vector_set (start, 0, 1); vector_set (start, 1, 3);
  vector_set (step,  0, 1); vector_set (step,  1, 1);

  nsteps = nmsimplex (Rosenbrock, start, step, (void*) &p, simplex_step_goal);
  fprintf_results    (fp, nsteps, start, (void*) &p, Rosenbrock);


  /*-------------------------------------------------------------------------*/
  // reset initial guess (starting point)
  vector_set (x, 0, 1); vector_set (x, 1, 3);

  fprintf (fp, "\n\n\nHIMMELBLAU'S FUNCTION\n");
  fprintf (fp, "Using qnewton\n");
  nsteps = qnewton (Himmelblau, x, NULL, eps);
  fprintf_results  (fp, nsteps, x, NULL, Himmelblau);

  // Reset initial guess-point
  vector_set (start, 0, 1); vector_set (start, 1, 3);
  vector_set (step, 0, 1) ; vector_set (step, 1, 1);
  fprintf (fp, "\nUsing nmsimplex\n");

  nsteps = nmsimplex (Himmelblau, start, step, NULL, simplex_step_goal);
  fprintf_results    (fp, nsteps, start, NULL, Himmelblau);

  // Free memory and close file pointer
  vector_free (x); vector_free (start); vector_free (step);
  fclose (fp);
}

