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
task_a (void)
{
  // Find a minimum of the Rosenbrock's valley function
  // Parameters (a=1, b=100)
  my_f_params p = {1, 100};

  // Initial guess-point
  vector* x = vector_alloc (2);
  vector_set (x, 0, 1);
  vector_set (x, 0, 3);

  double eps = 1e-12;

  printf ("ROSENBROCK'S VALLEY FUNCTION:\n");
  int nsteps = qnewton (Rosenbrock, x, (void*) &p, eps);
  printf ("nsteps = %i\n", nsteps);
  printf ("The minimum found is at:\n");
  vector_printf (x);
  printf ("Checking f(x,y)=%lg\n\n\n", Rosenbrock(x, (void*) &p) );


  // Initial guess-point
  vector_set (x, 0, 1);
  vector_set (x, 0, 3);

  printf ("HIMMELBLAU'S FUNCTION\n");
  nsteps = qnewton (Himmelblau, x, NULL, eps);
  printf ("nsteps = %i\n", nsteps);
  printf ("The minimum found is at:\n");
  vector_printf (x);
  printf ("Checking f(x,y)=%lg\n", Himmelblau(x, NULL) );



  vector_free (x);
}

