// Rosenbrock parameters
struct my_f_params {double a; double b;};
typedef struct my_f_params my_f_params;

// Gradient
void
Rosenbrock_gradient (vector* r, void* params, vector* fx)
{
  // f(x, y) = (a-x)*(a-x) + b*(y-x*x)*(y-x*x)
  // Unpack coordinate vector
  double x = vector_get (r, 0);
  double y = vector_get (r, 1);

  // Unpack parameters
  my_f_params* p = (my_f_params*) params;
  double a = p->a;
  double b = p->b;

  // Calculate gradient function
  vector_set (fx, 0, -2*(a-x) + b*2*(y-x*x)*(-2*x) );
  vector_set (fx, 1, b*2*(y-x*x) );
}

/* ............................................................
 * We search for the extremum(s) of the Rosenbrock valley function
 * by searching for the roots of its gradient
 * ............................................................ */
void
task_a (void)
{
  FILE* fp = fopen ("task_a.txt", "w");

  // Initialisation of coordinate and function vectors
  vector* r  = vector_alloc (2);
  vector* fx = vector_alloc (2);

  // Function parameters and initial guess of root coordinates
  my_f_params params = {1, 100}; // (a, b)
  vector_set (r, 0, 1);
  vector_set (r, 1, 2);

  fprintf(fp, "Initial guess / starting point\n");
  vector_fprintf (fp, r);

  // Acceptance tolerance
  double acc = 1e-7;
  newton (Rosenbrock_gradient, r, (void*) &params, acc);

  fprintf (fp, "Converged: root is found at\n");
  vector_fprintf (fp, r);
  fprintf (fp, "Checking the given root\n");
  Rosenbrock_gradient (r, (void*) &params, fx);

  fprintf (fp, "Evaluation at the root gives:\n");
  vector_fprintf(fp, fx);

  fclose(fp);
  vector_free (r);
}
