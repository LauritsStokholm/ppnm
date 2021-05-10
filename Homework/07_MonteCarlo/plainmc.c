#define RANDOM ( (double)rand()/RAND_MAX)

/* ............................................................
 * dim: dimensions of integration space
 * a,b: boundary of integration space
 * f  : integrand of vector x
 * N  : number of sample points
 * ............................................................ */
void
plainmc (int dim, double f(int dim, double* x), double* a, double* b, int N,
    double* result, double* error)
{
  // Calculate the volume of integration space
  double V = 1; 
  for (int i=0; i<dim; i++)
  {
    V*=b[i]-a[i];
  }

  // Setup sample points, calculate avg of f(xi)
  double sum = 0, sum2 = 0, x[dim];
  for (int i=0; i<N; i++)
  {
    // Calculate i'th element of sample point in 
    // dim-dimensional integration space
    for (int i=0; i<dim; i++)
    {
      x[i] = a[i] + RANDOM*(b[i] - a[i]); // Inbetween edges
    }
    // Evaluation of f at sample point
    double fx = f(dim, x);

    // Setup average and variation of f for result and error
    sum+=fx; sum2+=fx*fx;
  }
  double avg = sum/N, stddev = sqrt (sum2/N - avg*avg);
  *result = V*avg;
  *error  = V*stddev/sqrt(N);
  return;
}
