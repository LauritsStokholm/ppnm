// Stratified Sampling Monte Carlo (not used)

// More evenly than random sample points, so error is not given by variance.
// In reality, error is not trivial; we estimate it by the difference of the 
// result, given by two different sequences,
//

#define fracl(x) ((x)-floorl(x))
#define PI 3.1415926535897932384626433832795028841971693993751L
void lattice (int d, double* x)
{
  static int dim=0, n=-1;
  static long double *alpha;
  int i;

  if (d<0){ /* d<0 signals to (re-)initialise the lattice */
    dim=-d; n=0; alpha=(long double*)realloc(alpha, dim*sizeof(long double));
    for (i = 0; i<dim; i++)
    {
      alpha[i] = fracl(sqrtl(PI+i));
    }
  }
  else if (d>0)
  {
    n++; assert(d==dim && n>0);
    for (i=0; i<dim; i++)
    {
      x[i] = fracl(n*alpha[i]);
    }
  }
  else if (alpha!=NULL)
  {
    free (alpha);
  }
  return;
}

void latticemc (int dim, double f(int dim, double* x), double* a, double* b,
    int N, double* result, double* error)
{
  // Setup sample points
  double vol = 1, sum = 0, sum2 = 0, x[dim];

  // Calculate the volume of the integration space
  for (int i=0; i<dim; i++)
  {
    vol *= b[i]-a[i];
  }

  for (int i=0; i<N; i++)//Up to N sample points
  {
    lattice (i, x); // This sets all elements of x vector
    double fx = f(dim, x);
    sum  += fx;
    sum2 += fx*fx;
  }
  double avg = sum/N, stddev = sqrt (sum2/N - avg*avg);
  *result = vol*avg;
  *error = vol*stddev/sqrt(N);
  return;

}


