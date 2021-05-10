#define RANDOM ( (double)rand()/RAND_MAX)

/* ............................................................
 * dim: dimensions of integration space
 * a,b: boundary of integration space
 * f  : integrand of vector x
 * N  : number of sample points
 * ............................................................ */
void
strata (
    int dim,
    double f(int dim, double* x),
    double* a,
    double* b,
    double abstol,
    double reltol,
    int N,
    int n_reuse,
    double mean_reuse,
    double* result,
    double* error
    )
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
    for (int j=0; j<dim; j++)
    {
      x[j] = a[j] + RANDOM*(b[j] - a[j]); // Inbetween edges
    }
    // Evaluation of f at sample point
    double fx = f(dim, x);

    // Setup average and variation of f for result and error
    sum +=fx;
    sum2+=fx*fx;
  }
  double avg = sum/N, stddev = sqrt (sum2/N - avg*avg);

  *result = V*avg;
  *error  = V*stddev/sqrt(N);

  double tolerance = abstol + reltol * fabs(*result);
  if ( *error <  tolerance )
  { return; }
  else
  {
    // For each dimension: subdivide volume in two along the dimension
    // estimate sub-variances in the two sub-volumnes
    // and pick dimension with the largest sub-variance
    // This dimension is to be subdivided again;
    // dispatch two recursive calls to each of the sub-volumes;
    // and estimate the grand average and grand error.

    // Left / right handside of subdivision
    int N2 = 16*dim;
    int n_l[dim], n_r[dim];
    double x[dim], mean_l[dim], mean_r[dim], mean = 0;

    // Reset to default setting of left/right partition
    for (int k=0; k<dim; k++)
    {
      mean_l[k] =0; mean_r[k]=0; n_l[k]=0; n_r[k]=0;
    }

    // Set random sample point
    for (int i=0; i<N2; i++)// N: sample points
    {
      for (int k=0; k<dim; k++)
      {
        x[k] = a[k] + RANDOM*(b[k] - a[k]);
      }

      double fx = f(dim, x);

      // Do the subdivision
      for (int i=0; i<dim; i++)
      {
        // In which partition does the point belong?
        if (x[i] > (a[i] + b[i])/2.){ n_r[i]++; mean_r[i] +=fx; }
        else                        { n_l[i]++; mean_l[i] +=fx; }
      }
      mean += fx;
    }
    mean /= N2;

    // Calculate the mean values of left and right partition
    for (int k=0; k<dim; k++)
    {
      mean_l[k] /= n_l[k];
      mean_r[k] /= n_r[k];
    }


    // Determine which index has the greatest variance
    int index_division=0; double maxvariance=0;
    for (int k=0; k<dim; k++)
    {
      double variance = fabs (mean_r[k] - mean_l[k]);
      if (variance > maxvariance)
      {
        maxvariance = variance;
        index_division = k;
      }
    }

    double res = V*(mean*N2 + mean_reuse*n_reuse) / (N2 + n_reuse);
    double err = fabs ( mean_reuse - mean )*V;
    tolerance = abstol + fabs(res)*reltol;

    if (err<tolerance){*result = res; *error = err; return; }

    double a2[dim], b2[dim];
    for (int k=0; k<dim; k++) { a2[k] = a[k]; }
    for (int k=0; k<dim; k++) { b2[k] = b[k]; }

    a2[index_division] = (a[index_division] + b[index_division])/2.;
    b2[index_division] = (a[index_division] + b[index_division])/2.;

    double error_l, error_r, result_l, result_r;

    // Recursion (with reuse)
    strata(dim, f, a, b2, abstol/sqrt(2.), reltol, N2, n_l[index_division], mean_l[index_division], &result_l, &error_l);
    strata(dim, f, a2, b, abstol/sqrt(2.), reltol, N2, n_r[index_division], mean_r[index_division], &result_r, &error_r);

    *result = result_l + result_r;
    *error  = error_l  + error_r;

    return;
  }

}
