#ifndef HAVE_ACTIVATION_FUNCTIONS
#define HAVE_ACTIVATION_FUNCTIONS

double
costfunction (vector* annparams, void* params)
{
  // my_struct_params = (ann* network, vector* xdata, vector* ydata)
  my_struct_params* p = (my_struct_params*) params;

  ann* network  = p->ann;
  vector* xdata = p->x;
  vector* ydata = p->y;

  assert (xdata->size == ydata->size);
  double sum = 0, xk, yk, Fpxk;
  for (int k=0; k<xdata->size; k++)
  {
    xk = vector_get (xdata, k);
    yk = vector_get (ydata, k);

    Fpxk = network->f(xk, (void*) annparams);
    sum += (Fpxk - yk)*(Fpxk - yk) ;
  }
  printf ("sum=%lg\n", sum);
  return sum;
}

double gaussian (double x, void* params)
{
  vector* ANNparams = (vector*) params;

  int n = ANNparams->size / 3.;
  double ai, bi, wi, xi;

  double sum = 0;
  for (int i=0; i<n; i++)
  {
    ai = vector_get (ANNparams, 3*i);
    bi = vector_get (ANNparams, 3*i+1);
    wi = vector_get (ANNparams, 3*i+2);

    xi = (x - ai) / bi;
    sum += wi*exp(-xi*xi);
  }
  return sum;
}

double gaussian_derivative (double x, void* params)
{
  vector* ANNparams = (vector*) params;

  int n = ANNparams->size / 3.;
  double ai, bi, wi, xi;

  double sum = 0;
  for (int i=0; i<n; i++)
  {
    ai = vector_get (ANNparams, 3*i);
    bi = vector_get (ANNparams, 3*i+1);
    wi = vector_get (ANNparams, 3*i+2);

    xi = (x - ai) / bi;
    sum += (-2*xi/bi)*wi*exp(-xi*xi);
  }
  return sum;
}

double gaussian_integ (double xmin, double x, double epsabs, double epsrel, void* params)
{
  double result, error;
  adaptive_integrator (gaussian, params, xmin, x, epsabs, epsrel, &result, &error);
  return result;
}


double gaussian_wavelet (double x, void* params)
{
  vector* ANNparams = (vector*) params;

  int n = ANNparams->size / 3.;
  double ai, bi, wi, xi;

  double sum = 0;
  for (int i=0; i<n; i++)
  {
    ai = vector_get (ANNparams, 3*i);
    bi = vector_get (ANNparams, 3*i+1);
    wi = vector_get (ANNparams, 3*i+2);

    xi = (x - ai) / bi;
    sum += wi*xi*xi*exp(-xi*xi);
  }
  return sum;
}

double wavelet (double x, void* params)
{
  vector* ANNparams = (vector*) params;

  int n = ANNparams->size / 3.;
  double ai, bi, wi, xi;

  double sum = 0;
  for (int i=0; i<n; i++)
  {
    ai = vector_get (ANNparams, 3*i);
    bi = vector_get (ANNparams, 3*i+1);
    wi = vector_get (ANNparams, 3*i+2);

    xi = (x - ai) / bi;
    sum += wi*cos(5*xi)*xi*exp(-xi*xi);
  }
  return sum;
}


double my_wavelet (double x, void* params)
{
  vector* ANNparams = (vector*) params;

  int n = ANNparams->size / 3.;
  double ai, bi, wi, xi;

  double sum = 0;
  for (int i=0; i<n; i++)
  {
    ai = vector_get (ANNparams, 3*i);
    bi = vector_get (ANNparams, 3*i+1);
    wi = vector_get (ANNparams, 3*i+2);

    xi  = (x - ai) / bi;
    sum += wi*xi*xi*sin(xi);
  }
  return sum;
}

double my_wavelet_dif (double x, void* params)
{
  vector* ANNparams = (vector*) params;

  int n = ANNparams->size / 3.;
  double ai, bi, wi, xi, dxi;

  double sum = 0;
  for (int i=0; i<n; i++)
  {
    ai = vector_get (ANNparams, 3*i);
    bi = vector_get (ANNparams, 3*i+1);
    wi = vector_get (ANNparams, 3*i+2);

    xi  = (x - ai) / bi;
    dxi = 1./bi;
    sum += 2*wi*xi*dxi*sin(xi) + dxi*wi*xi*xi*cos(xi);
  }
  return sum;
}

double my_wavelet_int (double x, double xmin, double xmax, void* params)
{
  vector* ANNparams = (vector*) params;

  int n = ANNparams->size / 3.;
  double ai, bi, wi, xi;

  double sum = 0;
  for (int i=0; i<n; i++)
  {
    ai = vector_get (ANNparams, 3*i);
    bi = vector_get (ANNparams, 3*i+1);
    wi = vector_get (ANNparams, 3*i+2);

    xi  = (x - ai) / bi;
    sum += wi*

      wi*xi*xi*sin(xi);
  }
  return sum;
}



#endif
