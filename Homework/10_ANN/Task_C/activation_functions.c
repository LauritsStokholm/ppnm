#ifndef HAVE_ACTIVATION_FUNCTIONS
#define HAVE_ACTIVATION_FUNCTIONS


// ODE COST-FUNCTION
double
costfunction (vector* annparams, void* params)
{
  // my_struct_params = (ann* network, vector* xdata, vector* ydata)
  my_struct_params* p = (my_struct_params*) params;

  ann* network            = p->ann;
  vector* xdata           = p->x;
  double a                = p->a;
  double b                = p->b;
  double c                = p->c;
  double yc               = p->yc;
  double dyc              = p->dyc;
  double (*phi)(double y) = p->phii;

  assert (xdata->size == ydata->size);
  double sum = 0, xk, yk, Fpxk;

  int n = xdata->size;
  for (int k=0; k<n; k++)
  {
    xk = vector_get (xdata, k);

    double func  = network->f   (xk, (void*) annparams);
    double diff  = network->df  (xk, (void*) annparams);
    double ddiff = network->ddf (xk, (void*) annparams);

    sum += pow ( phi(x, func, diff, ddiff   ), 2);
  }
  sum *= (b-a) / (n-1);
  sum += pow (network->f (c, (void*) annparams) - yc , 2);
  sum += pow (network->df (c, (void*) annparams)- dyc, 2);

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

double gaussian_2derivative (double x, void* params)
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
    sum += (-2/bi)*wi*exp(-xi*xi) + 4*xi*xi/bi/bi*wi*exp(-xi*xi);
  }
  return sum;
}
