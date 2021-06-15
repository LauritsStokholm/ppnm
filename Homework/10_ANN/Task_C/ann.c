#ifndef HAVE_ANN_C
#define HAVE_ANN_C
#include "ann.h"
#include "activation_functions.c"

/*----------------------------------------------------------------------------*/
ann*
ann_alloc
  (
    int n,
    double (*f)   (double, void* params),
    double (*df)  (double, void* params),
    double (*ddf) (double, void* params)
  )
{
  ann* network = malloc ( sizeof(ann) );
  network -> nneurons = n;
  network -> parameters = vector_alloc (3*n);
  network -> f   = f  ;
  network -> df  = df ;
  network -> ddf = ddf;
  return network;
}

void
ann_free (ann* network)
{
  vector_free (network->parameters);
  free (network);
}

void
ann_init_p (ann* network, int xmin, int xmax)
{
  int n = network -> nneurons;
  vector* params = network -> parameters;

  // Three parameters per neuron
  for (int i=0; i<n; i++)
  {
    double ai = xmin + (xmax-xmin)*i/(n-1);
    double bi = 1, wi = 1;

    vector_set (params, 3*i,   ai);
    vector_set (params, 3*i+1, bi);
    vector_set (params, 3*i+2, wi);
  }
}

int
ann_train_nmsimplex
  (
    ann* network,
    double (*cost) (vector* x, void* params),
    void* params, // for the cost function
    double simplex_step_goal
  )
{
  // Assume initialised
  vector* start = vector_alloc (3*network->nneurons);
  vector* step  = vector_alloc (3*network->nneurons);

  // stepsize: (instead of unity)
  for (int i=0; i<(3*network->nneurons); i++)
  {
    vector_set (start, i, vector_get (network->parameters, i) );
    vector_set (step,  i, vector_get (network->parameters, i) / 5.);
  }
  int nsteps = nmsimplex (cost, start, step, (void*) params, simplex_step_goal);
  vector_memcpy (start, network->parameters);

  vector_free (start);
  vector_free (step);
  return 0;
}


void
ann_supervised_train
(
    ann* network,
    vector* xdata,
    double a,
    double b,
    double c,
    double yc,
    double dyc,
    double (*phi) (double y),
    double (*cost) (vector* x, void* params),
    double eps
)
{

  // For supervised learning, costfunction is deviation from data
  my_struct_params costparameters;
  costparameters.ann   = network;
  costparameters.x     = xdata;
  costparameters.a     = a;
  costparameters.b     = b;
  costparameters.c     = c;
  costparameters.yc    = yc;
  costparameters.dyc   = dyc;
  costparameters.phi   = phi;

  qnewton (cost, network->parameters, (void*) &costparameters, eps);
  //ann_train_nmsimplex (network, cost, (void*) &costparameters, eps);
}

double
ann_response (ann* network, double x)
{
  vector* params = network -> parameters;
  network->f (x, (void*) params);
}

#endif
