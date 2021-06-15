#include "main.h"    // Header functions
#include "library.c" // Vector library
#include "nmsimplex.c"
#include "qnewton.c"
#include "adapt.c"

#include "ann.h"
#include "ann.c"
#include "activation_functions.c"
#include "my_ode.c"

// Setting Supervised learning data
double a = 0, b = 2*M_PI;  // Domain of interest (interval

// Boundary condition for ODE
// Let c = 0, yc = 1, dyc = 0
double c = 0, yc = 1, dyc = 0;    // Boundary Conditions

int
main (int argc, char **argv)
{

 /* INITIALISATION */
  int n = 20; // Number of hidden-neurons
  double eps = 1e-12; // precision of minimisation
  ann* ann_gauss = ann_alloc (n, &gaussian, &gaussian_derivative, &gaussian_2derivative);
  int nx = 60; // grid number (partition)

  // Data for learning
  vector* xvals = vector_alloc (nx);

  for (int i=0; i<nx; i++)
  {
    double x = a + (b - a)*i/(nx-1); // (nx-1) as i=0 untill i=nx-1
    vector_set (xvals, i, x);
  }

  // Initialisation of parameters for ANN stucture
  ann_init_p (ann_gauss, a, b);
  ann_supervised_train (ann_gauss, xvals, a, b, c, yc, dyc,  &costfunction, eps);

  // iterate over data for plotting
  //FILE* fp = fopen ("datapoints.txt", "w");
  //for (int i=0; i<xvals->size; i++)
  //{
  //  double x = vector_get (xvals, i);
  //  double y = vector_get (yvals, i);
  //  double dy = function_to_fit_deriv(x);
  //  double Y  = function_to_fit_integ(x, xmin);
  //  fprintf (fp, "%lg\t%lg\t%lg\t%lg\n", x, y, dy, Y);
  //}

  //FILE* gp = fopen ("fit.txt", "w");
  //double dx = 1./64;
  //for (double x=xmin; x<=xmax; x+=dx)
  //{
  //  double y  = ann_gauss->f (x, (void*) ann_gauss->parameters);
  //  double dy = ann_gauss->df(x, (void*) ann_gauss->parameters);
  //  double Y  = ann_gauss->F (xmin, x, 1e-6, 1e-6, (void*) ann_gauss->parameters);
  //  fprintf (gp, "%lg\t%lg\t%lg\t%lg\n", x, y, dy, Y);
  //}

  // Free memory
  //fclose (fp); fclose (gp);
  vector_free (xvals); vector_free (yvals);
  ann_free (ann_gauss);

  return 0;
}

