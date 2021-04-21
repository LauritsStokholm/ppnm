
/* ............................................................
 *  Here are my implemenations of functions used for a
 *  linear interpolation.
 * ............................................................ */

void
linear_interp (int num_points, double* x, double* y, double resolution, FILE* fp_out)
{
  // We know x is ordered
  double xmin = x[0];
  double xmax = x[num_points-1];

  double lspline_eval, lspline_deriv, dx, dy;
  double lspline_integ = 0;
  int k;

  fprintf(fp_out, "x\teval\tderiv\tinteg\n");

  for (double w=xmin; w<=xmax; w+=resolution)
  {
  int i = binarysearch(num_points, x, w);
  dx = x[i+1] - x[i];
  dy = y[i+1] - y[i];

  lspline_eval = y[i] + (dy/dx) * (w-x[i]);
  lspline_deriv = dy/dx;

  k = 0;
  lspline_integ = y[i]*(w-x[i]) + (dy/dx)*(w-x[i])*(w-x[i])/2;
  while (k<i)
  {
    dx = x[k+1] - x[k];
    dy = y[k+1] - y[k];

    lspline_integ += y[k]*dx + (dy/dx)*dx*dx/2.;
    k++;
  }

  fprintf(fp_out, "%lg\t%lg\t%lg\t%lg\n", w, lspline_eval, lspline_deriv, lspline_integ);
  }
}



