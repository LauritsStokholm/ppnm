
void
gslInterp (int n, double* x, double* y, double resolution, FILE* output, const gsl_interp_type* method)
{

  gsl_interp* T = gsl_interp_alloc((const gsl_interp_type*) method, n);
  gsl_interp_init (T, x, y, n);

  int size_z = (x[n-1] - x[0])/resolution;

  double* z = malloc(sizeof(double)*size_z);
  gsl_interp_accel* acc = gsl_interp_accel_alloc();

  fprintf(output, "x\tinterp\tderiv\tinteg\n");
  for (double zval=x[0]; zval<x[n-1]; zval+=resolution)
  {
    //int idx = gsl_interp_bsearch (x, zval, x[0], x[n-1]);
    double w = gsl_interp_eval (T, x, y, zval, acc);
    double w_deriv = gsl_interp_eval_deriv (T, x, y, zval, acc);
    double w_integ = gsl_interp_eval_integ (T, x, y, x[0], zval, acc);

    fprintf(output, "%lg\t%lg\t%lg\t%lg\n", zval, w, w_deriv, w_integ);
  }


  free(z);
  gsl_interp_free(T);
  gsl_interp_accel_free(acc);
}

