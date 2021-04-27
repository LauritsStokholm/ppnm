void
recursive_part_with_reuse (double f(double), double a, double b,
                double abs, double rel, double f2, double f3,
                double* result, double* error)
{
  double f1 = f(a +   (b-a)/6);
  double f4 = f(a + 5*(b-a)/6);
  double Q = ((2*f1 + f2 + f3 + 2*f4)/6.) * (b-a);
  double q = ((  f1 + f2 + f3 + f4  )/4.) * (b-a);

  double tol = abs + rel * fabs(Q);
  double err = fabs(Q-q);

  if (err < tol) {*result = Q; *error = err; return;}
  else {
    double Q1, Q2, err1, err2;
    recursive_part_with_reuse (f, a, (a+b)/2., abs/sqrt(2.), rel, f1, f2, &Q1, &err1);
    recursive_part_with_reuse (f, (a+b)/2., b, abs/sqrt(2.), rel, f3, f4, &Q2, &err2);
    *result=Q1+Q2; *error=err1+err2; return;
  }
}

void
adaptive_integrator (double f(double), double a, double b,
                               double abs, double rel,
                               double* result, double* error)
{
  double f2 = f(a + 2*(b-a)/6);
  double f3 = f(a + 4*(b-a)/6);
  recursive_part_with_reuse (f, a, b, abs, rel, f2, f3, result, error);
}


// Clenshaw curtis part
static double A, B; // In view for this file only

// Auxilliary function
double F (double f(double), double t)
{
  return f( (A+B)/2. + ((A-B)/2.)*cos(t) )*sin(t)*(B-A)/2.;
}

void recursive_clenshaw (double f(double), double a, double b,
                         double abs, double rel, double f2, double f3,
                         double* result, double* error)
{
  double f1 = F(f, a +   (b-a)/6.);
  double f4 = F(f, a + 5*(b-a)/6.);
  double Q  = ((2*f1+f2+f3+2*f4)/6.)*(b-a);
  double q  = ((f1+f4+f2+f3)/4.) * (b-a);

  double tol = abs + rel*fabs(Q);
  double err = fabs(Q-q)/3.;

  if (err < tol) {*result = Q; *error = err; return;}
  else {
    double Q1, Q2, err1, err2;
    recursive_clenshaw (f, a, (a+b)/2., abs/sqrt(2), rel, f1, f2, &Q1, &err1);
    recursive_clenshaw (f, (a+b)/2., b, abs/sqrt(2), rel, f3, f4, &Q2, &err2);

    *result = Q1+Q2;
    *error  = err1+err2;
    return;
  }
}

void
clenshaw_curtis (double f(double), double a, double b,
    double abs, double rel, double* result, double* error)
{
  // Changing variables
  A = a, B=b;
  a=0, b=M_PI;

  double f2 = F(f, a+2*(b-a)/6.);
  double f3 = F(f, a+4*(b-a)/6.);

  recursive_clenshaw (f, a, b, abs, rel, f2, f3, result, error);
  return;
}

//////////////////////////////////////////////////////////////////////////////
// Implementation of infinite limits
static double (*inf_func) (double);

double inf_inf_wrapper (double x)
{// from -1 to +1
  return inf_func (x/(1-x*x)) * (1+x*x)/((1-x*x)*(1-x*x));
}

double a_inf_wrapper (double x)
{
  return inf_func (A + x/(1-x)) / ((1-x)*(1-x));
}

double inf_b_wrapper (double x)
{
  return inf_func (B + x/(1+x)) / ((1+x)*(1+x));
}

void adapt_inf (double f(double), double a, double b, double abs, double rel,
    double* result, double* error)
{
  inf_func = f;

  if( (-1)*isinf(a) && isinf(b) )
  {
    clenshaw_curtis(inf_inf_wrapper, -1, 1, abs, rel, result, error);
    return;
  }
  else if (isinf(b))
  {
    clenshaw_curtis(a_inf_wrapper, 0, 1, abs, rel, result, error);
    return;
  }
  else if ( (-1)*isinf(a) )
  {
    clenshaw_curtis(inf_b_wrapper, -1, 0, abs, rel, result, error);
    return;
  }
  clenshaw_curtis(f, a, b, abs, rel, result, error);
  return;
}

//////////////////////////////////////////////////////////////////////////////

