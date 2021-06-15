int equal(double a, double b, double tau, double epsilon){

  /* .......................................................
   * As part of exercise 3: Returns 1 if a and b are equal with
   * absolute precision tau or relative precision epsilon:
   *
   * |a-b|<tau
   * or 
   * |a-b|/|a+b| < epsilon/2
   * ....................................................... */
  double x = fabs(a-b);
  double y = fabs(a+b);
  double absolute_precision = x;
  double relative_precision = x/y;

  //printf("The absolute precision is %lg\n", absolute_precision);
  //printf("The relative precision is %lg\n", relative_precision);

  if(absolute_precision<tau)
    return 1;
  if(relative_precision<epsilon)
    return 1;
  else
    return 0;
}

int 
equality(void)
{
  /* .......................................................
   * Write a function with the signature; that returns 1 if the numbers a and b
   * are equal, iwth absolute precision 'tau' (|a-b|<tau).
   * Or similarly, with relative precision 'epsilon' (|a-b|/|a+b| < eps/2)
   * and returns 0 otherwise. 
   * ....................................................... */
  printf("We test our equal function:\n");

  double x = M_PI;
  double y = M_PI+DBL_EPSILON;
  double z = 3;

  double l = M_E;
  double m = M_E-DBL_EPSILON;
  double k = 2;

  double f = 9.82;
  double g = 10;
  double h = M_PI*M_PI;

  double tau0 = DBL_EPSILON;
  double tau1 = FLT_EPSILON;
  double tau2 = 1;

  double epsilon0 = tau0;
  double epsilon1 = tau1;
  double epsilon2 = tau2;


 /* PI */
  int axy = equal(x, y, tau0, epsilon0);
  int bxy = equal(x, y, tau1, epsilon1);
  int cxy = equal(x, y, tau2, epsilon2);

  int axz = equal(x, z, tau0, epsilon0);
  int bxz = equal(x, z, tau1, epsilon1);
  int cxz = equal(x, z, tau2, epsilon2);

 /* EULERS NUMBER */
  int alm = equal(l, m, tau0, epsilon0);
  int blm = equal(l, m, tau1, epsilon1);
  int clm = equal(l, m, tau2, epsilon2);

  int alk = equal(l, k, tau0, epsilon0);
  int blk = equal(l, k, tau1, epsilon1);
  int clk = equal(l, k, tau2, epsilon2);

 /* GRAVITATIONAL CONSTAN */
  int afg = equal(f, g, tau0, epsilon0);
  int bfg = equal(f, g, tau1, epsilon1);
  int cfg = equal(f, g, tau2, epsilon2);

  int afh = equal(f, h, tau0, epsilon0);
  int bfh = equal(f, h, tau1, epsilon1);
  int cfh = equal(f, h, tau2, epsilon2);

  printf("To test any equality, we test for the absolute and relative precisions\
that is, if |a-b|<tau and |a-b|/|a+b| < epsilon/2\n\n");
  printf("We have chosen tau and epsilon of three values (eps_i = tau_i), so:\n");
  printf("tau0=epsilon0=%lg\t\t(DBL_EPSILON)\n", tau0);
  printf("tau1=epsilon0=%lg\t\t(FLT_EPSILON)\n", tau1);
  printf("tau2=epsilon0=%lg\n", tau2);

  printf("Subjects for the tests are the following constants; PI, EULERS NUMBER and GRAVITATIONAL ACCELERATION (at earth)\n");
 /* PI */
  printf("PI, PI+DBL_EPSILON and 3\n");
  printf("Let x=%lg, y=%lg, z=%lg\n", x, y, z);
  printf("Is x=y? (by varying precision of 0, 1, 2)\n");
  printf("%i\n%i\n%i\n", axy, bxy, cxy);
  printf("Is x=z? (by varying precision of 0, 1, 2)\n");
  printf("%i\n%i\n%i\n", axz, bxz, cxz);

 /* EULERS NUMBER */
  printf("EULERS_NUMBER (E), E+DBL_EPSILON and 2\n");
  printf("Let l=%lg, m=%lg, k=%lg\n", l, m, k);
  printf("Is l=m? (by varying precision of 0, 1, 2)\n");
  printf("%i\n%i\n%i\n", alm, blm, clm);
  printf("Is l=k? (by varying precision of 0, 1, 2)\n");
  printf("%i\n%i\n%i\n", alk, blk, clk);

 /* GRAVITATIONAL NUMBER */
  printf("GRAVITATIONAL CONSTANT (G) G+DBL_EPSILON and PI*PI\n");
  printf("Let f=%lg, g=%lg, h=%lg \n", f, g, h);
  printf("Is f=g? (by varying precision of 0, 1, 2)\n");
  printf("%i\n%i\n%i\n", afg, bfg, cfg);
  printf("Is f=h? (by varying precision of 0, 1, 2)\n");
  printf("%i\n%i\n%i\n", afh, bfh, cfh);

  printf("This is as expected, so we continue..\n");
  printf("\n\n");
  return 0;
}
