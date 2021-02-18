#include "main.h"

int
main (int argc, char **argv)
{
  //Initialise variables
  double x, y, zero;
  complex z, zz;

/* Exercise 1: Mathematical functions */
  printf("Exercise 1: Evaluate following mathematical functions\n");

  printf("%-20s %32s %32s\n", "Function", "Argument", "Value");
  printf("%-32s %16s|%-16s %16s|%-16s\n", "", "Real", "Imag", "Real", "Imag");

  x = 5; zero = 0;
  y = tgamma(x);
  printf("%-32s %16lg|%-16lg %16lg|%-16lg \n", "True gamma (tgamma)", x, zero , y, zero);

  x = 0.5;
  y = j1(x);
  printf("%-32s %16lg|%-16lg %16lg|%-16lg\n", "Bessel 1st kind 1st order (j1) ", x, zero , y, zero);

  x = -2;
  z = csqrt(x);
  printf("%-32s %16lg|%-16lg %16lg|%-16lg\n", "Complex square root", x, zero, creal(z), cimag(z));

  z = I*M_PI;
  zz = cexp(z);
  printf("%-32s %16lg|%-16lg %16lg|%-16lg\n", "Complex exponential", creal(z), cimag(z), creal(zz), cimag(zz));

  z = I;
  zz = cexp(z);
  printf("%-32s %16lg|%-16lg %16lg|%-16lg\n", "Complex exponential", creal(z), cimag(z), creal(zz), cimag(zz));

  zz = cpow(z, M_E);
  printf("%-32s %16lg|%-16lg %16lg|%-16lg\n", "Complex power (base i)", M_E, zero, creal(zz), cimag(zz));

  zz = cpow(z, z);
  printf("%-32s %16lg|%-16lg %16lg|%-16lg\n", "Complex power (base i)", creal(I), cimag(I), creal(zz), cimag(zz));

  printf("\n");

/* Exercise 2: Significant digits */
  float x_g;
  long double x_Lg;

  printf("Exercise 2: Here we test the significant digits of variable types\n");
  x_g = 1.f/9;
  x = 1./9;
  x_Lg = 1.L/9;

  printf("Float:       %.25g\n", x_g);
  printf("Double:      %.25lg\n", x);
  printf("Long double: %.25Lg\n", x_Lg);
  printf("              .123456789|123456789|123456789|\n");


  return 0;
}

