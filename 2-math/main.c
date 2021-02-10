#include "main.h"

int
main (int argc, char **argv)
{
  //Initialise variables
  double x, y;
  complex z, zz;

/* Exercise 1: Mathematical functions */
  // Exercise 1.1: True gamma function
  x = 5;
  y = tgamma(x);
  printf("Exercise 1.1: True gamma function; gamma(%lg) = %lg\n", x, y);

  // Exercise 1.2: Bessel function
  x = 0.5;
  y = j1(x);
  printf("Exercise 1.2: Bessel function of first kind; j1(%lg) = %lg\n", x, y);

  // Exercise 1.3: Complex square root function
  x = -2;
  z = csqrt(x);
  printf("Exercise 1.3: The complex square root; sqrt(%lg) = %lg + i %lg\n", x, creal(z), cimag(z));

  // Exercise 1.4: Complex exponential function
  z = I*M_PI;
  zz = cexp(z);
  printf("Exercise 1.4: The complex exponential; exp(%lg+i%lg) = %lg + i %lg\n", creal(z), cimag(z), creal(zz), cimag(zz));

  // Exercise 1.5: Complex exponential function
  z = I;
  zz = cexp(z);
  printf("Exercise 1.5: The complex exponential; exp(%lg+i%lg) = %lg = i %lg\n", creal(z), cimag(z), creal(zz), cimag(zz));

  // Exercise 1.6: Complex power function
  zz = cpow(z, M_E);
  printf("Exercise 1.6: The complex power; cpow(%lg+i%lg, %lg) = %lg = i %lg\n", creal(z), cimag(z), M_E, creal(zz), cimag(zz));

  // Exercise 1.7: Complex power function
  zz = cpow(z, z);
  printf("Exercise 1.7: The complex power; cpow(%lg+i%lg, %lg+i%lg) = %lg = i %lg\n", creal(z), cimag(z), creal(z), cimag(z), creal(zz), cimag(zz));

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

