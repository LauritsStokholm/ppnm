#include "main.h"
#include "1_intmax.h"
#include "2_epsilon.h"
#include "3_sum.h"
#include "4_equal.h"
#include "5_digits.h"

int
main (int argc, char **argv)
{
  printf ("Task 1: Using while/for/dowhile loops for INT_MIN and INT_MAX.\n");
  intmax();
  printf ("Task 2: The machine epsilon comparison of FLOAT, DOUBLE, LONGDOUBLE.\n");
  epsilon();
  printf ("Task 3: Consider the harmonic series of (1/n); compare numerically by summing up or summing down.\n");
  sum();
  printf ("Task 4: Testing double precision constants for equality.\n");
  equality();
  printf ("Task 5: Testing the name_digit(int i) function.\n");
  digits();
  return 0;
}

