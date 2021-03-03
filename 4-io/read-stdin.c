#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

int
main (int argc, char **argv)
{
  printf("B) Build a program that reads a set of numbers from the standard \
input and prints these numbers together with function values in a table form \
to the standard output\n");
  double x;
  printf("x\tcos(x)\tsin(x)\n");
  while(scanf("%lg", &x) != EOF )
  {
    printf("%lg\t%lg\t%lg\n",x,cos(x),sin(x));
  }
  return 0;
}

