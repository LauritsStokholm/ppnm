#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

int
main (int argc, char* argv[])
{
  printf(
  "A) Build a program that reads a set of numbers from the command-line and \
prints these numbers together with function values in a table form to \
the standard output.\n");

  // Given numbers x, program will print x, cos(x), sin(x).
  assert(argc > 1);

  // i=0 is name of program, we iterate over given numbers
  printf("x\tcos(x)\tsin(x)\n");
  for (int i=1; i<argc; i++)
  {
    printf("%i\t%lg\t%lg\n",atoi(argv[i]),cos(i),sin(i));
  }

  return 0;
}
