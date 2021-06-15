#include <stdio.h>
#include <math.h>
#include <assert.h>

int
main (int argc, char **argv)
{
  assert(argc==3); // We expect ./main, input.txt, output.txt

  FILE* input  = fopen(argv[1], "r");
  FILE* output = fopen(argv[2], "w");

  fprintf(output,
"C) Build a program that reads a set of numbers from an inputfile \
and writes them together with function values in a table form to an \
outputfile. The program must read the names of the input/output files from \
the command-line and then use fopen.\n");
  double x;
  fprintf(output, "x\tcos(x)\tsin(x)\n");
  while(fscanf(input, "%lg", &x) != EOF )
  {
    fprintf(output, "%lg\t%lg\t%lg\n",x,cos(x),sin(x));
  }

  fclose(input);
  fclose(output);
  return 0;
}

