#include <stdio.h>
void
generate_data (void)
{
  FILE* datafile = fopen("Rutherford.data", "w");
  double time[] = {1, 2, 3, 4, 6, 9, 10, 13, 15};
  double y[] = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};

  fprintf (datafile, "t\ty\tdy\n");
  for (int i=0; i<9; i++)
  {
    fprintf (datafile, "%lg\t%lg\t%lg\n", time[i], y[i], y[i]/20.);
  }
  fclose(datafile);
}

int
main (int argc, char* argv[])
{
  generate_data();
  return 0;
}
