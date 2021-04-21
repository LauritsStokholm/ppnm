#include "main.h"
#include "gsl.c"
#include "linear.c"
#include "quadratic.c"
#include "cubic.c"
#include "utilities.c"

#include <math.h>


//double linterp_integ (double* x, double* y, double z)
//{
//
//  return integ;
//}


int
main (int argc, char **argv)
{
  /* ............................................................
   * Parameters:
   * datfile: file containing (x, y) data. If not specified, the
   * program will use "data.txt" as default.
   * 
   * Resolution is for interpolation step.
   * ............................................................ */

  /* Set default datafile */
  //char* datfile = "data.txt";
  double resolution = 1e-2;


  // In case you want to change data radically
  FILE* fp = fopen("data.txt", "w");
  for (double x=0; x<10; x+=1e-2){fprintf(fp, "%lg\t%lg\n", x, cos(x));}
  fclose(fp);
  char* datfile = "data.txt";
  //



  /* Check for specified parameters */
  while (1)
  {
    int opt = getopt (argc, argv, "n:r");
    if (opt == -1) break;
    switch (opt)
    {
      case 'n': datfile = optarg; break;
      case 'r': resolution = atof (optarg); break;
      default:
      fprintf (stderr, "Usage: %s [-n name of datafile] [-r resolution]\n", argv[0]);
      exit (EXIT_FAILURE);
    }
  }

  /* Prepare streamlines */
  FILE* fp_data= fopen(datfile, "r");
  FILE* fp_lout = fopen("lspline.txt", "w");
  FILE* fp_qout = fopen("qspline.txt", "w");
  FILE* fp_cout = fopen("cspline.txt", "w");
  FILE* fp_lgsl = fopen("lgsl.txt", "w");
  FILE* fp_cgsl = fopen("cgsl.txt", "w");


  /* Read datafile, save size and allocate list for storage of data */
  double xi = 0; double yi = 0;
  int num_points = 0;
  while (fscanf(fp_data, "%lg\t%lg", &xi, &yi) != EOF) {num_points++;}
  rewind(fp_data); // Necessary, as fscanf is at EOF

  printf("Number of datapoints: %i\n", num_points);
  double* x = malloc(sizeof(double)*num_points);
  double* y = malloc(sizeof(double)*num_points);

  int idx = 0;
  while (fscanf(fp_data, "%lg\t%lg", &xi, &yi) != EOF){
    x[idx] = xi;
    y[idx] = yi;
    idx++;
  }


 /* Interpolation routines  */
  linear_interp    (num_points, x, y, resolution, fp_lout);
  quadratic_interp (num_points, x, y, resolution, fp_qout);
  cubic_interp     (num_points, x, y, resolution, fp_cout);
  gslInterp (num_points, x, y, resolution, fp_lgsl, gsl_interp_linear);
  gslInterp (num_points, x, y, resolution, fp_cgsl, gsl_interp_cspline);

  free(x); free(y);
  return 0;
}
