#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double my_exp(double x){
if(x<0)return 1/my_exp(-x);
if(x>1./8)return pow(my_exp(x/2),2);
return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}

int
main (int argc, char **argv)
{
  FILE* fp = fopen ("output.txt", "w");
  double dx = 1e-3, xmin = -2, xmax = 5;
  fprintf(fp, "x\texp\tmyexp\n");
  for (double x=xmin; x<xmax; x+=dx){
    fprintf (fp, "%lg\t%lg\t%lg\n", x, exp(x), my_exp(x));
  }

  fclose(fp);


  return 0;
}

