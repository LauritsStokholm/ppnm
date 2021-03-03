#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf_gamma.h>
#include <assert.h>

// Stirling approx (toggle: 0 truegamma, 1 lngamma)
double my_tgamma(double x){
  ///single precision gamma function (Gergo Nemes, from Wikipedia)
  if(x<0)return M_PI/sin(M_PI*x)/my_tgamma(1-x);
  if(x<9)return my_tgamma(x+1)/x;
  double log_gamma = (log(2*M_PI) - log(x))/2. + x*(log(x + 1./(12*x - 1./(10*x))) - 1);
  return exp(log_gamma);
}

double 
my_lgamma (double x)
{
  return (log(2*M_PI) - log(x))/2. + x*(log(x + (1./(12*x - (1./(10*x))))) - 1);
}


double complex my_cgamma(double complex z){

  if( creal(z) < 0 ){
    return M_PI / csin( M_PI*z ) / my_cgamma( 1 - z );
  }
  if( creal(z) < 9 ){
    return my_cgamma( z + 1 ) / z;
  }

  double complex lnGamma  =  z * clog( z + 1 / ( 12 * z - 1 /(z*10) ) ) - z + clog( 2 * M_PI / z ) / 2;
  double complex result   =  cexp(lnGamma);

  return result;
}


int
main (int argc, char **argv)
{
  // Initialise variables
  double x, y;
  double complex z;

  // Steamlines
  FILE* input1  = fopen("my_tgamma_in.txt", "r");
  FILE* output1 = fopen("my_tgamma_out.txt", "w");

  FILE* input2  = fopen("my_lgamma_in.txt", "r");
  FILE* output2 = fopen("my_lgamma_out.txt", "w");

  FILE* input3  = fopen("my_cgamma_in.txt", "r");
  FILE* output3 = fopen("my_cgamma_out.txt", "w");

  // Headers for plotting
  fprintf(output1, "x\tmath\tgsl\tmy\n");
  fprintf(output2, "x\tmath\tgsl\tmy\n");
  fprintf(output3, "x\ty\tz\n");

  //Calculations
  while(fscanf(input1, "%lg",&x)!=EOF)
  {
    fprintf(output1, "%lg\t%lg\t%lg\t%lg\n", x, tgamma(x), gsl_sf_gamma(x), my_tgamma(x));
  }
  while(fscanf(input2, "%lg",&x)!=EOF)
  {
    fprintf(output2, "%lg\t%lg\t%lg\t%lg\n", x, lgamma(x), gsl_sf_lngamma(x), my_lgamma(x));
  }

  while(fscanf(input3, "%lg\t%lg", &x, &y) != EOF)
  {
    z = x + y*I;
    fprintf(output3, "%lg\t%lg\t%lg\n", x, y, fmin( 6, cabs(my_cgamma(z) ) ));
  }

  fclose(input1);
  fclose(input2);
  fclose(input3);
  fclose(output1);
  fclose(output2);
  fclose(output3);


  return 0;
}

