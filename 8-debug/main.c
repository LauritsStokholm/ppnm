// I have changed quotations to angle brackets; since " " is for local
// header files, and these are standard libraries.
#include <stdio.h>
#include <gsl/gsl_matrix.h>

int print_half_00(gsl_matrix* m)
{

// Here I have changed the division, which previously was integer division
// (resulting in zero), whereas this 1./2 is double division resulting in (0.5)
  double half = 1./2.;

// The asterisk (*) was removed, since status is not a pointer to an integer,
// but just an integer. Also the ampersand (&) was removed at m, since
// gsl_matrix_get takes a gsl_matrix* object and that is what m is.
// Also gsl_matrix_get returns a double, not an integer, so %i was changed to
// %lg.
  int status = printf( "half m_{00} = %lg\n", gsl_matrix_get(m,0,0)*half ); // gsl_matrix_get returns double

//Freeing the matrix should not be done here.
//gsl_matrix_free(m);

// Changed status to 0 upon sucessfull completeness for comparison tests
	return 0;
}

int main(void)
{
  // Changed (0, 0) to (1, 1) matrix_alloc takes n>0 only. This might have been
  // from a confusion with counting indices from 0. But dimensions are strictly
  // positive (non-zero)
  gsl_matrix* m = gsl_matrix_alloc(1,1);

  // Removed ampersand, as m is already pointer to matrix.
  gsl_matrix_set(m,0,0,66);

  printf("half m_{00} (should be 33):\n");

  // Removed asterisk as print returns integer and not pointer to integer. Also
  // ampersand was removed again.
  int status = print_half_00(m);

  // Removed asterisks, since status is an integer now. Also %g was changed to
  // %i since status is not a double, but an integer.
  if(status>0){
    printf("status=%i : SOMETHING WENT TERRIBLY WRONG (status>0)\n",status);
  }
  else{
    printf("status=%i : everything went just fine (status=0)\n",status);
  }

  // Removed ampersand
  gsl_matrix_free(m);
return 0;
}
