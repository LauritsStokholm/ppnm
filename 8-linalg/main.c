// Standard libraries
#include<stdio.h>
#include<stdlib.h>

// For linear algebra:
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>

#include"myheaders.h"

#include"task1_linalg.c"
#include"task2_eigenv.c"


int main(int argc, char* argv[])
{
  task1_linalg();
  task2_eigenv();
  return 0;
}
