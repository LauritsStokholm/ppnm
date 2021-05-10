#include<stdio.h>
#include<math.h>

#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>

#include"task1.c"
#include"task2.c"

int
main(void)
{
  task1();
  task2();

  return 0;
}
