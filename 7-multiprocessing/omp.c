#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>
#include "omp.h"                  // CFLAGS += -fopenmp -lgomp

void throw_dart(int n, int seedp, int* lp)
{
  double x, y;
  for (int i=0; i<n; i++)
  {
    // Generate point (x, y)
    x = (double) rand_r(&seedp) / RAND_MAX;
    y = (double) rand_r(&seedp) / RAND_MAX;
    // Is it within circle ?
    if (x*x + y*y < 1) {*lp += 1;}
  }
}


int
main(int argc, char** argv)
{

  int n = 1000, l0 = 0, l1 = 0, l2 = 0, seedp0 = 1, seedp1 = 13, seedp2 = 42;
  while (1)
  {
    int opt = getopt (argc, argv, "n:");
    if (opt == -1) break;
    switch (opt)
    {
      case 'n': n = atoi (optarg); break;
      default:
      fprintf (stderr, "Usage: %s [-n Iteration number]\n", argv[0]);
    exit (EXIT_FAILURE);
    }
  }

  #pragma omp parallel sections
  // the following sections wil be run parallelly in separate threads
   {
   #pragma omp section  // first thread will run this block of code
      {  
        throw_dart(n, seedp0, &l0);
      }  
   #pragma omp section  // second thread will run this block of code
      {  
        throw_dart(n, seedp1, &l1);
      }  
   #pragma omp section  // second thread will run this block of code
      {  
        throw_dart(n, seedp2, &l2);
      }  
   }
   int L = l0+l1+l2;
   int N = 3*n;
   double pi = 4 * (double) L / (double) N;
   double absolute_error = fabs(M_PI - pi);
   printf("%i\t%lg\t%lg\n", N, absolute_error, pi);
}
