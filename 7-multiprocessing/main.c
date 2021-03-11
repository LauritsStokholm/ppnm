#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>

typedef struct {int n; int* lp; unsigned int seedp;} my_params;


void* my_function (void* arg)
{
  my_params* params = (my_params*) arg;
  int  n = params -> n; // Number of points to scatter
  int* lp = params -> lp; // Pointer to number within area
  unsigned int seedp = params -> seedp; //

  int x, y;

  // Number of points withing circle
  for (int i=0; i<n; i++)
  {
    // Generate point (x, y)
    double x = (double) rand_r(&seedp) / RAND_MAX;
    double y = (double) rand_r(&seedp) / RAND_MAX;
    // Is it within circle ?
    if (x*x + y*y < 1) {*lp += 1;}
  }
  return NULL;
}

int
main (int argc, char **argv)
{

  // Unpacking parameters
  // Default
  int n = 1000, l0 = 0, l1 = 0, l2 = 0;
  unsigned int seedp0 = 1, seedp1 = 13, seedp2 = 42;

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


  my_params p0 = {.n = n, .lp = &l0, .seedp = seedp0}; // thread 0
  my_params p1 = {.n = n, .lp = &l1, .seedp = seedp1}; // thread 1
  my_params p2 = {.n = n, .lp = &l2, .seedp = seedp2}; // thread 1

  pthread_t t0, t1, t2; /* administration of threads */
  pthread_attr_t* attributes = NULL; /* default attributes */
  int flag0 = pthread_create(&t0, attributes, my_function, (void*)&p0);
  int flag1 = pthread_create(&t1, attributes, my_function, (void*)&p1);
  int flag2 = pthread_create(&t2, attributes, my_function, (void*)&p2);

  //my_function( (void*) &x2); /* meanwhile in the master thread... */

  void* retval = NULL; /* no return value from thread */
  flag0 = pthread_join(t0, retval); /* join the other thread */
  flag1 = pthread_join(t1, retval); /* join the other thread */
  flag2 = pthread_join(t2, retval); /* join the other thread */

  // Now total sum of each thread (for best estimation of pi)
  int L = l0 + l1 + l2;
  int N = 3*n;
  double pi = 4 * (double) L / (double) N;
  double absolute_error = fabs(M_PI - pi);

//  printf("Points within circle:%i\nTotal points:%i\n\n", 4*L, N);
//  printf("M_PI:%lg\nEstimate:%lg\n", M_PI, pi);
//  printf("To try new number of points type ./main -n [(int) your number]\n");
//  }
    printf("%i\t%lg\t%lg\n", N, absolute_error, pi);

  return 0;
}

