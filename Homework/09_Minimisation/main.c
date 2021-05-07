#include "main.h"    // Header functions
#include "library.c" // Vector library
#include "qnewton.c"
#include "mysimplex.c"
#include "task_a.c"
#include "task_b.c"

int
main (int argc, char **argv)
{
  task_a();
  task_b();
  // task_c is to implement downhill simplex method (see mysimplex.c)
  return 0;
}

