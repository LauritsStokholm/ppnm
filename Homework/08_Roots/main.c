#include "main.h"    // Header functions

#include "library.c"
#include "GS_decomp.c"
#include "newton.c"
#include "ode.c"
#include "task_a.c"
#include "task_bc.c"

int
main (int argc, char **argv)
{
  task_a();
  task_b();
  task_c();
  return 0;
}

