#include "main.h"
#include "library.h"
#include "utilities.h"
#include "task_a.c"
#include "task_b.c"

int
main (int argc, char* argv[])
{
  task_a();
  task_b();

  // Task C: We do timing seperately (see my_timer.c and gsl_timer.c)
  return 0;
}
