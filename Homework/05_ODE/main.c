#include "main.h"
#include "library.h"
#include "ode.c"

#include "task_a.c"
#include "threebody.c"

int
main (int argc, char* argv[])
{
  task_a();
  threebody();
  return 0;
}
