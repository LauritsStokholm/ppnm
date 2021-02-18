int
intmax_dowhile (int arg)
{
  int i=0;
  if (arg == 0)
    do i++; while (i+1>i);
  else
    do i--; while (i>i-1);
  return i;
}

int
intmax_while (int arg)
{
  int i=0;
  if (arg == 0)
  {
    while (i+1>i) i++;
  }
  else
  {
    while (i-1<i) i--;
  }
  return i;
}

int
intmax_for (int arg)
{
  int i=0;
  if (arg == 0)
  {
  for (int j=0; j+1>j;j++)
  {
    i++;
  }
  }
  else
  {
    for (int j=0; j-1<j;j--)
    {
      i--;
    }
  }
  return i;
}

void
intmax (void)
{
  int I = intmax_while(0);
  int i = intmax_while(1);
  int J = intmax_dowhile(0);
  int j = intmax_dowhile(1);
  int K = intmax_for(0);
  int k = intmax_for(1);
  printf("Greatest integer, INT_MAX, by limits.h    %d\n", INT_MAX);
  printf("Greatest integer, I, by while-loop is:    %d\n", I);
  printf("Greatest integer, J, by do-while-loop is: %d\n", J);
  printf("Greatest integer, K, by for-loop is:      %d\n", K);
  printf("We can test the succesor function on these suprema\n");
  printf("by adding one to each the of integers\n");
  printf("INT_MAX+1 = %d\n", INT_MAX+1);
  printf("I+1       = %d\n", I+1);
  printf("J+1       = %d\n", J+1);
  printf("K+1       = %d\n", K+1);
  printf("\n\n\n");
  printf("Smallest integer, INT_MIN, by limits.h    %d\n", INT_MIN);
  printf("Smallest integer, i, by while-loop is:    %d\n", i);
  printf("Smallest integer, j, by do-while-loop is: %d\n", j);
  printf("Smallest integer, k, by for-loop is:      %d\n", k);
  printf("We can test the inverse succesor function on these suprema\n");
  printf("by subtracting one to each of the integers\n");
  printf("INT_MIN-1 = %d\n", INT_MIN-1);
  printf("i-1       = %d\n", i-1);
  printf("j-1       = %d\n", j-1);
  printf("k-1       = %d\n", k-1);
  printf("Realise, INT_MAX = (2^32)/2 - 1, and INT_MIN = (2^32)/2. This is as expected on a 32-bit signed binary integer.\n");
  printf("\n\n\n");
}
