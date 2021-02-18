int
sum (void)
{
/* .......................................................
 * Using (x)n=INT_MAX/n for integer n's, calculate the series;
 * a) 1+1/2+1/3+1/4+...+1/(x)n
 * b) 1/(x)n + 1/((x)n - 1) + ... + 1/3 + 1/2 + 1
 *
 * Compare and explain the difference. Does the sum converge
 * as a function of n? [(x)n going towards INT_MAX.]
 *
 * Try different data types.
 *
 * ....................................................... */

 /* We will use the two integers i and j for summations */
  int i=1;  /* Sum up (0 is avoided) */
  int j=0;  /* Sum down (xn - xn) is avoided */

  int xn = INT_MAX/3;

  printf("IF ONE TRIES WITH xn = 10 or xn = 1000, one see the effect on \
      floating point precision series. They actually do diverge, depending\
      on summing up or down. Amazing!"); 

 /* Do-While-Loops */
  float sum_up_float_dowhile = 0;
  do {sum_up_float_dowhile+=1./i; i++;}
  while(i<=xn);

  float sum_down_float_dowhile = 0;
  do {sum_down_float_dowhile+=1./(xn-j); j++;}
  while(j<xn);

 /* Reset */
  i = 1;
  j = 0;

  double sum_up_double_dowhile = 0;
  do {sum_up_double_dowhile+=1./i; i++;}
  while(i<=xn);

  double sum_down_double_dowhile = 0;
  do {sum_down_double_dowhile+=1./(xn-j);
  j++;}
  while(j<xn);

 /* Reset */
  i = 1;
  j = 0;

  long double sum_up_longdouble_dowhile = 0;
  do {sum_up_longdouble_dowhile +=1./i; i++;}
  while(i<=xn);

  long double sum_down_longdouble_dowhile = 0;
  do {sum_down_longdouble_dowhile+=1./(xn-j);
  j++;}
  while(j<xn);

 /* Reset */
  i = 1;
  j = 0;


 /* For-Loops */
  float sum_up_float_for = 0;
  for(int k=1; k<=xn; k++){sum_up_float_for+=1./k;}

  float sum_down_float_for = 0;
  for(int k=xn; k>0; k--){sum_down_float_for+=1./k;}

  double sum_up_double_for = 0;
  for(int k = 1; k<=xn; k++){sum_up_double_for+=1./k;}

  double sum_down_double_for = 0;
  for(int k=xn; k>0; k--){sum_down_double_for+=1./k;}

  long double sum_up_longdouble_for = 0;
  for(int k=1; k<=xn; k++){sum_up_longdouble_for+=1./k;}

  long double sum_down_longdouble_for = 0;
  for(int k=xn; k>0; k--){sum_down_longdouble_for+=1./k;}


 /* While-Loops */
  float sum_up_float_while = 0;
  while(i<=xn){sum_up_float_while+=1./i; i++;}

  float sum_down_float_while = 0;
  while(j<xn){sum_down_float_while+=1./(xn-j); j++;}

 /* Reset */
  i = 1;
  j = 0;

  double sum_up_double_while = 0;
  while(i<=xn){sum_up_double_while+=1./i; i++;}

  double sum_down_double_while = 0;
  while(j<xn){sum_down_double_while+=1./(xn-j); j++;}

 /* Reset */
  i = 1;
  j = 0;

  long double sum_up_longdouble_while = 0;
  while(i<=xn){sum_up_longdouble_while+=1./i; i++;}

  long double sum_down_longdouble_while = 0;
  while(j<xn){sum_down_longdouble_while+=1./(xn-j); j++;}


  printf("The calculated series values for the data types float, double and \
longdouble are :\n");

  printf("Upward Summation\n");
  printf("Floating values:\n");
  printf("%f\t (dowhile) \n%f\t (while) \n%f\t (for) \n", sum_up_float_dowhile, sum_up_float_while, sum_up_float_for);
  printf("Double values:\n");
  printf("%lg\t (dowhile) \n%lg\t (while) \n%lg (for) \n", sum_up_double_dowhile, sum_up_double_while, sum_up_double_for);
  printf("Long double values:\n");
  printf("%Lg\t (dowhile) \n%Lg\t (while) \n%Lg (for) \n", sum_up_longdouble_dowhile, sum_up_longdouble_while, sum_up_longdouble_for);
  printf("\n");
  printf("Downward Summation\n");
  printf("Floating values:\n");
  printf("%f\t (dowhile) \n%f\t (while) \n%f\t (for) \n", sum_down_float_dowhile, sum_down_float_while, sum_down_float_for);
  printf("Double values:\n");
  printf("%lg\t (dowhile) \n%lg\t (while) \n%lg (for) \n", sum_down_double_dowhile, sum_down_double_while, sum_down_double_for);
  printf("Long double values:\n");
  printf("%Lg\t (dowhile) \n%Lg\t (while) \n%Lg (for) \n", sum_down_longdouble_dowhile, sum_down_longdouble_while, sum_down_longdouble_for);

  return 0;
}
