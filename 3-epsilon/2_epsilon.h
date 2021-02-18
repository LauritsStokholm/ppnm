#ifndef HAVE_EPSILON_H
#define HAVE_EPSILON_H
#include <stdio.h>
#include <float.h>
#include <limits.h>
#endif

/* FLT_EPSILON: Float data types: */
void
epsilon (void)
{
  float a=1;
  while(1+a!=1){a/=2;}
  a*=2;

  float b=1;
  do {b/=2;}
  while(1+b!=1);
  b*=2;

  float c;
  for(c=1; 1+c!=1; c/=2){}
  c*=2;


 /* DBL_EPSILON: Double data types */
  double x=1;
  while(1+x!=1){x/=2;}
  x*=2;

  double y=1;
  do {y/=2;}
  while(1+y!=1);
  y*=2;

  double z;
  for(z=1; 1+z!=1; z/=2){}
  z*=2;

 /* LDBL_EPSILON: Long double data types */
  long double l=1;
  while(1+l!=1){l/=2;}
  l*=2;

  long double m=1;
  do {m/=2;}
  while(1+m!=1);
  m*=2;

  long double n;
  for(n=1; 1+z!=1; z/=2){}
  n*=2;


  printf("FLT_EPSILON calculated\n");
  printf("%.8f (while)\n", a);
  printf("%.8f (do-while)\n",b);
  printf("%.8f (for-loop)\n",c);
  printf("%.8f (floats.h)\n\n", FLT_EPSILON);

  printf("DBL_EPSILON calculated\n");
  printf("%.8lg (while)\n", x);
  printf("%.8lg (do-while)\n",y);
  printf("%.8lg (for-loop)\n",z);
  printf("%.8lg (floats.h)\n\n", DBL_EPSILON);

  printf("LDBL_EPSILON calculated\n");
  printf("%.8Lg (while-loop)\n",l);
  printf("%.8Lg (do-while-loop)\n",m);
  printf("%.8Lg (for-loop)\n",n);
  printf("%.8Lg (floats.h)\n\n", LDBL_EPSILON);
}
