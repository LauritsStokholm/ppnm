#include<stdio.h>
#include"komplex.h"
#include<math.h>

void komplex_print (char *s, komplex a) {
  printf ("%s (%g,%g)\n", s, a.re, a.im);
}

void komplex_set (komplex* z, double x, double y) {
  (*z).re = x;
  (*z).im = y;
}

komplex komplex_new (double x, double y) {
  komplex z = { x, y };
  return z;
}

komplex komplex_add (komplex a, komplex b) {
  komplex result = { a.re + b.re , a.im + b.im };
  return result;
}

komplex komplex_sub (komplex a, komplex b) {
  komplex result = { a.re - b.re, a.im - b.im };
  return result;
}

komplex komplex_sub_alternative (komplex a, komplex b) {
  komplex x = komplex_new(-b.re, -b.im);
  komplex result = komplex_add(a, x);
  return result;
}

komplex komplex_mul (komplex a, komplex b) {
  komplex result = {a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re};
  return result;
}

komplex komplex_conjugate (komplex z) {
  komplex x = komplex_new(z.re, -z.im);
  return x;
} 

double komplex_abs (komplex z) {
  komplex z_conj = komplex_conjugate(z);
  komplex abs_val = komplex_mul(z, z_conj);
  double result = sqrt(abs_val.re);
  return result;
}

komplex komplex_div (komplex a, komplex b) {
  komplex result = komplex_new((a.re*b.re+a.im*b.im)/(pow(komplex_abs(b), 2)), (a.im*b.re - a.re*b.im)/pow(komplex_abs(b), 2));
  return result;
}

komplex komplex_exp (komplex z) {
  komplex result = komplex_new(exp(z.re)*cos(z.im), exp(z.re)*sin(z.im));
  return result;
}

komplex komplex_cos (komplex z) {
  komplex result = komplex_new(cos(z.re)*cosh(z.im), -sin(z.re)*sinh(z.im));
  return result;
}

komplex komplex_sin (komplex z) {
  komplex result = komplex_new(sin(z.re)*cosh(z.im), cos(z.re)*sinh(z.im));
  return result;
}

komplex komplex_sqrt (komplex z) {
  double real = sqrt((komplex_abs(z) + z.re)/2);
  double imag = sqrt((komplex_abs(z) - z.re)/2);

  komplex result = komplex_new(real, imag);
  return result;
}



int real_equal(double a, double b, double tau, double epsilon){

  /* .......................................................
   * As part of exercise 3: Returns 1 if a and b are equal with
   * absolute precision tau or relative precision epsilon:
   *
   * |a-b|<tau
   * or 
   * |a-b|/|a+b| < epsilon/2
   * ....................................................... */
  double x = fabs(a-b);
  double y = fabs(a+b);
  double absolute_precision = x;
  double relative_precision = x/y;

  if(absolute_precision<tau)
    return 1;
  if(relative_precision<epsilon)
    return 1;
  else
    return 0;
}


int komplex_equal (komplex a, komplex b, double abs_eps, double rel_eps)
{
  if (real_equal(a.re, b.re, abs_eps, rel_eps) && real_equal(a.im, b.im, abs_eps, rel_eps))
  {
    return 1;
  }
  else return 0;
}

