#include "komplex.h"
#include <stdio.h>
#define TINY 1e-6

/* ............................................................
 * This exercises introduces structures; specifically we define
 * our own customisable "complex-number" structure komplex.
 * ............................................................ */
int main(void)
{
  // Initialize structures
  komplex a = {1,2}, b = {3,4};
  komplex z; komplex* ptr_z = &z;
printf("We test komplex_print...\n");
  komplex_print("a=",a);
  komplex_print("b=",b);

  printf("\nWe test komplex_set...\n");
  komplex_set(ptr_z, a.re, b.re);
  komplex x = {a.re, b.re};
  komplex* ptr_x = &x;
  komplex_print("z should=",x);
  komplex_print("z actually=",z);

  printf("\nWe test komplex_new...\n");
  z = komplex_new(a.im, b.im);
  komplex_set(ptr_x, a.im, b.im);
  komplex_print("z should=",x);
  komplex_print("z actually=",z);

  printf("\nWe test komplex_add...\n");
  z = komplex_add(a,b);
  komplex_set(ptr_x, 4, 6);
  komplex_print("a+b should =", x);
  komplex_print("a+b actually =", z);

  printf("\nWe test komplex_sub...\n");
  z = komplex_sub(a,b);
  komplex_set(ptr_x, -2, -2);
  komplex_print("a-b should =", x);
  komplex_print("a-b actually =", z);

  printf("\nWe test komplex_sub_alternative...\n");
  z = komplex_sub_alternative(a,b);
  komplex_set(ptr_x, -2, -2);
  komplex_print("a-b should =", x);
  komplex_print("a-b actually =", z);

  printf("\nWe test komplex_mul...\n");
  z = komplex_mul(a, b);
  komplex_set(ptr_x, -5, 10);
  komplex_print("a*b should =", x);
  komplex_print("a*b actually =", z);

  printf("\nWe test komplex_conjugate...\n");
  komplex za = komplex_conjugate(a);
  komplex zb = komplex_conjugate(b);
  komplex_set(ptr_x, 1, -2);
  komplex_print("a_conj should =", x);
  komplex_print("a_conj actually=", za);
  komplex_set(ptr_x, 3, -4);
  komplex_print("b_conj should =", x);
  komplex_print("a_conj actually=", zb);

  printf("\nWe test komplex_abs...\n");
  double abs_za = komplex_abs(a);
  double abs_zb = komplex_abs(b);
  double abs_x = 2.2360679775;
  printf("abs(a) should=%lg\n", abs_x);
  printf("abs(a) actually=%lg\n", abs_za);
  abs_x = 5;
  printf("abs(b) should=%lg\n", abs_x);
  printf("abs(b) actually=%lg\n", abs_zb);

  printf("\nWe test komplex_div...\n");
  za = komplex_div(a, b);
  zb = komplex_div(b, a);
  komplex_set(ptr_x, 0.44 ,0.08);
  komplex_print("a/b should =", x);
  komplex_print("a/b actually =", za);
  komplex_set(ptr_x, 2.2, -0.4);
  komplex_print("b/a should =", x);
  komplex_print("b/a actually =", zb);

  printf("\nWe test komplex_exp...\n");
  za = komplex_exp(a);
  zb = komplex_exp(b);
  komplex_set(ptr_x, -1.13120438 ,2.47172667);
  komplex_print("exp(a) should=", x);
  komplex_print("exp(a) actually=", za);
  komplex_set(ptr_x, -13.1287831 ,-15.2007845);
  komplex_print("exp(b) should=", x);
  komplex_print("exp(b) actually=", zb);

  printf("\nWe test komplex_sin...\n");
  za = komplex_sin(a);
  zb = komplex_sin(b);
  komplex_set(ptr_x, 3.16577851 ,1.95960104);
  komplex_print("sin(a) should=", x);
  komplex_print("sin(a) actually=", za);
  komplex_set(ptr_x, 3.85373804 ,-27.0168133);
  komplex_print("sin(b) should=", x);
  komplex_print("sin(b) actually=", zb);

  printf("\nWe test komplex_cos...\n");
  za = komplex_cos(a);
  zb = komplex_cos(b);
  komplex_set(ptr_x, 2.03272301 ,-3.0518978);
  komplex_print("cos(a) should=", x);
  komplex_print("cos(a) actually=", za);
  komplex_set(ptr_x, -27.0168133,-3.85373804 );
  komplex_print("cos(b) should=", x);
  komplex_print("cos(b) actually=", zb);

  printf("\nWe test komplex_sqrt...\n");
  za = komplex_sqrt(a);
  zb = komplex_sqrt(b);
  komplex_set(ptr_x, 1.27201965 , 0.786151378);
  komplex_print("sqrt(a) should=", x);
  komplex_print("sqrt(a) actually=", za);
  komplex_set(ptr_x, 2, 1);
  komplex_print("sqrt(b) should=", x);
  komplex_print("sqrt(b) actually=", zb);





/* the following is optional */

  printf("\nWe test komplex_equal...\n");
  komplex r = {1, 2}, dr = {TINY, TINY}, R = komplex_add(r, dr);

  if( komplex_equal(r,R,TINY,TINY) )
    printf("test 'komplex_equal' passed :) \n");
  else
    printf("test 'komplex_equal' failed: debug me, please... \n");
}
