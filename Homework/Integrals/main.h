#include <stdio.h>
#include <math.h>


/////////////////////////// TASK A ///////////////////////////////
/* ............................................................
 * This is the recursive part
 * ............................................................ */
void
recursive_part_with_reuse (
    double f(double), double a, double b,
    double abs, double rel, double f2, double f3,
    double* result, double* error);


/* ............................................................
 * f: integrand
 * [a, b]: integral interval
 * abs, rel: tolerance (absolute/relative)
 * result, error: pointers to result and the error
 * ............................................................ */
void
adaptive_integrator (
    double f(double), double a, double b,
    double abs, double rel, double* result, double* error);


/////////////////////////// TASK B ///////////////////////////////
void
clenshaw_curtis (
    double f(double), double a, double b, double abs, double rel,
    double* result, double* error);


void
recursive_clenshaw (
    double f(double), double a, double b, 
    double abs, double rel, double f2, double f3, double* result, double* error);


