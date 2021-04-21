#ifndef HAVE_MYHEADER_H
#define HAVE_MYHEADER_H

#include <stdio.h>
#include <assert.h>
#include <getopt.h>

// gsl routines
#include <gsl/gsl_interp.h>

// My constructed functions 
// (one could just use cubic spline and set d=0 for the quadratic spline)
typedef struct {int n; double *x, *y, *b, *c;} quadratic_spline;
typedef struct {int n; double *x, *y, *b, *c, *d;} cubic_spline;

int    binarysearch  (int     num_points, double* x, double  z);
void   linear_interp (int     num_points, double* x, double* y, double resolution, FILE* fp_out);

void quadratic_interp (int num_points, double* x, double* y, double resolution, FILE* fp_out);
quadratic_spline* quadratic_spline_alloc (int n, double* x, double* y);
double quadratic_spline_eval (int idx, quadratic_spline* s, double z);
double quadratic_spline_deriv (int idx, quadratic_spline* s, double z);
double quadratic_spline_integ (int idx, quadratic_spline* s, double z);
void quadratic_spline_free (quadratic_spline*s);

void cubic_interp (int num_points, double* x, double* y, double resolution, FILE* fp_out);
cubic_spline* cubic_spline_alloc (int n, double *x, double* y);
double cubic_spline_eval (int idx, cubic_spline* s, double z);
double cubic_spline_deriv (int idx, cubic_spline* s, double z);
double cubic_spline_integ (int idx, cubic_spline* s, double z);
void cubic_spline_free (cubic_spline* s);

void gslInterp (int n, double* x, double* y, double resolution, FILE* output, const gsl_interp_type* method);

#endif
