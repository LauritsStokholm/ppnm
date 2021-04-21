#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

//typedef struct {int n; double *x, *y, *b, *c, *d;} cubic_spline;

cubic_spline* cubic_spline_alloc (int n, double *x, double* y)
{
  cubic_spline* s = (cubic_spline*) malloc(sizeof(cubic_spline));
  s->x = (double*) malloc(n*sizeof(double));
  s->y = (double*) malloc(n*sizeof(double));
  s->b = (double*) malloc(n*sizeof(double));
  s->c = (double*) malloc((n-1)*sizeof(double));
  s->d = (double*) malloc((n-1)*sizeof(double));
  s->n = n;

  for (int i=0; i<n; i++)
  {
    s->x[i] = x[i];
    s->y[i] = y[i];
  }

  double h[n-1], p[n-1];

  for (int i=0; i<n-1; i++)
  {
    h[i] = x[i+1]-x[i];
    assert(h[i]>0);
  }

  for (int i=0; i<n-1; i++)
  {
    p[i] = (y[i+1]-y[i])/h[i];
  }

  double D[n], Q[n-1], B[n]; // building the tridiagonal system

  D[0] = 2; D[n-1]=2;
  for (int i=0; i<n-2; i++)
  {
    D[i+1] = 2*h[i]/h[i+1] + 2;
  }

  Q[0] = 1;
  for (int i=0; i<n-2; i++)
  {
    Q[i+1] = h[i]/h[i+1];
  }

  B[0] = 3*p[0]; B[n-1] = 3*p[n-2]; 
  for (int i=0; i<n-2; i++)
  {
    B[i+1] = 3*(p[i]+p[i+1]*h[i]/h[i+1]);
  }

  // Gauss Elemination
  for (int i=1; i<n; i++)
  {
    D[i] -= Q[i-1]/D[i-1];
    B[i] -= B[i-1]/D[i-1];
  }

  s->b[n-1] = B[n-1]/D[n-1]; // back substitution

  for (int i=n-2; i>=0; i--)
  {
    s->b[i] = (B[i]-Q[i]*s->b[i+1])/D[i];
  }

  for(int i=0; i<n-1; i++)
  {
    s->c[i]=(-2*s->b[i] - s->b[i+1] + 3*p[i]) / h[i];
    s->d[i]=(   s->b[i] + s->b[i+1] - 2*p[i]) / h[i]/h[i];
  }

  return s;
}

void cubic_spline_free (cubic_spline* s)
{
  free(s->x); free(s->y); free(s->b); free(s->c); free(s->d); free(s);
}

double cubic_spline_eval (int idx, cubic_spline* s, double z)
{
  double h = z - s->x[idx];
  return s->y[idx] + h*(s->b[idx] + h*(s->c[idx] + h*s->d[idx]));
}

//////////
double cubic_spline_deriv (int idx, cubic_spline* s, double z)
{
  double h = z - s->x[idx];
  double result = s->b[idx] + h*(2*s->c[idx] + 3*h*s->d[idx]);
  return result;
}

//////////
double cubic_spline_integ (int idx, cubic_spline* s, double z)
{
  double h = z - s->x[idx];
  double integ = h*(s->y[idx] + s->b[idx]*h/2. + s->c[idx]*h*h/3. + s->d[idx]*h*h*h/4.);

  int  k = 0;
  while (k<idx)
  {
    h = s->x[k+1] - s->x[k];
    integ += h*(s->y[k] + s->b[k]*h/2. + s->c[k]*h*h/3. + s->d[k]*h*h*h/4.);
    k++;
  }

  return integ;
}


void cubic_interp (int num_points, double* x, double* y, double resolution, FILE* fp_out)
{
  cubic_spline* s = cubic_spline_alloc (num_points, x, y);

  fprintf(fp_out, "x\teval\tderiv\tinteg\n");
  for (double w=x[0]; w<=x[num_points-1]; w+=resolution)
  {
  int i = binarysearch(num_points, x, w);

  double eval  = cubic_spline_eval (i, s, w);
  double deriv = cubic_spline_deriv(i, s, w);
  double integ = cubic_spline_integ(i, s, w);

  fprintf(fp_out, "%lg\t%lg\t%lg\t%lg\n", w, eval, deriv, integ);
  }

  cubic_spline_free(s);
}
