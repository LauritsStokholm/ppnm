#include <stdlib.h>
#include <assert.h>

//typedef struct {int n; double *x, *y, *b, *c;} quadratic_spline;

quadratic_spline*
quadratic_spline_alloc (int n, double* x, double* y)
{
  quadratic_spline* s = (quadratic_spline*) malloc(sizeof(quadratic_spline));
  s->b = (double*) malloc((n-1)*sizeof(double));
  s->c = (double*) malloc((n-1)*sizeof(double));
  s->x = (double*) malloc(n*sizeof(double));
  s->y = (double*) malloc(n*sizeof(double));
  s->n = n;
  for (int i=0; i<n; i++)
  {
    s->x[i] = x[i];
    s->y[i] = y[i];
  }

  int i;
  double p[n-1], h[n-1];
  for (i=0; i<n-1; i++)
  {
    h[i]=x[i+1]-x[i];
    p[i]=(y[i+1]-y[i])/h[i];
  }

  s->c[0]=0;
  for (i=0; i<n-2; i++)
  {
    s->c[i+1]=(p[i+1]-p[i] - s->c[i]*h[i])/h[i+1];
    s->c[n-2]/=2;
  }

  for (i=n-3; i>=0; i--)
  {
    s->c[i] = (p[i+1]-p[i] - s->c[i+1]*h[i+1])/h[i];
  }

  for (i=0; i<n-1; i++)
  {
    s->b[i]=p[i] - s->c[i]*h[i];
  }

  return s;
}

double quadratic_spline_eval (int idx, quadratic_spline* s, double z)
{
//  assert( (z>= s->x[0]) && (z<= s->x[s->n-1]) );
//  int idx_min=0, idx_max=s->n-1;
//  while (idx_max - idx_min >1)
//  {
//    int mid = (idx_min+idx_max)/2.;
//    if (z> s->x[mid]) idx_min = mid; else idx_max=mid;
//  }
  double h = z - s->x[idx];
  return s->y[idx] + h*(s->b[idx] + h*s->c[idx]);
}

double quadratic_spline_deriv (int idx, quadratic_spline* s, double z)
{
  double h = z - s->x[idx];
  return s->b[idx] + 2*s->c[idx]*h;
}

double quadratic_spline_integ (int idx, quadratic_spline* s, double z)
{
  //double c0 =  s->y[i] + s->b[i]*s->x[i] + s->c[i]*s->x[i]*s->x[i];
  //double c1 = s->b[i]/2. - s->c[i];
  //double c2 = s->c[i]/3.;
  //return z*(c0 + c1*z + c2*z*z);

  double h = z - s->x[idx];
  double integ = s->y[idx]*h + s->b[idx]*h*h/2. + s->c[idx]*h*h*h/3.;

  int  k = 0;
  while (k<idx)
  {
    h = s->x[k+1] - s->x[k];
    integ += s->y[k]*h + s->b[k]*h*h/2. + s->c[k]*h*h*h/3.;
    k++;
  }

  return integ;
}

void quadratic_spline_free (quadratic_spline*s)
{
  free(s->x); free(s->y); free(s->b); free(s->c); free(s);
}


void quadratic_interp (int num_points, double* x, double* y, double resolution, FILE* fp_out)
{
  quadratic_spline* s = quadratic_spline_alloc (num_points, x, y);

  fprintf(fp_out, "x\teval\tderiv\tinteg\n");
  for (double w=x[0]; w<=x[num_points-1]; w+=resolution)
  {
  int i = binarysearch(num_points, x, w);

  double eval = quadratic_spline_eval (i, s, w);
  double deriv = quadratic_spline_deriv(i, s, w);
  double integ = quadratic_spline_integ(i, s, w);

  fprintf(fp_out, "%lg\t%lg\t%lg\t%lg\n", w, eval, deriv, integ);
  }

  quadratic_spline_free(s);
}
