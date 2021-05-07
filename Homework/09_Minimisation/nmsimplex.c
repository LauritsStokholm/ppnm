
/* ............................................................
 * Downhill simplex method / Nelder-Mead / amoeba
 *
 * advantages:    stability and lack of use of derivatives
 * disadvantages: relatively slow convergence compared to newton
 *
 * simplex: polytope/figure n+1 points (vertices) {p1, .. pn+1}
 *          (where pi is an n-dimensional vector)
 *
 * Highest point
 *
 * ............................................................ */

// Operations to be performed on the simplex (polytope)
void reflection
( double* highest, double* centroid, int dim, double* reflected)
{
  for (int i=0; i<dim; i++)
  {
    reflected[i] = 2*centroid[i] - highest[i];
  }
}

void expansion
( double* highest, double* centroid, int dim, double* expanded)
{
  for (int i=0; i<dim; i++)
  {
    expanded[i] = 3*centroid[i] - 2*highest[i];
  }
}

void contraction
( double* highest, double* centroid, int dim, double* contracted)
{
  for (int i=0; i<dim; i++)
  {
    contracted[i] = 0.5*centroid[i] + 0.5*highest[i];
  }
}

void reduction
( double* simplex[], int dim, int lo)
{
  for (int k=0; k<dim+1; k++)
  {
    if (k != lo)
    {
      for (int i=0; i<dim; i++)
      {
        simplex[k][i] = 0.5*(simplex[k][i] + simplex[lo][i]);
      }
    }
  }
}


// Calculate norm of vector ab (distance between points)
double distance (double* a, double* b, int dim)
{
  double s=0; 
  for (int i=0; i<dim; i++)
  {
    s+= pow(b[i] - a[i], 2);
  }
  return sqrt(s);
}

double size (double* simplex[], int dim)
{
  double s=0;
  for (int k=1; k<dim+1; k++)
  {
    double dist = distance(simplex[0], simplex[k], dim);
    if (dist>s) s=dist;
  }
  return s;
}

void
simplex_update
(
  double* simplex[],
  double* f_values,
  int d,
  int* hi,
  int* lo,
  double* centroid
)
{
  *hi=0; *lo=0;
  double highest = f_values[0], lowest = f_values[0];
  for (int k=1; k<d+1; k++)
  {
    double next = f_values[k];
    if (next>highest) {highest=next; *hi=k;}
    if (next<lowest ) {lowest =next; *lo=k;}
  }

  for (int i=0; i<d; i++)
  {
    double s=0;
    for (int k=0; k<d+1; k++)
    {
      if (k != *hi)
      {
        s += simplex[k][i];
      }
    }
    centroid[i] = s/d;
  }
}


void
simplex_initiate 
  (
    double fun (double*), double* simplex[],
    double* f_values, int d,
    int* hi, int* lo, double* centroid
  )
{
  for (int k=0; k<d+1; k++)
  {
    f_values[k] = fun(simplex[k]);
  }
  simplex_update (simplex, f_values, d, hi, lo, centroid);
}


int 
downhill_simplex
(
  double f(double*),       // objective
  double* simplex[],       // objective
  int d,                   // objective
  double simplex_size_goal // objective
)
{
  int hi, lo, k = 0;
  double centroid[d], f_value[d+1], p1[d], p2[d];
  simplex_initiate ( f, simplex, f_value, d, &hi, &lo, centroid);

  while (size(simplex, d)>simplex_size_goal)
  {
    simplex_update (simplex, f_value, d, &hi, &lo, centroid);
    reflection(simplex[hi], centroid, d, p1);
    double f_re = f(p1);

    if (f_re < f_value[lo]) // reflection looks good: try expansion
    {
      expansion (simplex[hi], centroid, d, p2);
      double f_ex = f(p2);
      if (f_ex<f_re) // accept expansion
      {
        for (int i=0; i<d; i++)
        {
          simplex[hi][i] = p2[i];
          f_value[hi] = f_ex;
        }
      }
      else {// reject expansion and accept reflection
        for (int i=0; i<d; i++)
        {
          simplex[hi][i] = p1[i];
          f_value[hi] = f_re;
        }
      }
    }
    else {//reflection was not good
      if (f_re < f_value[hi]){// ok, accept reflection
        for (int i=0; i<d; i++)
        {
          simplex[hi][i] = p1[i];
          f_value[hi] = f_re;
        }
      }
      else{ //try contraction
        contraction (simplex[hi], centroid, d, p1);
        double f_co = f(p1);
        if (f_co < f_value[hi]){ // accept contration
          for (int i=0; i<d; i++)
          {
            simplex[hi][i] = p1[i];
            f_value[hi] = f_co;
          }
        }
        else{// do reduction
          reduction (simplex, d, lo);
          simplex_initiate (f, simplex, f_value, d, &hi, &lo, centroid);
        }
      }
    }
    k++;
  }
  return k;
}
