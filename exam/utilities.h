#ifndef HAVE_UTILITIES
#define HAVE_UTILITIES


// One-sided Jacobi diagonalisation for Singular Value Decomposition (SVD)
// with cyclic sweeps

// Calculate A*J
// J(theta, p, q) = diag(1, 1, ..., cos(theta), 1 .., sin(theta), 1...1)
// with off-diag = 0, i,j
//               = -sin(theta), if i,j = p,q
//               =  cos(theta), if i,j = q,p
void
timesJ (matrix* M, int p, int q, double theta)
{
  double c = cos(theta);
  double s = sin(theta);

  for (int i=0; i<M->size1; i++)
  {
    double new_mip = c*matrix_get (M, i, p) - s*matrix_get (M, i, q);
    double new_miq = s*matrix_get (M, i, p) + c*matrix_get (M, i, q);
    matrix_set (M, i, p, new_mip);
    matrix_set (M, i, q, new_miq);
  }
}


/* ............................................................
 * J is square, so for given A of size (n, m), then
 * A' = A*J -> J size (m, m)
 * ............................................................ */
int
SVD (matrix* A, matrix* V)
{
  int n = A->size1, m = A->size2, changed = 1, sweeps = 0;
  matrix_set_identity (V); // For safety

  // Column vectors
  vector* ap = vector_alloc (n);
  vector* aq = vector_alloc (n);

  while (changed)
  {
    sweeps++;
    changed=0;
    for (int p=0; p<n-1; p++)
    {
      for (int q=p+1; q<n; q++)
      {

        // The A matrix has (for p<q) a substructure on the form:
        //    app..apq
        //    .. .. ..
        //    aqp..aqq
        // for aqp=apq by symmetry
        matrix_column_view_vector (A, p, ap);
        matrix_column_view_vector (A, q, aq);

        double apq = vector_innerproduct (ap, aq);
        double aqq = vector_innerproduct (aq, aq);
        double app = vector_innerproduct (ap, ap);

        double theta = 0.5*atan2(2*apq, aqq-app);

        double c = cos(theta);
        double s = sin(theta);

        // For two-sided Jacobi
        //double new_app = c*c*app - 2*s*c*apq + s*s*aqq;
        //double new_aqq = s*s*app + 2*s*c*apq + c*c*aqq;

        //For one-sided Jacobi
        double new_app = c*app - s*apq;
        double new_aqq = s*apq + c*aqq;

        // || is an logical or operation
        if (new_app!=app || new_aqq!=aqq) // do rotation
        {
          changed = 1;
          //vector_set (D_diag, p, new_app);
          //vector_set (D_diag, q, new_aqq);

          timesJ (A, p, q, theta); // A←A*J
          timesJ (V, p, q, theta); // V←V*J
        }
      }
    }
  }
  vector_free (ap);
  vector_free (aq);
  return sweeps;
}

#endif
