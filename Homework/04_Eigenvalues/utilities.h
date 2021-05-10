#ifndef HAVE_UTILITIES
#define HAVE_UTILITIES
// Jacobi diagonalisation with cyclic sweeps

// Calculate A*J
void
timesJ (matrix* A, int p, int q, double theta)
{
  double c = cos(theta);
  double s = sin(theta);

  for (int i=0; i<A->size1; i++)
  {
    double new_aip = c*matrix_get (A, i, p) - s*matrix_get (A, i, q);
    double new_aiq = s*matrix_get (A, i, p) + c*matrix_get (A, i, q);
    matrix_set (A, i, p, new_aip);
    matrix_set (A, i, q, new_aiq);
  }
}

// Calculate J*A
void
Jtimes (matrix* A, int p, int q, double theta)
{
  double c = cos(theta);
  double s = sin(theta);

  for (int j=0; j<A->size2; j++)
  {
    double new_apj =  c*matrix_get (A, p, j) + s*matrix_get (A, q, j);
    double new_aqj = -s*matrix_get (A, p, j) + c*matrix_get (A, q, j);
    matrix_set (A, p, j, new_apj);
    matrix_set (A, q, j, new_aqj);
  }
}

void
sweep (matrix* A, matrix* V)
{
  int n = A->size1; // A is symmetric (n=m)
  int changed = 1;
  while (changed)
  {
    changed=0;
    for (int p=0; p<n-1; p++)
    {
      for(int q=p+1;q<n;q++)
      {
        double apq = matrix_get (A, p, q);
        double app = matrix_get (A, p, p);
        double aqq = matrix_get (A, q, q);
        double theta = 0.5*atan2(2*apq, aqq-app);
        double c = cos(theta);
        double s = sin(theta);
        double new_app = c*c*app - 2*s*c*apq + s*s*aqq;
        double new_aqq = s*s*app + 2*s*c*apq + c*c*aqq;
        if (new_app!=app || new_aqq!=aqq) // do rotation
        {
          changed=1; 
          // NOTICE: J^T = J_inverse = J(-theta)
          timesJ (A, p, q,  theta); // A←A*J
          Jtimes (A, p, q, -theta); // A←J^T*A*J
          timesJ (V, p, q,  theta); // V←V*J
        }
      }
    }
  }
}


int
jacobi_diag (matrix* A, vector* eigenvals, matrix* eigenvecs)
{
  int n = A->size1; // A is symmetric (n=m)
  int changed = 1;
  int sweeps = 0;
  matrix_set_identity (eigenvecs); // For safety

  // If A already is diagonal
  for (int i=0; i<n; i++)
  {
    vector_set (eigenvals, i, matrix_get (A, i, i));
  }

  while (changed)
  {
    sweeps++;
    changed=0;
    for (int p=0; p<n-1; p++)
    {
      for (int q=p+1; q<n; q++)
      {
        double apq = matrix_get (A, p, q);
        double app = matrix_get (A, p, p);
        double aqq = matrix_get (A, q, q);

        double theta = 0.5*atan2(2*apq, aqq-app);
        double c = cos(theta);
        double s = sin(theta);

        double new_app = c*c*app - 2*s*c*apq + s*s*aqq;
        double new_aqq = s*s*app + 2*s*c*apq + c*c*aqq;

        // || is an logical or operation
        if (new_app!=app || new_aqq!=aqq) // do rotation
        {
          changed = 1;
          vector_set (eigenvals, p, new_app);
          vector_set (eigenvals, q, new_aqq);

//          for (int i=0; i<p; i++)
//          {
//            double aip = matrix_get (A, i, p);
//            double aiq = matrix_get (A, i, q);
//            matrix_set (A, i, p, c*aip - s*aiq);
//            matrix_set (A, i, p, s*aip + c*aiq);
//          }
//          for (int i=p+1; i<q; i++)
//          {
//            double api = matrix_get (A, p, i);
//            double aiq = matrix_get (A, i, q);
//            matrix_set (A, p, i, c*api - s*aiq);
//            matrix_set (A, i, p, s*api + c*aiq);
//          }
//          for (int i=q+1; i<n; i++)
//          {
//            double api = matrix_get (A, p, i);
//            double aqi = matrix_get (A, q, i);
//            matrix_set (A, p, i, c*api - s*aqi);
//            matrix_set (A, q, i, s*api + c*aqi);
//          }
//          for (int i=0; i<n; i++)
//          {
//            double vip = matrix_get (eigenvecs, i, p);
//            double viq = matrix_get (eigenvecs, i, q);
//            matrix_set (eigenvecs, i, p, c*vip - s*viq);
//            matrix_set (eigenvecs, i, q, s*vip + c*viq);
//          }
          // NOTICE: J^T = J_inverse = J(-theta)
          timesJ (A, p, q,  theta); // A←A*J
          Jtimes (A, p, q, -theta); // A←J^T*A*J
          timesJ (eigenvecs, p, q,  theta); // V←V*J
        }
      }
    }
  }
  return sweeps;
}

#endif
