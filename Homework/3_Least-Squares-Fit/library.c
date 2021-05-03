/* ............................................................
 * Within this libraryfile, essential constructable functions as 
 * well as simple algebraic functions for the customised structures 
 * vector and matrix (see main.h) are defined
 * ............................................................ */

/* Allocation */
//funlst* funlst_alloc (int n)
//{
//  funlst* F = malloc(sizeof(funlst));
//  (*F).size = n;
//  (*F).f = malloc(n*sizeof(double (*func)(double)));
//  return F;
//}

vector* vector_alloc (int n)
{
  vector* v = malloc(sizeof(vector));
  (*v).size = n;
  (*v).data = malloc(n*sizeof(double));
  return v;
}

matrix* matrix_alloc (int n, int m)
{
  matrix* A = malloc(sizeof(matrix));
  (*A).size1 = n; (*A).size2 = m;
  (*A).data = malloc(2*n*m*sizeof(double)); // WHY IS THIS NECESARRY WITH 2??
  return A;
}

/* Free memory */
//void funlst_free (funlst* F) { free((*F).f); free (F); }
void vector_free (vector* v) { free((*v).data); free(v); }
void matrix_free (matrix* A) { free((*A).data); free(A); }

/* Set data */
//void funlst_set (funlst* F, int i, double(*func)(double)) {(*F).f[i]=func;}
void vector_set (vector* v, int i, double x){(*v).data[i]=x;}
void matrix_set (matrix* A, int i, int j, double x)
{ (*A).data[i + j*(*A).size1] = x; }

/* Get data */
//double (*f)(double) funlst_get (funlst* F, int i){return (*F).f[i];}
double vector_get (vector* v, int i){return (*v).data[i];}
double matrix_get (matrix* A, int i, int j){return (*A).data[i + j*(*A).size1];}
// REMEMBER size1 instead of size2?


/* Generators */
void vector_set_random (vector* v)
{
  srand(time(NULL)); // Initialisation of RNG
  for (int i=0; i<v->size; i++)
  {
    double val = rand() / (double) RAND_MAX;// in [0, 1]
    vector_set (v, i, val);
  }
}

void vector_set_zero (vector* v)
{
  for (int i=0; i<v->size; i++) { vector_set(v, i, 0); }
}

void vector_set_basis (vector* v, int i)
{
  vector_set_zero(v);
  vector_set (v, i, 1);
}

void matrix_set_random (matrix* A)
{
  srand(time(NULL)); // Initialisation of RNG
  for (int i=0; i<A->size1; i++)
  {
    for (int j=0; j<A->size2; j++)
    {
      double val = rand() / (double) RAND_MAX;// in [0, 1]
      matrix_set (A, i, j, val);
    }
  }
}

void matrix_set_zero (matrix* A)
{
  for (int i=0; i<A->size1; i++)
  {
    for (int j=0; j<A->size2; j++)
    {
      matrix_set(A, i, j, 0);
    }
  }
}

void matrix_set_identity (matrix* A)
{
  for (int i=0; i<A->size1; i++)
  {
    for(int j=0; j<A->size2; j++)
    {
      if (i==j) matrix_set(A,i,j,1);
      else matrix_set(A,i,j,0);
    }
  }
}

/* Print */
void vector_fprintf (FILE* fp, vector* v)
{
  for (int i=0; i<v->size; i++)
  {
    fprintf(fp, "%.4lg\t", vector_get(v, i));
    if (i+1 == v->size){fprintf(fp, "\n");}
  }
}

void matrix_fprintf (FILE* fp, matrix* A)
{
  for (int i=0; i<A->size1; i++)
  {
    int j = 0;
    while (j < A->size2)
    {
      fprintf(fp, "%.4lg\t", matrix_get(A, i, j));
      j++;
      if (j == A->size2){fprintf(fp, "\n");}
    }
  }
}

/* Memory copy */
void matrix_memcpy (matrix* original, matrix* copy)
{
  assert((original->size1 == copy->size1) && (original->size2 == copy->size2));
  for (int i=0; i<original->size1; i++)
  {
    for (int j=0; j<original->size2; j++)
    {
      matrix_set(copy, i, j, matrix_get(original, i, j));
    }
  }
}

void vector_memcpy (vector* original, vector* copy)
{
  assert(original->size == copy->size);
  for (int i=0; i<original->size; i++)
  {
    vector_set (copy, i, vector_get(original, i));
  }
}


void matrix_column_view_vector (matrix* A, int j, vector* v)
{
  assert (v->size == A->size1);
  for (int i=0; i<A->size1; i++)
  {
    vector_set (v, i, matrix_get (A, i, j));
  }
}

void matrix_row_view_vector (matrix* A, int j, vector* v)
{
  assert (v->size == A->size2);
  for (int i=0; i<A->size2; i++)
  {
    vector_set (v, i, matrix_get (A, j, i));
  }
}

void vector_view_matrix_column (matrix* A, int j, vector* v)
{
  assert (v->size == A->size1);
  for (int i=0; i<v->size; i++)
  {
    matrix_set (A, i, j, vector_get(v, i));
  }
}

/* Algebraic functions */
/* Norm */
double matrix_column_norm (matrix* A, int j)
{
  double norm2 = 0;
  for (int i=0; i<A->size1; i++)
  {
    double aij = matrix_get(A, i, j);
    norm2 += aij*aij;
  }
  return sqrt(norm2);
}

/* Products */
// Calculates A*b and stores result in b;
void matrix_vector_product (matrix* A, vector* b, vector* Ab)
{
  assert (A->size2 == b->size);
  double sumk, aik, bk;
  for (int i=0; i<A->size1; i++)
  {
    sumk = 0;
    for (int k=0; k<A->size2; k++)
    {
      aik = matrix_get (A, i, k);
      bk = vector_get (b, k);
      sumk += aik*bk;
    }
    vector_set (Ab, i, sumk);
  }
}


// Calculates A*B and stores it in B;
void matrix_matrix_product (matrix* A, matrix* B, matrix* AB)
{
  assert (A->size2 == B->size1);

  double sumk, aik, bkj;
  // cij = sum_k aik*bkj
  for (int i=0; i< A->size1; i++)
  {
    for (int j=0; j<B->size2; j++)
    {

      sumk = 0;
      for (int k=0; k<A->size2; k++)
      {
        aik = matrix_get (A, i, k);
        bkj = matrix_get (B, k, j);
        sumk += aik*bkj;
      }
      matrix_set (AB, i, j, sumk);
    }
  }
}

void matrix_transpose_memcpy (matrix* src, matrix* dest)
{
  for (int i=0; i<src->size1; i++)
  {
    for (int j=0; j<src->size2; j++)
    {
      matrix_set(dest, j, i, matrix_get(src, i, j));
    }
  }
}


// Calculates the inverse of an upper triangular matrix
void matrix_inverse_upper_triangular (matrix* R, matrix* R_inv)
{
  assert (R->size1 == R->size2); // necesarry condition of invertibility

  // Used for ith standard basis element
  vector* ei = vector_alloc (R->size1);

  // Allocate inverse and set zero
  matrix_set_zero(R_inv);

  // Solve the system R*x = ei, which gives the i'th column of R_inverse
  for (int i=0; i<R->size1; i++)
  {
    // Set ei (0, 0 ,0 .... , 1, 0, 0.. .0)
    vector_set_basis (ei, i);
    backsub (R, ei);
    vector_view_matrix_column(R_inv, i, ei);
  }
}

void backsub (matrix* A, vector* b)
{
  for (int i =b->size-1; i>=0; i--)
  {
    double s = vector_get(b, i);
    for (int k=i+1; k<A->size1; k++)
    {
      s-= matrix_get(A, i, k)*vector_get(b, k);
    }
    vector_set (b, i, s/matrix_get (A, i, i));
  }
}

int vector_vector_test_equal (vector* v, vector* w, double epsabs, double epsrel)
{
  assert (v->size == w->size);
  double sum = 0;
  for (int i=0; i<v->size; i++)
  {
    sum += real_equal (vector_get (v, i), vector_get (w, i), epsabs, epsrel);
  }

  if (sum != v->size){return 1;}
  return 0;
}

void fprintf_result_vector_vector_test_equal (FILE* fp, vector* v, vector* w, double epsabs, double epsrel)
{
    if (vector_vector_test_equal (v, w, epsabs, epsrel) == 0)
    {
      fprintf(fp, "Result: Yes, they are equal.\n");
    }
    else {fprintf(fp, "Result, No, they are not equal.\n");}
}



// Test if A = B within relative and absolute tolerance
// return 0 if equal, and 1 if not equal.
void fprintf_result_matrix_matrix_test_equal (FILE* fp, matrix* A, matrix* B, double epsabs, double epsrel)
{
    if (matrix_matrix_test_equal(A, B, epsabs, epsrel) == 0)
    {
      fprintf(fp, "Result: Yes, they are equal.\n");
    }
    else {fprintf(fp, "Result, No, they are not equal.\n");}
}


int matrix_matrix_test_equal (matrix* A, matrix* B, double epsabs, double epsrel)
{
  assert (A->size1 == B->size1 && A->size2 == B->size2);
  double sum = 0;
  for (int i=0; i<A->size1; i++)
  {
    for (int j=0; j<A->size2; j++)
    {
      sum += real_equal (matrix_get (A, i, j), matrix_get (B, i, j), epsabs, epsrel);
    }
  }

  if (sum != A->size1*A->size2)
  {
    return 1;
  }

  return 0;
}

int real_equal (double a, double b, double tau, double epsilon){

  /* .......................................................
   * Returns 1 if a and b are equal with
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

