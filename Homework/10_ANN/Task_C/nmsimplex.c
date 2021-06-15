
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
 * REMEMBER: simplex is a matrix of size (dim, dim+1);
 * we conside its column vectors as edges of polytope and there are dim+1
 * edges!
 *
 * ............................................................ */

// Operations to be performed on the simplex (polytope)
void reflection
( vector* highest, vector* centroid, vector* reflected)
{
  assert (highest->size == centroid->size && centroid->size == reflected->size);
  for (int i=0; i<highest->size; i++)
  {
    double val = 2*vector_get (centroid, i) - vector_get (highest, i);
    vector_set (reflected, i, val);
  }
}

void expansion
( vector* highest, vector* centroid, vector* expanded)
{
  assert (highest->size == centroid->size && centroid->size == expanded->size);
  for (int i=0; i<highest->size; i++)
  {
    double val = 3*vector_get (centroid, i) - 2*vector_get (highest, i);
    vector_set(expanded, i, val);
  }
}

void contraction
( vector* highest, vector* centroid, vector* contracted)
{
  assert (highest->size == centroid->size && centroid->size == contracted->size);
  for (int i=0; i<highest->size; i++)
  {
    double val = 0.5*vector_get (centroid, i) + 0.5*vector_get (highest, i);
    vector_set (contracted, i, val);
  }
}

void reduction (matrix* simplex, int lo)
{
  for (int k=0; k<simplex->size2+1; k++) // n+1 edges
  {
    if (k != lo)
    {
      for (int i=0; i<simplex->size1; i++) // n dimensional space
      {
        double val = 0.5*(matrix_get (simplex, k, i) + matrix_get(simplex, lo, i));
        matrix_set (simplex, k, i, val);
      }
    }
  }
}


// Calculate norm of vector ab (distance between points)
double distance (vector* a, vector* b)
{
  assert (a->size == b->size);
  vector* c = vector_alloc (a->size);
  vector_memcpy (b, c); // c = b
  vector_sub (c, a);    // c = b - a

  double result = vector_norm (c);
  vector_free (c);
  return result;
}

// Calculates largest distance between any edge and the lowest point
double diameter (matrix* simplex, int lo)
{
  int d = simplex->size1;

  vector* simplex_lo = vector_alloc (d);
  vector* simplex_k  = vector_alloc (d);

  matrix_column_view_vector (simplex, lo, simplex_lo);

  double s = 0;
  // Iterate over all edges (vectors of simplex)
  for (int k=1; k<d+1; k++)
  {
    matrix_column_view_vector (simplex, k, simplex_k);
    double dist = distance (simplex_lo, simplex_k);
    if (dist>s) {s = dist;}
  }
  vector_free (simplex_lo);
  vector_free (simplex_k);
  return s;
}

void
simplex_update
(
  matrix* simplex,  // matrix of columnvectors = edges of polytope
  vector* f_values, // vector of evaluated values of objective function at edges
  int* hi,          // index of vector for highest function value edge
  int* lo,          // index of vector for lowest  function value edge
  vector* centroid  // center of gravity (point) of polytope
)
{
  int d = simplex->size1;

  // Before we do a linear search for indices
  *hi=0; *lo=0;
  double highest = vector_get (f_values, 0);
  double lowest  = vector_get (f_values, 0);

  // Set index values of vectors (0<index<simplex->size2)
  for (int k=1; k<d+1; k++) // dont compare the 0th with itself
  {
    double next = vector_get (f_values, k);
    if (next>highest) {highest=next; *hi=k;}
    if (next<lowest ) {lowest =next; *lo=k;}
  }


  for (int i=0; i<d; i++)// d dimensions
  {
    double s=0;
    for (int k=0; k<d+1; k++)// d+1 edges
    {
      if (k != *hi) // center of gravity (not counting the highest point)
      {
        s += matrix_get (simplex, i, k);
      }
    }
    vector_set (centroid, i, s/d);
  }
}


void
simplex_initiate
  (
    double fun (vector* x, void* params), //objective function
    void* params, // pointer to struct of parameter for objective function
    matrix* simplex, // matrix of polytope edges as column vectors
    vector* f_values, // vector of function evaluations at edges
    int* hi, int* lo, // pointers to indices (to be changed) of highest/lowest edge
    vector* centroid // center of gravity point (vector)
  )
{
  int d = simplex->size1; // simplex of d+1 edges in d dimensional space

  // k'th edge of polytope
  vector* simplex_k = vector_alloc (d);
  for (int k=0; k<d+1; k++)
  {
    matrix_column_view_vector (simplex, k, simplex_k);
    vector_set (f_values, k, fun(simplex_k, params) );
  }

  simplex_update (simplex, f_values, hi, lo, centroid);
  vector_free (simplex_k);
}



/* ............................................................
 * downhill simplex:
 * This is the function to be called outside this file
 * creates the simplex from starting point and step (taken in all directions)
 * on return, the start guess is changed to lowest edge
 *
 *  matrix* simplex, // matrix of polytope edges as column vectors
 * ............................................................ */
int
nmsimplex
(
  double fun (vector* x, void* params), // objective function
  vector* start,// a single edge
  vector* step, // step in each direction to create simplex of edges from start
  void* params, // pointer to struct of parameter for objective function
  double simplex_size_goal // epsilon
)
{

  // indices of highest, lowest and k'th edge
  int hi, lo, k = 0;

  int d = start->size; // d dimensional space and d+1 edges

  // create simplex from start point and step
  matrix* simplex = matrix_alloc (d, d+1);

  // set all edges to the start point
  for (int i=0; i<d+1; i++)
  {
    vector_view_matrix_column (simplex, i, start);
  }

  // Take step in each direction (diagonal)
  for (int i=0; i<d; i++)
  {
    matrix_set (simplex, i, i, vector_get(start, i) + vector_get (step, i) );
  }

  vector* centroid = vector_alloc (d);   // Center of gravity (point)
  vector* f_values = vector_alloc (d+1); // evaluation vals at d+1 edges
  vector* p1 = vector_alloc (d);
  vector* p2 = vector_alloc (d);

  vector* simplex_hi = vector_alloc (d);
  vector* simplex_lo = vector_alloc (d);

  // Calculate f_values, centroid, hi, lo
  simplex_initiate (fun, params, simplex, f_values, &hi, &lo, centroid);

  // As long as polytope has a diameter greater than our condition
  while (diameter (simplex, lo) > simplex_size_goal)
  {
    // Calculate hi, lo, and centroid
    simplex_update (simplex, f_values, &hi, &lo, centroid);

    // Set edges of highest and lowest points
    matrix_column_view_vector ( simplex, hi, simplex_hi);
    matrix_column_view_vector ( simplex, lo, simplex_lo);

    reflection (simplex_hi, centroid, p1);
    double f_re = fun (p1, params); // reflection of highest point
    double f_lo = vector_get (f_values, lo);
    double f_hi = vector_get (f_values, hi);

    if (f_re < f_lo) // reflection looks good: try expansion
    {
      expansion (simplex_hi, centroid, p2);
      double f_ex = fun(p2, params); // expandsion of highest point

      if (f_ex<f_re) // accept expansion
      {
        // set simplex_hi (higher point) to p2
        vector_memcpy (p2, simplex_hi); // update point
        vector_view_matrix_column (simplex, hi, simplex_hi); // update simplex
        vector_set (f_values, hi, f_ex); // update f_values vector
      }
      else {// reject expansion and accept reflection
        vector_memcpy (p1, simplex_hi); // update point
        vector_view_matrix_column (simplex, hi, simplex_hi); // update simplex
        vector_set (f_values, hi, f_re); // update f_values vector
      }
    }

    else {//reflection was not good
      if (f_re < f_hi)// ok, accept reflection
      {
        vector_memcpy (p1, simplex_hi); // update point
        vector_view_matrix_column (simplex, hi, simplex_hi); // update simplex
        vector_set (f_values, hi, f_re); // update f_values vector
      }
      else{ //try contraction
        contraction (simplex_hi, centroid, p1);
        double f_co = fun(p1, params); // contraction of highest point
        if (f_co < f_hi) // accept contration
        {
          vector_memcpy (p1, simplex_hi); // update point
          vector_view_matrix_column (simplex, hi, simplex_hi); // update simplex
          vector_set (f_values, hi, f_co); // update f_values vector
        }
        else{// do reduction
          reduction (simplex, lo);
          simplex_initiate (fun, params, simplex, f_values, &hi, &lo, centroid);
        }
      }
    }
    k++;
  }

  // Completed
  matrix_column_view_vector (simplex, lo, simplex_lo);
  vector_memcpy (simplex_lo, start); // set starting guess to lowest edge

  vector_free (centroid);
  vector_free (f_values);
  vector_free (p1);
  vector_free (p2);
  vector_free (simplex_hi);
  vector_free (simplex_lo);
  matrix_free (simplex);

  return k;
}
