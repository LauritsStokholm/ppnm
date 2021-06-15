// Assume rescaled variables, so we use machine_epsilon for tolerance of
// step_size as well as jacobian step

#define delta sqrt(DBL_EPSILON)
double alpha      = 1e-4; // Parameter for Armijo condition
double lambda_min = 1e-4; // If minimal lambda reached; take step unconditionally


/* ............................................................
 * Derivative calculated numerically using finite differences
 * ............................................................ */
void
numeric_gradient
  (
    double f (vector* x, void* params), // Objective function
    vector* x,                          // Variables (x, y, z...)
    void* params,                       // In case of parameters/coefficients
    vector* gradient                    // On return, calculated gradient
  )
{
 /* Preparation of vectors */
  double z0, z1, diffquotient;
  vector* xplus  = vector_alloc (x->size);  vector_set_zero (xplus);
  vector* dx_k   = vector_alloc (x->size);  vector_set_zero (dx_k);

  z0 = f (x, params); // z0 = f(x)
  for (int k=0; k<x->size; k++) // dxk
  {
    //set dx to (0, 0, ..., dx, 0, 0,...)
    vector_set_basis (dx_k, k);
    vector_scale     (delta, dx_k);

    // set xplus
    vector_memcpy (x, xplus);  // xplus = x
    vector_add (xplus, dx_k);  // xplus = x+dx

    // Calculate df/dxi = f(x1, x2..., xi+dx,..) - f(x1, fx2,... xi,..) / dx
    z1 = f (xplus, params);

    diffquotient = (z1 - z0) / delta;
    vector_set (gradient, k, diffquotient);
  }
  vector_free (xplus);
  vector_free (dx_k);
}



/* ............................................................
  * Avoid recalculations of the Hessian matrix at each iteration.
  * Instead: update based on gradient vector analysis
  * Let:
  * s = dx*lambda (newton step scaled by 0<lambda<1 optimised by linear
  *                backtracking)
  * B = inverse_Hessian (zeroth order approx = identity, and then updated
  *                      instead of calculated)
  * y = grad_f(x+s) - grad_f(x)
  * u = s - B*y
  * then dB (update of B for ith iteration) satisfy 
  * dB*y = u, which can be solved
  * ............................................................ */
int
qnewton
  (
    double  f(vector* x, void* params), // Objective function
    vector* x,        // Input: starting point; Return: approximation to root
    void*   params,   // Possible parameters / coefficients for f
    double  eps       // Accuracy goal, on exit |gradient| < eps
  )
{
  int n = x->size, nsteps = 0;

  // Prepare vectors
  vector* DX        = vector_alloc (n); // newtons step   s = lambda*DX
  vector* xplus     = vector_alloc (n); // next iteration xplus = x+s
  vector* gradient0 = vector_alloc (n); // gradF(x)
  vector* gradient1 = vector_alloc (n); // gradF(xplus)
  vector* y         = vector_alloc (n); // y = gradient1 - gradient0
  vector* u         = vector_alloc (n); // u = s - By
  vector* By        = vector_alloc (n); // By = B*y

  // Prepare matrices
  matrix* B  = matrix_alloc (n, n); // Approximation of inverse Hessian matrix
  matrix* dB = matrix_alloc (n, n); // Update of B

  // Prepare scalars
  double z0, z1, lambda, condition;

  // Set vectors and matrices for safety
  vector_set_zero (gradient0); vector_set_zero (gradient1);
  vector_set_zero (xplus); vector_set_zero (DX);
  matrix_set_identity (B); matrix_set_zero (dB);

  // Let the algorithm begin
  while (nsteps < 1000)
  {
    nsteps++;
    // Calculate DX [newtons step]
    numeric_gradient (f, x, params, gradient0); // Calculate grad_f(x)

    //printf ("Gradient\n");
    //vector_printf (gradient0);

    matrix_scale (B, -1); // Set B for the calculation dx = -B * grad_f(x)
    matrix_vector_product (B, gradient0, DX); // Calculate Newtons step in DX
    matrix_scale (B, -1); // Set B correctly back


    //printf ("DX:\n");
    //vector_printf (DX);
    // Modify Newtons step by scale 0<lambda<1 by backtracking
    lambda = 1.; // (initially) so vector_scale (lambda, DX) unnecessary

    z0 = f (x, params);      // f(x)
    //printf ("z0 = %lg\n", z0);

    vector_memcpy (x, xplus); // xplus = x
    vector_add (xplus, DX);   // xplus += lambda*DX;

    //printf ("xplus=\n"); vector_printf (xplus);
    z1 = f (xplus, params);   // f(x + s)

    //printf ("z1 = %lg\n", z1);

    double DXgradient0;
    vector_dot_product (DX, gradient0, &DXgradient0);
    //printf ("DXgradient0=%lg\n" ,DXgradient0);

    // Armijo condition
    condition = z0 + alpha*DXgradient0;
    printf ("z1: %lg\n", z1);
    printf ("cond: %lg\n", condition);
    while (z1 >= condition)
    {
      // Do scale Newtons step
      lambda /= 2.;
      vector_scale (1./2., DX);

      // Update xplus
      vector_memcpy (x, xplus); // xplus = x
      vector_add (xplus, DX);   // xplus += lambda*DX;

      // Update conditions
      z1 = f (xplus, params);
      //printf ("z1=%lg\n", z1);
      vector_dot_product (DX, gradient0, &DXgradient0);
      condition = z0 + alpha*DXgradient0;

      // If Armijo condition can not be fulfilled
      // (Reset DX and take full step)
      if (lambda < lambda_min)
      {
        matrix_set_identity (B);
        break;
      }
      // when departing, optimised n: s = (1/2)^n*DX been found
      // where lambda = (1/2)^n
    }

    // We now have xplus = x + s = x + lambda*DX for optimised s
    // Now calculate the update dB for the inverse Hessian matrix B
    // Symmetric rank-1 update (SR1) dB = v*vT for
    // dB = (uuT) / (uT*y)
    numeric_gradient (f, xplus, params, gradient1); // Calculate grad_f(x)
    vector_memcpy (gradient1, y);      // set y = gradient1
    vector_sub    (y, gradient0);      // set y = gradient1 - gradient0
    matrix_vector_product (B, y, By);  // Calculate By

    // Calculate u = s - By
    vector_set_zero (u);     // set u = 0
    vector_memcpy   (DX, u); // set u = s
    vector_sub      (u, By); // set u = s - By

    // Calculate uTy
    double uTy;
    vector_dot_product (u, y, &uTy);

    // Calculate vvT
    vector_outerproduct (u, u, dB);
    matrix_scale (dB, 1./uTy);

    // Only perform update if denominator is not too small
    if ( fabs (uTy) > 1e-6 ) // eps = 1e-6 is taken from the notes
    {
      matrix_sum (B, dB); // Update B = B+dB
    }

    // Update point x (take step)
    vector_memcpy (xplus, x);
    z0 = z1; // f(x) = f(x+s)

    // Conditions for breaking
    if (vector_norm (gradient1) < eps)
    {
      fprintf(stderr, "qnewton: |grad|<eps\n");
      break;
    }

    if ( vector_norm (DX) < delta*vector_norm(x) )
    {
      fprintf (stderr, "qnewton: |s|<delta*|x|\n");
      break;
    }

  } // end while-loop

  // Free memory
  vector_free (DX); vector_free (xplus); vector_free (gradient0);
  vector_free (gradient1); vector_free (y); vector_free (u);
  vector_free (By); matrix_free (B); matrix_free (dB);

  return nsteps;
}
