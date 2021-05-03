// Assume rescaled variables, so we use machine_epsilon for tolerance of
// step_size as well as jacobian step

#define delta sqrt(DBL_EPSILON)

// Derivative calculated numerically using finite differences
void
jacobian (void f (vector* x, void* params, vector* fx), vector* x, void* params, matrix* J)
{
  vector* fx     = vector_alloc(x->size);  vector_set_zero (fx);
  vector* fxplus = vector_alloc(x->size);  vector_set_zero (fxplus);
  vector*  xplus = vector_alloc(x->size);  vector_set_zero (xplus);
  vector* dx_k   = vector_alloc(x->size);  vector_set_zero (dx_k);

  f (x, params, fx);                 // calculate f(x)
  //fx = (f1(x), f2(x))

  for (int k=0; k<x->size; k++) // dxk
  {
    //set dx to (0, 0, ..., dx, 0, 0,...)
    vector_set_basis (dx_k, k);
    vector_scale     (delta, dx_k);

    // set xplus=x to xplus = x + dx_k
    vector_memcpy (x, xplus);  // Set xplus to x
    vector_add (xplus, dx_k);

    // Calculate (fi(x1, x2, ... xk+dx,..) - fi(x1, x2, .., xk,..) / dx)
    // (f1(xk+dx) - f1(xk), f2(xk+dx)-f2(x)....)
    f (xplus, params, fxplus);

    //fxplus = (f1(x+dx), f2(x+dx))
    vector_sub (fxplus, fx); // fxplus = fxplus - fx
    vector_scale (1./delta, fxplus);

    for (int i=0; i<x->size; i++) //dfi
    {
      matrix_set (J, i, k, vector_get(fxplus, i));
    }
  }
  vector_free (fx);
  vector_free (fxplus);
  vector_free (xplus);
  vector_free (dx_k);
}

void
newton (void f(vector* x, void* params, vector* fx), vector* x, void* params, double eps)
{
  /* ............................................................
   * function f: takes input vector x and fills fx = f(x)
   * vector x: on entry: the starting point; upon exit: aprox root
   * int n   : dimension of f = (f1, f2,.. fn);
   * double epsilon: accuracy goal: on exi the condition
   * norm(f(x)) < epsilon should be satisfied.
   * ............................................................ */

  int n = x->size;
  vector* fx = vector_alloc(n); //--- (..., fi(x1,..,xk+dx,..), ...)
  vector* fxplus = vector_alloc(n); //--- (..., fi(x1,..,xk+dx,..), ...)
  vector*  xplus = vector_alloc(n);
  vector* DX = vector_alloc (n);
  matrix* J  = matrix_alloc (n, n); //--- Jacobian
  matrix* Q  = matrix_alloc (n, n); //--- QR decomp
  matrix* R  = matrix_alloc (n, n); //--- QR decomp

  while (1)
  {
    // Calculate f(x)
    f (x, params, fx); // Set fx = (f1(x1, x2,..), f2(x1, x2,..),...)

    // Calculate df/dx
    jacobian (f, x, params, J); // Set jacobian; J_ik = dfi/dxk

    // Find direction, DX, by a type of gradient descent
    vector_scale (-1., fx);
    GS_decomp (J, Q, R);      // Ruins jacobian, but sets QR=J
    GS_solve  (Q, R, fx, DX);  // Solve the linear system J*DX = -f(x)
    vector_scale (-1., fx);    // Reset fx after scale in solve

    // Calculate xplus
    vector_memcpy (x, xplus); // stores x into xplus
    vector_add (xplus, DX);   // xplus = x + DX

    // Calculate f(x+DX)
    f (xplus, params, fxplus);

    // simple backtracking linesearch algorithm
    double step_size = 1;

    // If converged within tolerance or step_size is too small; return;
    double tolerance_x = 1./32;
    double tolerance_y = (1 - step_size/2.)*vector_norm(fx);
    double norm_fxplus = vector_norm (fxplus);

    while (step_size > tolerance_x && norm_fxplus > tolerance_y )
    {
      step_size /= 2.;
      vector_scale (step_size, DX);// DX = DX*step_size
      vector_add (xplus, DX);      // xplus = xplus + DX*step_size
      f (xplus, params, fxplus);
      tolerance_y = (1 - step_size/2.)*vector_norm(fx);
      norm_fxplus = vector_norm (fxplus);
    }

    // After minimised step
      //vector_scale (step_size, DX); (this is saved in the loop)
      vector_add (x, DX);

      if ( (vector_norm (DX) < delta) || (vector_norm(fxplus) < eps) )
      {
        printf("ended: final root is:\n");
        vector_printf (x);
        break;
      }
  }

  vector_free (fx); vector_free (fxplus);
  vector_free (xplus); vector_free (DX);
  matrix_free (Q); matrix_free (R);
  matrix_free (J);
  return;
}
