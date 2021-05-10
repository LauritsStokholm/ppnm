#include "main.h"    // Header functions
#include "ode.h"     // ODE library
#include "library.h" // Vector library

void default_setting_of_y (vector* y)
{
  vector_set (y, 0 , -0.97000436) ; vector_set (y, 1,  0.24308753);/*(x0, y0)*/
  vector_set (y, 2 ,  0         ) ; vector_set (y, 3,  0         );/*(x1, y1)*/
  vector_set (y, 4 ,  0.97000436) ; vector_set (y, 5, -0.24308753);/*(x2, y2)*/
  vector_set (y, 6 ,  0.4662036850); vector_set (y, 7,  0.4323657300);/*(vx0, vy0)*/
  vector_set (y, 8 , -0.93240737) ; vector_set (y, 9, -0.86473146);/*(vx1, vy1)*/
  vector_set (y, 10,  0.4662036850) ; vector_set (y, 11, 0.4323657300);/*(vx2, vy2)*/
}



void
threebody_equation (double t, vector* y, vector* dydt)
{
  /* ............................................................
   * in XY-plane, y = ((x0, y0), (x1, y1), (x2, y2),
   *                   (vx0, vy0), (vx1, vy1), (vx2, vy2))
   *
   *                = (r0, r1, r2, v0, v1, v2)
   *                = (R, V)
   *
   * So dydt = (dr0dt, dr1dt, dr2dt, dv0dt, dv1dt, dv2dt)
   *         = (y[3],  y[4],  y[5],  ...  , .... ,  ... )
   * where the dots are filled according to the newtonian formula of
   * the threebody problem
   *
   * mi*dvi/dt = sum_j(not equal i) (G*mi*mj / |rj-ri|^2) * rj-ri / |rj-ri|
   * ............................................................ */

  double G = 1., m0 = 1., m1 = 1., m2 = 1.;

  // Positions of the planets
  double x0 = vector_get (y, 0); double y0 = vector_get (y, 1);
  double x1 = vector_get (y, 2); double y1 = vector_get (y, 3);
  double x2 = vector_get (y, 4); double y2 = vector_get (y, 5);

  // Velocities of the planets
  double vx0 = vector_get (y, 6);  double vy0 = vector_get (y, 7);
  double vx1 = vector_get (y, 8);  double vy1 = vector_get (y, 9);
  double vx2 = vector_get (y, 10); double vy2 = vector_get (y, 11);

  // Calculate rji = rj - ri (rij = -rji)
  double x01 = (x1 - x0); double y01 = (y1 - y0);
  double x02 = (x2 - x0); double y02 = (y2 - y0);
  double x12 = (x2 - x1); double y12 = (y2 - y1);

  //double r01 = x01 + y01;
  //double r02 = x02 + y02;
  //double r12 = x12 + y12;

  // Calculate the norm of rij (norm(rij) = norm(rji))
  double r01 = sqrt ( x01*x01 + y01*y01 );
  double r02 = sqrt ( x02*x02 + y02*y02 );
  double r12 = sqrt ( x12*x12 + y12*y12 );

  // Calculate dvxdt and dvydt for each

  // 1st planet
  double dvx0dt = G*( (m1*x01 / (r01*r01*r01) ) + (m2*x02 / (r02*r02*r02) ));
  double dvy0dt = G*( (m1*y01 / (r01*r01*r01) ) + (m2*y02 / (r02*r02*r02) ));

  // 2nd planet (Negative signs due to r01 = -r10)
  double dvx1dt = G*(-(m0*x01 / (r01*r01*r01) ) + (m2*x12 / (r12*r12*r12) ));
  double dvy1dt = G*(-(m0*y01 / (r01*r01*r01) ) + (m2*y12 / (r12*r12*r12) ));

  // 3rd planet
  double dvx2dt = G*(-(m0*x02 / (r02*r02*r02) ) - (m1*x12 / (r12*r12*r12) ));
  double dvy2dt = G*(-(m0*y02 / (r02*r02*r02) ) - (m1*y12 / (r12*r12*r12) ));

  // Setting dydt
  vector_set (dydt, 0,  vx0);    vector_set (dydt, 1,  vy0);
  vector_set (dydt, 2,  vx1);    vector_set (dydt, 3,  vy1);
  vector_set (dydt, 4,  vx2);    vector_set (dydt, 5,  vy2);
  vector_set (dydt, 6,  dvx0dt); vector_set (dydt, 7,  dvy0dt);
  vector_set (dydt, 8,  dvx1dt); vector_set (dydt, 9,  dvy1dt);
  vector_set (dydt, 10, dvx2dt); vector_set (dydt, 11, dvy2dt);
}

void
threebody (void)
{
  // System interval
  double a = 0, b = 4*M_PI;
  // h (step-size), abs and eps (absolute, relative precision)
  double h = 0.1, abs = 1e-3, eps = 1e-3;
  int n = 12, max = 10000000;
  int rkXY_method;
  vector* y = vector_alloc (n);


  rkXY_method = 12;
  default_setting_of_y (y);
  const char* string_name = "threebody.txt";
  my_ode_driver (threebody_equation, y, a, b, h, abs, eps, max, rkXY_method, string_name);

  rkXY_method = 23;
  default_setting_of_y (y);
  my_ode_driver (threebody_equation, y, a, b, h, abs, eps, max, rkXY_method, string_name);

  vector_free (y);
}

