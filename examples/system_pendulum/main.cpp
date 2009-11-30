#include "hermes1d.h"
#include "solver_umfpack.h"

// ********************************************************************
// This example solves the mathematical pendulum equation 
// y'' + k**2 * sin(y) = 0 in an interval (A, B), equipped with the 
// initial conditions y(A) = Init_angle, y'(0) = Init_vel. The 
// system is decomposed into two first order ODE and solved via 
// the Newton's method starting from zero initial condition.
// Note that the method diverges for longer time intervals, 
// depending on the interval length, number of elements, and 
// the initial polynomial degree.
//
// Derivation:
// m*l*u'' = -m*g*sin(u)
// so:
// u'' + k^2 * sin(u) = 0
// with k^2 = g/l
// so we have to solve a system of two nonlinear second-order equations
// v' + k^2 sin u = 0
// u' - v = 0
// in an interval (0, 2*pi) equipped with Dirichlet bdy conditions
// u(0) = 0, v(0) = k
// The approximate (linearized) solution is u(x) = sin(k*x), v(x) = k*cos(k*x)

// General input:
static int N_eq = 2;
int N_elem = 1292;          // number of elements
double A = 0, B = 10;     // domain end points
int P_init = 1;            // initial polynomal degree
double k = 0.5;

// Tolerance for the Newton's method
double TOL_NEWTON = 1e-5;

// Boundary conditions
double Init_angle = M_PI/2.;      // initial angle
double Init_vel = 0;              // initial velocity

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create mesh
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Init_angle);
  mesh->set_bc_left_dirichlet(1, Init_vel);
  int N_dof = mesh->assign_dofs();
  printf("N_dof = %d\n", N_dof);

  // Register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_0_0);
  dp->add_matrix_form(0, 1, jacobian_0_1);
  dp->add_matrix_form(1, 0, jacobian_1_0);
  dp->add_matrix_form(1, 1, jacobian_1_1);
  dp->add_vector_form(0, residual_0);
  dp->add_vector_form(1, residual_1);

  // Allocate vector y_prev
  double *y_prev = new double[N_dof];
  if (y_prev == NULL) error("res or y_prev could not be allocated in main().");

  // Set zero initial condition for the Newton's method
  for(int i=0; i<N_dof; i++) y_prev[i] = 0; 

  // Newton's loop
  int success, iter_num;
  success = newton(dp, mesh, y_prev, TOL_NEWTON, iter_num);
  if (!success) error("Newton's method did not converge."); 
  printf("Finished Newton's iteration (%d iter).\n", iter_num);

  // Plot the solution
  Linearizer l(mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  printf("Done.\n");
  delete [] y_prev;
  return 1;
}
