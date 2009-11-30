#include "hermes1d.h"
#include "solver_umfpack.h"

// ********************************************************************
// This example solves a system of two linear second-order equations 
// - u'' + v - f_0 = 0
// - v'' + u - f_1 = 0
// in an interval (A, B) equipped with Dirichlet bdy conditions 
// u(A) = exp(A), u(B) = exp(B), v(A) = exp(-A), v(B) = exp(-B). 
// The exact solution is u(x) = exp(x), v(x) = exp(-x). 

// General input:
static int N_eq = 2;
int N_elem = 2;          // number of elements
double A = 0, B = 1;     // domain end points
int P_init = 2;          // initial polynomal degree

// Tolerance for the Newton's method
double TOL_NEWTON = 1e-5;

// Boundary conditions
double Val_dir_left_0 = exp(A);
double Val_dir_right_0 = exp(B);
double Val_dir_left_1 = exp(-A);
double Val_dir_right_1 = exp(-B);

// Function f_0(x)
double f_0(double x) {
  return -exp(x) + exp(-x);
}

// Function f_1(x)
double f_1(double x) {
  return -exp(-x) + exp(x);
}

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create mesh
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left_0);
  mesh->set_bc_right_dirichlet(0, Val_dir_right_0);
  mesh->set_bc_left_dirichlet(1, Val_dir_left_1);
  mesh->set_bc_right_dirichlet(1, Val_dir_right_1);
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
