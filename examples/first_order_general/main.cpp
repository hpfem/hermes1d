#include "hermes1d.h"

// ********************************************************************

// This example solves the general first-order equation 
// y' = f(y, x) in an interval (A, B), equipped with the 
// initial condition y(A) = YA. The function f can be linear
// or nonlinear in 'y', as long as it is differentiable
// with respect to this variable (needed for the Newton's method). 

// General input:
static int N_eq = 1;                    // number of equations
int N_elem = 10;                        // number of elements
double A = 0, B = 10;                   // domain end points
double YA = 1;                          // equation parameter
int P_init = 2;                         // initial polynomal degree

// Tolerance for the Newton's method
double TOL_NEWTON = 1e-5;

// Function f(y, x)
double f(double y, double x) {
  return -y;
}

// Function dfdy(y, x)
double dfdy(double y, double x) {
  return -1;
}

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // create mesh
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, YA);
  int N_dof = mesh->assign_dofs();
  printf("N_dof = %d\n", N_dof);

  // register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian);
  dp->add_vector_form(0, residual);

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
  return 1;
}
