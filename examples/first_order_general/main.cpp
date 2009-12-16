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
  // Create coarse mesh, set Dirichlet BC, enumerate 
  // basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, YA);
  printf("N_dof = %d\n", mesh->assign_dofs());

  // Register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian);
  dp->add_vector_form(0, residual);

  // Newton's loop
  int success, iter_num;
  success = newton(0, dp, mesh, TOL_NEWTON, iter_num);
  if (!success) error("Newton's method did not converge."); 
  printf("Finished Newton's iteration (%d iter).\n", iter_num);

  // Plot the solution
  Linearizer l(mesh);
  l.plot_solution("solution.gp");

  printf("Done.\n");
  return 1;
}
