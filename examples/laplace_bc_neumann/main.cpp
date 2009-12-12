#include "hermes1d.h"
#include "solver_umfpack.h"

// ********************************************************************

// This example solves the Poisson equation -u'' - f = 0 in
// an interval (A, B), equipped with a Dirichlet boundary
// condition on the left and a Neumann BC on the right. 

// General input:
static int N_eq = 1;
int N_elem = 30;                       // number of elements
double A = 0, B = 2*M_PI;              // domain end points
int P_init = 1;                        // initial polynomal degree

// Boundary conditions
double Val_dir_left = 0;               // Dirichlet condition left
double Val_neum_right = 0;             // Neumann condition right
                                       // (derivative at end point)
// Tolerance for the Newton's method
double TOL_NEWTON = 1e-5;

// Function f(x)
double f(double x) {
  return sin(x);
  //return 1;
}

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate 
  // basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left);
  printf("N_dof = %d\n", mesh->assign_dofs());

  // Register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_vol);
  dp->add_vector_form(0, residual_vol);
  dp->add_vector_form_surf(0, residual_surf_right, BOUNDARY_RIGHT);

  // Newton's loop
  int success, iter_num;
  success = newton(dp, mesh, TOL_NEWTON, iter_num);
  if (!success) error("Newton's method did not converge."); 
  printf("Finished Newton's iteration (%d iter).\n", iter_num);

  // Plot the solution
  Linearizer l(mesh);
  l.plot_solution("solution.gp");

  printf("Done.\n");
  return 1;
}
