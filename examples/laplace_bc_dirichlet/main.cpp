#include "hermes1d.h"
#include "solver_umfpack.h"

// ********************************************************************

// This example solves the Poisson equation -u'' - f = 0 in
// an interval (A, B), equipped with Dirichlet boundary
// conditions on both end points. 

// General input:
static int N_eq = 1;
int N_elem = 3;                         // number of elements
double A = 0, B = 2*M_PI;               // domain end points
int P_init = 3;                         // initial polynomal degree

// Boundary conditions
double Val_dir_left = 1;                // Dirichlet condition left
double Val_dir_right = 1;               // Dirichlet condition right

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
  // create mesh
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left);
  mesh->set_bc_right_dirichlet(0, Val_dir_right);
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
  delete [] y_prev;
  return 1;
}
