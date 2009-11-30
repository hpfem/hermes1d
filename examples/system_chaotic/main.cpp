#include "hermes1d.h"
#include "solver_umfpack.h"

// ********************************************************************
// This example solves a nonlinear system of four first-order equations
// x1' - DAMPING*(1 - x2^2)*x1   +   x2 = 0
// x2'         -x1   +   x3 = 0
// x3'         -x2   +   x4 = 0
// x4'         -x3          = 0

// in an interval (0, 10) equipped with Dirichlet bdy conditions
// x1(0) = 1, x2(0) = 0, x3(0) = 0, x4(0) = 0

// General input:
static int N_eq = 4;
int N_elem = 500;           // number of elements
double A = 0, B = 20;       // domain end points
int P_init = 2;             // initial polynomal degree

// Damping parameter
int DAMPING_STEPS = 20;  // Number of damping steps. The entire problem
                         // will be run repeatedly, with the DAMPING parameter 
                         // increased from 0 to 1 in DAMPING_STEPS. Every time, 
                         // the last result is used as initial cond. for the 
                         // new computation.   
double DAMPING;          // DAMPING is an artificial param. used to 
                         // reduce the strength of the nonlinearity. 
                         // (The nonlinearity is multiplied with it.)

// Error tolerance
double TOL_NEWTON = 1e-5;// tolerance for the Newton's method 

// Boundary conditions
double Val_dir_left_1 = 1;
double Val_dir_left_2 = 0;
double Val_dir_left_3 = 0;
double Val_dir_left_4 = 0;

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // create mesh
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left_1);
  mesh->set_bc_left_dirichlet(1, Val_dir_left_2);
  mesh->set_bc_left_dirichlet(2, Val_dir_left_3);
  mesh->set_bc_left_dirichlet(3, Val_dir_left_4);
  int N_dof = mesh->assign_dofs();
  printf("N_dof = %d\n", N_dof);

  // register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_1_1);
  dp->add_matrix_form(0, 1, jacobian_1_2);
  dp->add_matrix_form(1, 0, jacobian_2_1);
  dp->add_matrix_form(1, 1, jacobian_2_2);
  dp->add_matrix_form(1, 2, jacobian_2_3);
  dp->add_matrix_form(2, 1, jacobian_3_2);
  dp->add_matrix_form(2, 2, jacobian_3_3);
  dp->add_matrix_form(2, 3, jacobian_3_4);
  dp->add_matrix_form(3, 2, jacobian_4_3);
  dp->add_matrix_form(3, 3, jacobian_4_4);
  dp->add_vector_form(0, residual_1);
  dp->add_vector_form(1, residual_2);
  dp->add_vector_form(2, residual_3);
  dp->add_vector_form(3, residual_4);

  // Allocate vector y_prev
  double *y_prev = new double[N_dof];
  if (y_prev == NULL) error("res or y_prev could not be allocated in main().");

  // Set y_prev zero
  for(int i=0; i<N_dof; i++) y_prev[i] = 0; 

  // Damping loop
  for(int damp_step = 1; damp_step < DAMPING_STEPS+1; damp_step++) {
    DAMPING = sin(damp_step*(1./DAMPING_STEPS)*M_PI/2.);

    printf("Damping: %g\n", DAMPING);

    // Newton's loop
    int success, iter_num;
    success = newton(dp, mesh, y_prev, TOL_NEWTON, iter_num);
    if (!success) error("Newton's method did not converge."); 
    printf("Finished Newton's iteration (%d iter).\n", iter_num);
  }

  // Plot the solution
  Linearizer l(mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  printf("Done.\n");
  delete[] y_prev;
  return 1;
}
