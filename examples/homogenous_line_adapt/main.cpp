#include "hermes1d.h"

// ********************************************************************
// This example shows solution of harmonic steady state on the homogeneous 
// transmision line. Wave propagation is described by two linear 
// differential equation with complex coefficients.
// dU(x)/dx = -(R+j\omega L)I(x)
// dI(x)/dx = -(G+j\omega C)U(x)
// These equations are rewrited into four equations with real coefficients
// Ur' - R*Ir + \omega*L*Ii = 0
// Ui' - R*Ii - \omega*L*Ir = 0
// Ir' - G*Ur + \omega*C*Ui = 0
// Ii' - G*Ui - \omega*C*Ur = 0

// in an interval (0, 10) equipped with Dirichlet bdy conditions
// x1(0) = 1, x2(0) = 0, x3(0) = 0, x4(0) = 0

// General input:
int N_eq = 4;                           // Number of equations
int N_elem = 1000;                      // Number of elements
double A = 0, B = 10;                   // Domain end points
int P_init = 2;                         // Initial polynomial degree

// Equation parameters
double L = 25e-9;                       // induktance [H/m]
double C = 10e-12;                      // capacitance [F/m]
double G = 1e-9;                        // conductance [S/m]
double R = 1e-3;                        // resistance [Ohm/m]
double omega = 2*M_PI*3e8;      
double Zl = 60;                         // load impedance[Ohm]

// Error tolerances
double TOL_NEWTON_COARSE = 1e-4;        // Coarse mesh
double TOL_NEWTON_REF = 1e-4;           // Reference mesh

// Adaptivity
const int ADAPT_TYPE = 0;               // 0... hp-adaptivity
                                        // 1... h-adaptivity
                                        // 2... p-adaptivity
const double THRESHOLD = 0.7;           // Refined will be all elements whose error
                                        // is greater than THRESHOLD*max_elem_error
const double TOL_ERR_REL = 1e-2;        // Tolerance for the relative error between 
                                        // the coarse mesh and reference solutions
const int NORM = 1;                     // To measure errors:
                                        // 1... H1 norm
                                        // 0... L2 norm

// Boundary conditions
double Val_dir_left_1 = 1;       // real part of the voltage at the beginnig of the line
double Val_dir_left_2 = 0;       // imaginary part of the voltage at the beginnig of the line
double Val_dir_left_3 = 0;
double Val_dir_left_4 = 0;

// At the end of the line, there is an indirect boundary condition U(l) = I(l)*Zl (see below)

// Exact solution:
const int EXACT_SOL_PROVIDED = 0;
double exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
}

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left_1);
  mesh->set_bc_left_dirichlet(1, Val_dir_left_2);
  mesh->set_bc_left_dirichlet(2, Val_dir_left_3);
  mesh->set_bc_left_dirichlet(3, Val_dir_left_4);
  int N_dof = mesh->assign_dofs();
  printf("N_dof = %d\n", N_dof);

  // Create discrete problem on coarse mesh
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_1_1);
  dp->add_matrix_form(0, 2, jacobian_1_3);
  dp->add_matrix_form(0, 3, jacobian_1_4);
  dp->add_matrix_form(1, 1, jacobian_2_2);
  dp->add_matrix_form(1, 2, jacobian_2_3);
  dp->add_matrix_form(1, 3, jacobian_2_4);
  dp->add_matrix_form(2, 0, jacobian_3_1);
  dp->add_matrix_form(2, 1, jacobian_3_2);
  dp->add_matrix_form(2, 2, jacobian_3_3);
  dp->add_matrix_form(3, 0, jacobian_4_1);
  dp->add_matrix_form(3, 1, jacobian_4_2);
  dp->add_matrix_form(3, 3, jacobian_4_4);
  dp->add_vector_form(0, residual_1);
  dp->add_vector_form(1, residual_2);
  dp->add_vector_form(2, residual_3);
  dp->add_vector_form(3, residual_4);
  dp->add_matrix_form_surf(0, 0, jacobian_surf_right_U_Re, BOUNDARY_RIGHT);
  dp->add_matrix_form_surf(0, 2, jacobian_surf_right_U_Im, BOUNDARY_RIGHT);
  dp->add_matrix_form_surf(1, 1, jacobian_surf_right_I_Re, BOUNDARY_RIGHT);
  dp->add_matrix_form_surf(1, 3, jacobian_surf_right_I_Im, BOUNDARY_RIGHT);

  // Allocate vectors res and y_prev
  double *y_prev = new double[N_dof];
  if (y_prev == NULL) error("y_prev could not be allocated in main().");

  // Set y_prev zero
  for(int i=0; i<N_dof; i++) y_prev[i] = 0; 

  // Initial Newton's loop on coarse mesh
  int success, iter_num;
  success = newton(dp, mesh, y_prev, TOL_NEWTON_COARSE, iter_num);
  if (!success) error("Newton's method did not converge."); 
  printf("Finished initial coarse mesh Newton's iteration (%d iter).\n", 
         iter_num);

  // Replicate coarse mesh including dof arrays
  Mesh *mesh_ref = mesh->replicate();

  // Refine entire mesh_ref uniformly in 'h' and 'p'
  int start_elem_id = 0; 
  int num_to_ref = mesh_ref->get_n_active_elem();
  mesh_ref->reference_refinement(start_elem_id, num_to_ref);
  int N_dof_ref = mesh_ref->get_n_dof();
  printf("Fine mesh created (%d DOF).\n", N_dof_ref);

  // Allocate vector y_prev_ref
  double *y_prev_ref = new double[N_dof_ref];
  if (y_prev_ref == NULL) 
    error("y_prev_ref could not be allocated in main().");

  // Transfer coarse mesh solution to the fine mesh
  transfer_solution(mesh, mesh_ref, y_prev, y_prev_ref);
  printf("Coarse mesh solution copied to fine mesh.\n");

  // Convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Convergence History", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");

  // Main adaptivity loop
  int adapt_iterations = 1;
  while(1) {
    printf("============ Adaptivity step %d ============\n", adapt_iterations); 

    // Newton's loop on fine mesh
    success = newton(dp, mesh_ref, y_prev_ref, TOL_NEWTON_REF, iter_num);
    if (!success) error("Newton's method did not converge."); 
    printf("Finished fine mesh Newton's iteration (%d iter).\n", 
           iter_num);

    // Starting with second adaptivity step, obtain new coarse 
    // mesh solution via Newton's method. Initial condition is 
    // the last coarse mesh solution.
    if (adapt_iterations > 1) {

      // Newton's loop on coarse mesh
      success = newton(dp, mesh, y_prev, TOL_NEWTON_COARSE, iter_num);
      if (!success) error("Newton's method did not converge."); 
      printf("Finished coarse mesh Newton's iteration (%d iter).\n", 
             iter_num);
    }

    // In the next step, estimate element errors based on 
    // the difference between the fine mesh and coarse mesh solutions. 
    double err_est_array[MAX_ELEM_NUM]; 
    double err_est_total = calc_elem_est_errors(NORM, 
           mesh, mesh_ref, y_prev, y_prev_ref, err_est_array);

    // Calculate the norm of the fine mesh solution
    double ref_sol_norm = calc_approx_sol_norm(NORM, mesh_ref, y_prev_ref);

    // Calculate an estimate of the global relative error
    double err_est_rel = err_est_total/ref_sol_norm;
    printf("Relative error (est) = %g %%\n", 100.*err_est_rel);

    // If exact solution available, also calculate exact error
    if (EXACT_SOL_PROVIDED) {
      // Calculate element errors wrt. exact solution
      int order = 20; // heuristic parameter
      double err_exact_total = calc_exact_sol_error(NORM, 
              mesh, y_prev, exact_sol, order);
     
      // Calculate the norm of the exact solution
      // (using a fine subdivision and high-order quadrature)
      int subdivision = 500; // heuristic parameter
      double exact_sol_norm = calc_exact_sol_norm(NORM, 
                   exact_sol, N_eq, A, B, subdivision, order);
      // Calculate an estimate of the global relative error
      double err_exact_rel = err_exact_total/exact_sol_norm;
      printf("Relative error (exact) = %g %%\n", 100.*err_exact_rel);
      graph.add_values(0, N_dof, 100 * err_exact_rel);
    }

    // add entry to DOF convergence graph
    graph.add_values(1, N_dof, 100 * err_est_rel);

    // Decide whether the relative error is sufficiently small
    if(err_est_rel*100 < TOL_ERR_REL) break;

    // debug
    //if (adapt_iterations == 2) break;

    // Returns updated coarse and fine meshes, with the last 
    // coarse and fine mesh solutions on them, respectively. 
    // The coefficient vectors and numbers of degrees of freedom 
    // on both meshes are also updated. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, err_est_array,
          mesh, mesh_ref, y_prev, y_prev_ref, N_dof, N_dof_ref);

    adapt_iterations++;
  }

  // Plot meshes, results, and errors
  adapt_plotting(mesh, mesh_ref, y_prev, y_prev_ref,
           NORM, EXACT_SOL_PROVIDED, exact_sol);

  // Save convergence graph
  graph.save("conv_dof.gp");

  printf("Done.\n");
  delete [] y_prev;
  delete [] y_prev_ref;
  return 1;
}
