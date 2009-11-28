#include "hermes1d.h"

// ********************************************************************

// This test makes sure that an exact function 
// sin() is approximated with the relative 
// error 1e-5 with no more than 20 DOF.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

// General input:
static int N_eq = 1;
int N_elem = 2;                         // Number of elements
double A = -M_PI, B = M_PI;             // Domain end points
int P_init = 1;                         // Initial polynomal degree

// Stopping criteria for Newton
double TOL_NEWTON_COARSE = 1e-5;        // Coarse mesh
double TOL_NEWTON_REF = 1e-3;           // Reference mesh

// Adaptivity
const int ADAPT_TYPE = 0;               // 0... hp-adaptivity
                                        // 1... h-adaptivity
                                        // 2... p-adaptivity
const double THRESHOLD = 0.7;           // Refined will be all elements whose error
                                        // is greater than THRESHOLD*max_elem_error
const double TOL_ERR_REL = 1e-5;        // Tolerance for the relative error between 
                                        // the coarse mesh and reference solutions
const int NORM = 0;                     // To measure errors:
                                        // 1... H1 norm
                                        // 0... L2 norm

// Boundary conditions
double Val_dir_left = 0;                // Dirichlet condition left
double Val_dir_right = 0;               // Dirichlet condition right

// Function f(x)
double f(double x) {
  return sin(x);
  //return 2;
}

// Exact solution:
// When changing exact solution, do not 
// forget to update interval accordingly
const int EXACT_SOL_PROVIDED = 1;
double exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
  u[0] = sin(x);
  dudx[0] = cos(x);
  //u[0] = 1. - x*x;
  //dudx[0] = -2.*x;
}

// ********************************************************************

// bilinear form for the Jacobi matrix 
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...basis function
// v...test function
// u_prev...previous solution (all solution components)
double jacobian(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM], double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i]*dvdx[i]*weights[i];
  }
  return val;
};

// (nonlinear) form for the residual vector
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...approximate solution
// v...test function
// u_prev...previous solution (all solution components)
double residual(int num, double *x, double *weights, 
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM], double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (du_prevdx[0][i]*dvdx[i] - f(x[i])*v[i])*weights[i];
  }
  return val;
};

/******************************************************************************/
int main() {
  Mesh *mesh = NULL;                // coarse mesh 
  Mesh *mesh_ref = NULL;            // reference mesh
  Matrix *mat = NULL;               // Jacobi matrix (coarse mesh)
  Matrix *mat_ref = NULL;           // Jacobi matrix (reference mesh)
  double *y_prev = NULL;            // vector of unknown coefficients (coarse mesh)
  double *y_prev_ref = NULL;        // vector of unknown coefficients (reference mesh)
  double *res = NULL;               // residual vector (coarse mesh)
  double *res_ref = NULL;           // residual vector (reference mesh)
  DiscreteProblem *dp = NULL;       // discrete problem

  // convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Convergence History", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");


  // Create coarse mesh, set Dirichlet BC, enumerate basis functions
  mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left);
  mesh->set_bc_right_dirichlet(0, Val_dir_right);
  int N_dof = mesh->assign_dofs();
  printf("N_dof = %d\n", N_dof);

  // Create discrete problem on coarse mesh
  dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian);
  dp->add_vector_form(0, residual);

  // Main adaptivity loop
  int adapt_iterations = 1;
  while(1) {
 
    printf("============ Adaptivity step %d ============\n", adapt_iterations); 
    printf("------ Newton iteration on coarse mesh ----\n"); 
    printf("N_dof = %d\n", N_dof);

    // (Re)allocate Jacobi matrix and vectors y_prev and res
    if (mat != NULL) delete mat;
    mat = new CooMatrix(N_dof);
    if (res != NULL) delete res;
    res = new double[N_dof];
    if (y_prev != NULL) {
      delete y_prev;
      y_prev = new double[N_dof];
      // FIXME: y_prev should be defined here as the projection
      // of y_prev_ref onto the coarse mesh!
      for(int i=0; i<N_dof; i++) y_prev[i] = 0; 
    }
    else {
      // initially y_prev is seto to zero vector
      y_prev = new double[N_dof];
      for(int i=0; i<N_dof; i++) y_prev[i] = 0; 
    }

    // Obtain coarse mesh solution via Newton's method
    int newton_iterations = 0;
    while (1) {
      // Erase the matrix:
      mat->zero();

      // Construct residual vector
      dp->assemble_matrix_and_vector(mesh, mat, res, y_prev); 

      // Calculate norm of residual vector
      double res_norm = 0;
      for(int i=0; i<N_dof; i++) res_norm += res[i]*res[i];
      res_norm = sqrt(res_norm);

      // If residual norm less than TOL_NEWTON_COARSE, quit
      // latest solution is in y_prev
      printf("Residual norm (coarse): %.15f\n", res_norm);
      if(res_norm < TOL_NEWTON_COARSE) break;

      // Change sign of vector res
      for(int i=0; i<N_dof; i++) res[i]*= -1;

      // Solve the matrix system
      solve_linear_system_umfpack((CooMatrix*)mat, res);

      // Update y_prev by new solution which is in res
      for(int i=0; i<N_dof; i++) y_prev[i] += res[i];

      newton_iterations++;
      printf("Finished coarse mesh Newton iteration: %d\n", newton_iterations);
    }
    // Update y_prev by new solution which is in res
    for(int i=0; i<N_dof; i++) {
      y_prev[i] += res[i];
      //printf("y_prev[%d] = %g\n", i, y_prev[i]);
    }

    // Create reference mesh
    printf("Creating reference mesh.\n");
    if (mesh_ref != NULL) {
      // Adjust the reference mesh according to refinements 
      // that were done in coarse mesh.
      // FIXME: the deletion and replication below is tamporary
      delete mesh_ref;
      mesh_ref = mesh->replicate();
    }
    else {
      // First time: replicate the mesh
      mesh_ref = mesh->replicate();
    }

    // Perform refinements in the reference mesh
    // Refines 'num_to_ref' elements starting with element 'start_elem_id'
    // For now, refine entire mesh uniformly in 'h' and 'p'
    int start_elem_id = 0; 
    int num_to_ref = mesh->get_n_active_elem();
    //mesh_ref->reference_refinement(0, 2);
    mesh_ref->reference_refinement(start_elem_id, num_to_ref);

    // Enumerate DOF in the reference mesh
    int N_dof_ref = mesh_ref->assign_dofs();
    
    // (Re)allocate Jacobi matrix mat_ref and vectors 
    // y_prev_ref and res_ref on reference mesh
    if (mat_ref != NULL) delete mat_ref;
    mat_ref = new CooMatrix(N_dof_ref);
    if (res_ref != NULL) delete res_ref;
    res_ref = new double[N_dof_ref];
    if (y_prev_ref != NULL) delete y_prev_ref;
    y_prev_ref = new double[N_dof_ref];

    // transfer previous reference solution onto the new 
    // reference mesh
    // FIXME: so far the new coarse mesh solution is used
    printf("Transfering solution to reference mesh.\n");
    transfer_solution(mesh, mesh_ref, y_prev, y_prev_ref);

    printf("--- Newton iteration on reference mesh ----\n"); 
    printf("N_dof_ref = %d\n", N_dof_ref);

    // Obtain reference solution via Newton's method
    int newton_iterations_ref = 0;
    while(1) {
      // Zero the matrix:
      mat_ref->zero();

      // Construct residual vector
      dp->assemble_matrix_and_vector(mesh_ref, 
                                   mat_ref, res_ref, y_prev_ref); 

      // Calculate norm of residual vector
      double res_ref_norm = 0;
      for(int i=0; i<N_dof_ref; i++) res_ref_norm += res_ref[i]*res_ref[i];
      res_ref_norm = sqrt(res_ref_norm);

      // If residual norm less than TOL_NEWTON_REF, quit
      // latest solution is in y_prev
      printf("Residual norm (ref): %.15f\n", res_ref_norm);
      if(res_ref_norm < TOL_NEWTON_REF) break;

      // Change sign of vector res_ref
      for(int i=0; i<N_dof_ref; i++) res_ref[i]*= -1;

      // Solve the matrix system
      solve_linear_system_umfpack((CooMatrix*)mat_ref, res_ref);

      // Update y_prev by new solution which is in res
      for(int i=0; i<N_dof_ref; i++) y_prev_ref[i] += res_ref[i];

      newton_iterations_ref++;
      printf("Finished fine mesh Newton iteration: %d\n", newton_iterations_ref);
    }
    // Update y_prev by the increment stored in res
    for(int i=0; i<N_dof_ref; i++) {
      y_prev_ref[i] += res_ref[i];
      //printf("y_prev_ref[%d] = %g\n", i, y_prev_ref[i]);
    }

    // Estimate element errors (squared)
    double err_est_squared_array[MAX_ELEM_NUM]; 
    double err_est_total = calc_elem_est_errors_squared(NORM, mesh, mesh_ref, y_prev, 
                                 y_prev_ref, err_est_squared_array);

    // Calculate the norm of the reference solution
    double ref_sol_norm = calc_approx_sol_norm(NORM, mesh_ref, y_prev_ref);

    // Calculate an estimate of the global relative error
    double err_est_rel = err_est_total/ref_sol_norm;
    printf("Relative error (est) = %g %%\n", 100.*err_est_rel);

    // If exact solution available, also calculate exact error
    if (EXACT_SOL_PROVIDED) {
      // Calculate element errors wrt. exact solution (squared)
      int order = 20; // heuristic parameter
      double err_exact_total = calc_exact_sol_error(NORM, mesh, y_prev, exact_sol, order);
     
      // Calculate the norm of the exact solution
      // (using a fine subdivision and high-order quadrature)
      int subdivision = 500; // heuristic parameter
      double exact_sol_norm = calc_exact_sol_norm(NORM, exact_sol, N_eq, A, B,
                                                  subdivision, order);

       // Calculate an estimate of the global relative error
      double err_exact_rel = err_exact_total/exact_sol_norm;
      printf("Relative error (exact) = %g %%\n", 100.*err_exact_rel);
      graph.add_values(0, N_dof, 100 * err_exact_rel);
    }

    // add entry to DOF convergence graph
    graph.add_values(1, N_dof, 100 * err_est_rel);

    // Decide whether the relative error is sufficiently small
    if(err_est_rel*100 < TOL_ERR_REL) break;

    // extra code for this test:
    if (N_dof > 20) return ERROR_FAILURE;

    // Refine elements in the id_array list whose id_array >= 0
    mesh->adapt(NORM, ADAPT_TYPE, THRESHOLD, mesh_ref, y_prev, 
                y_prev_ref, err_est_squared_array);
    N_dof = mesh->assign_dofs();

    adapt_iterations++;
  };

  printf("Success!\n");
  return ERROR_SUCCESS;
}


















