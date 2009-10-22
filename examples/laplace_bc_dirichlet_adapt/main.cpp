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
int P_init = 2;                         // initial polynomal degree

// Error tolerance
double TOL_NEWTON_BASIC = 1e-5;  // tolerance for the Newton's method on basic mesh
double TOL_NEWTON_REF = 1e-5;    // tolerance for the Newton's method on reference mesh
double TOL_ADAPT = 1e-5;         // tolerance for the adaptivity loop      

// Boundary conditions
double Val_dir_left = 1;                // Dirichlet condition left
double Val_dir_right = 1;               // Dirichlet condition right

// Tolerance for the Newton's method
double TOL = 1e-5;

// Function f(x)
double f(double x) {
  return sin(x);
  //return 1;
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
  // create mesh
  Mesh mesh(A, B, N_elem, P_init, N_eq);
  mesh.set_bc_left_dirichlet(0, Val_dir_left);
  mesh.set_bc_right_dirichlet(0, Val_dir_right);
  int N_dof_basic = mesh.assign_dofs();
  printf("N_dof_basic = %d\n", N_dof_basic);

  // register weak forms
  DiscreteProblem dp(&mesh);
  dp.add_matrix_form(0, 0, jacobian);
  dp.add_vector_form(0, residual);

  // Allocate Jacobi matrix and residual
  Matrix *mat = new CooMatrix(N_dof_basic);
  double *y_prev = new double[N_dof_basic];
  double *res = new double[N_dof_basic];

  // Zero initial condition for the Newton's method
  // on the basic mesh
  for(int i=0; i<N_dof_basic; i++) y_prev[i] = 0; 

  // Main adaptivity loop
  int adapt_iterations = 1;
  Mesh *mesh_ref;
  double *y_prev_ref;
  while(1) {
 
    printf("------------------ Adaptivity iteration %d --------------------\n", 
           adapt_iterations); 
    printf("----------- Starting Newton iterations on basic mesh ----------\n"); 

    // Obtain basic mesh solution via Newton's iteration
    int newton_iterations = 0;
    while (1) {
      // Erase the matrix:
      mat->zero();

      // Construct residual vector
      dp.assemble_matrix_and_vector(mat, res, y_prev); 

      // Calculate L2 norm of residual vector
      double res_norm = 0;
      for(int i=0; i<N_dof_basic; i++) res_norm += res[i]*res[i];
      res_norm = sqrt(res_norm);

      // If residual norm less than TOL_NEWTON_BASIC, quit
      // latest solution is in y_prev
      printf("Residual L2 norm: %.15f\n", res_norm);
      if(res_norm < TOL_NEWTON_BASIC) break;

      // Change sign of vector res
      for(int i=0; i<N_dof_basic; i++) res[i]*= -1;

      // Solve the matrix system
      solve_linear_system_umfpack((CooMatrix*)mat, res);

      // Update y_prev by new solution which is in res
      for(int i=0; i<N_dof_basic; i++) y_prev[i] += res[i];

      newton_iterations++;
      printf("Finished coarse Newton iteration: %d\n", newton_iterations);
    }
    // Update y_prev by new solution which is in res
    for(int i=0; i<N_dof_basic; i++) y_prev[i] += res[i];

    printf("--------- Starting Newton iterations on reference mesh  ---------\n"); 

    // Replicate current mesh
    mesh_ref = mesh.replicate();

    // Refines 'num_to_ref' elements starting with element 'start_elem_id'
    // For now, refine entire mesh uniformly in 'h' and 'p'
    int start_elem_id = 0; 
    int num_to_ref = mesh.get_n_active_elem();
    //mesh_ref->refine_elems(0, 1);
    mesh_ref->refine_elems(start_elem_id, num_to_ref);

    // Enumerate DOF in the reference mesh
    int N_dof_ref = mesh_ref->assign_dofs();
    printf("N_dof_ref = %d\n", N_dof_ref);
    
    // Register weak forms
    DiscreteProblem dp_ref(mesh_ref);
    dp_ref.add_matrix_form(0, 0, jacobian);
    dp_ref.add_vector_form(0, residual);

    // Allocate Jacobi matrix and residual
    Matrix *mat_ref = new CooMatrix(N_dof_ref);
    y_prev_ref = new double[N_dof_ref];
    double *res_ref = new double[N_dof_ref];


    /*
    for (int i=0; i < N_dof_ref; i++)
        y_prev_ref[i] = 0;
        */
    transfer_solution(&mesh, mesh_ref, y_prev, y_prev_ref);
    /*
    for (int i=0; i < N_dof_ref; i++)
        printf("y_prev_ref[%d] = %f\n", i, y_prev_ref[i]);
        */

    // Use the basic solution as the initial condition 
    // for the Newton's iteration for the reference solution
    //mesh_ref.project(mesh, y_prev, 0, mesh.get_n_active_elems(), y_ref_prev);

    // Obtain basic mesh solution via Newton's iteration
    int newton_iterations_ref = 0;
    while(1) {
      // Zero the matrix:
      mat_ref->zero();

      // Construct residual vector
      dp_ref.assemble_matrix_and_vector(mat_ref, res_ref, y_prev_ref); 

      // Calculate L2 norm of residual vector
      double res_ref_norm = 0;
      for(int i=0; i<N_dof_ref; i++) res_ref_norm += res_ref[i]*res_ref[i];
      res_ref_norm = sqrt(res_ref_norm);

      // If residual norm less than TOL_NEWTON_REF, quit
      // latest solution is in y_prev
      printf("ref: Residual L2 norm: %.15f\n", res_ref_norm);
      if(res_ref_norm < TOL_NEWTON_REF) break;

      // Change sign of vector res_ref
      for(int i=0; i<N_dof_ref; i++) res_ref[i]*= -1;

      // Solve the matrix system
      solve_linear_system_umfpack((CooMatrix*)mat_ref, res_ref);

      // Update y_prev by new solution which is in res
      for(int i=0; i<N_dof_ref; i++) y_prev_ref[i] += res_ref[i];

      newton_iterations_ref++;
      printf("Finished coarse Newton iteration: %d\n", newton_iterations_ref);
    }
    // Update y_prev by new solution which is in res
    for(int i=0; i<N_dof_ref; i++) y_prev_ref[i] += res_ref[i];

    //debugging prints for y_prev:
    //printf("y_prev_ref = \n");
    //for(int i=0; i<N_dof_ref; i++) printf("%g, ", y_prev_ref[i]);
    //printf("\n");

    /* TO BE USED FOR ADAPTIVITY
    // Use the difference between the two solutions to determine 
    // list of elements to be refined (plus the polynomial degrees
    // on sons). Also calculate the total error estimate. 
    unsigned num_elems_to_be_refined; 
    int id_array[1000];
    int2 p_pair_array[1000];
    double err_adapt = 0;
    ...
    */

    // Check whether error is less than tolerance
    double err_adapt = 0; 
    if(err_adapt < TOL_ADAPT) break;

    // Perform the hp-refinements in the basic mesh
    //mesh.refine_multi_elems(num_elems_to_be_refined, id_array, p_pair_array);

    adapt_iterations++;
  };

  // plotting the basic solution
  Linearizer l(&mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  // plotting the reference solution
  Linearizer lxx(mesh_ref);
  const char *out_filename2 = "solution_ref.gp";
  lxx.plot_solution(out_filename2, y_prev_ref);

  printf("Done.\n");
  return 1;
}
