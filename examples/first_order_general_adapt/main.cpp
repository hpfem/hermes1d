#include "hermes1d.h"
#include "solver_umfpack.h"
#include "adapt.h"

// ********************************************************************

// This example uses automatic hp-adaptivity to solve the general 
// first-order equation y' = f(y, x) in an interval (A, B), equipped 
// with the initial condition y(A) = YA. The function f can be linear
// or nonlinear in 'y', as long as it is differentiable
// with respect to this variable (needed for the Newton's method). 

// General input:
const int N_eq = 1;                     // Number of equations
const int N_elem = 5;                   // Number of elements
const double A = 0, B = 10;             // Domain end points
const double YA = 1;                    // Equation parameter
const int P_init = 1;                   // Initial polynomal degree

// Tolerances for Newton methods
const double TOL_NEWTON_COARSE = 1e-5;  // Coarse mesh
const double TOL_NEWTON_REF = 1e-3;     // Reference mesh

// Adaptivity
const double THRESHOLD = 0.3;           // Refined will be all elements whose error
                                        // is greater than THRESHOLD*max_elem_error
const double TOL_ERR_REL = 1e-5;        // Tolerance for the relative error between 
                                        // the coarse mesh and reference solutions
 
// Right-hand side function f(y, x)
double f(double y, double x) {
  //return -y; // with y(0)=1, exact solution is y=exp(-x)
  return -y*y; // with y(0)=1, exact solution is y=1/(x+1)
}

// y-derivative of dfdy(y, x)
double dfdy(double y, double x) {
  //return -1;
  return -2*y;
}

// ********************************************************************

// Bilinear form for the Jacobi matrix 
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...basis function
// v...test function
// u_prev...previous solution (all solution components)
double jacobian(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dudx[i]*v[i] - dfdy(u_prev[0][i], x[i])*u[i]*v[i])*weights[i];
  }
  return val;
};

// (Nonlinear) form for the residual vector
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...approximate solution
// v...test function
// u_prev...previous solution (all solution components)
double residual(int num, double *x, double *weights, 
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (du_prevdx[0][i]*v[i] - f(u_prev[0][i], x[i])*v[i])*weights[i];
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
  DiscreteProblem *dp = NULL;       // discrete problem (coarse mesh)
  DiscreteProblem *dp_ref = NULL;   // discrete problem (reference mesh)

  // Create coarse mesh, set Dirichlet BC, enumerate basis functions
  mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, YA);
  int N_dof = mesh->assign_dofs();

  // Create discrete problem on coarse mesh
  dp = new DiscreteProblem(mesh);
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
      dp->assemble_matrix_and_vector(mat, res, y_prev); 

      // Calculate L2 norm of residual vector
      double res_norm = 0;
      for(int i=0; i<N_dof; i++) res_norm += res[i]*res[i];
      res_norm = sqrt(res_norm);

      // If residual norm less than TOL_NEWTON_COARSE, quit
      // latest solution is in y_prev
      printf("Residual L2 norm: %.15f\n", res_norm);
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
    for(int i=0; i<N_dof; i++) y_prev[i] += res[i];

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
    
    // Register weak forms
    if (dp_ref != NULL) delete dp_ref;
    dp_ref = new DiscreteProblem(mesh_ref);
    dp_ref->add_matrix_form(0, 0, jacobian);
    dp_ref->add_vector_form(0, residual);

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
      dp_ref->assemble_matrix_and_vector(mat_ref, res_ref, y_prev_ref); 

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
      printf("Finished fine mesh Newton iteration: %d\n", newton_iterations_ref);
    }
    // Update y_prev by the increment stored in res
    for(int i=0; i<N_dof_ref; i++) y_prev_ref[i] += res_ref[i];

    // Calculate element errors (squared)
    double err_squared_array[MAX_ELEM_NUM]; //FIXME - change this to dynamic allocation
    double err_total = calc_elem_L2_errors_squared(mesh, mesh_ref, y_prev, 
                       y_prev_ref, err_squared_array);

    // Calculate the L2 norm of the reference solution
    double ref_sol_L2_norm = calc_solution_L2_norm(mesh_ref, y_prev_ref);

    // Decide whether the relative error is sufficiently small
    double err_rel = err_total/ref_sol_L2_norm;
    printf("Estimated relative error = %g \%\n", 100.*err_rel);
    if(err_rel < TOL_ERR_REL) break;

    // Sort err_squared_array[] in decreasing order
    int id_array[MAX_ELEM_NUM];
    sort_element_errors(mesh->get_n_active_elem(), 
                        err_squared_array, id_array);

    // Decide which elements will be refined
    double max_elem_error = sqrt(err_squared_array[0]);
    for (int i=0; i<mesh->get_n_active_elem(); i++) {
      if(sqrt(err_squared_array[i]) < THRESHOLD*max_elem_error) id_array[i] = -1; 
    }
    
    // Print list of elements to be refined
    printf("Elements to be refined:\n");
    for (int i=0; i<mesh->get_n_active_elem(); i++) {
      if (id_array[i] >= 0) printf("Elem[%d], error = %g\n", id_array[i], 
                                   sqrt(err_squared_array[i]));
    }

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


    // Perform the hp-refinements in the coarse mesh
    //mesh->refine_multi_elems(num_elems_to_be_refined, id_array, p_pair_array);

    adapt_iterations++;

    if (adapt_iterations == 3) break;
  };

  // plotting the coarse mesh solution
  Linearizer l(mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  // plotting the reference solution
  Linearizer lxx(mesh_ref);
  const char *out_filename2 = "solution_ref.gp";
  lxx.plot_solution(out_filename2, y_prev_ref);

  printf("Done.\n");
  return 1;
}
