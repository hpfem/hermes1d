#include "hermes1d.h"

// ********************************************************************
// This example shows solution of harmonic steady state on the homoegenous transmision line.
// Wave propagation is described by two linear differencial equation with complex coefficients.
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

// ********************************************************************

void plotting(Mesh *mesh, Mesh *mesh_ref, double *y_prev, double *y_prev_ref) 
{
  // Plot the coarse mesh solution
  Linearizer l(mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  // Plot the reference solution
  Linearizer lxx(mesh_ref);
  const char *out_filename2 = "solution_ref.gp";
  lxx.plot_solution(out_filename2, y_prev_ref);

  // Plot the coarse and reference mesh
  const char *mesh_filename = "mesh.gp";
  mesh->plot(mesh_filename);
  const char *mesh_ref_filename = "mesh_ref.gp";
  mesh_ref->plot(mesh_ref_filename);

  // Plot the error estimate (difference between 
  // coarse and reference mesh solutions)
  const char *err_est_filename = "error_est.gp";
  mesh->plot_error_est(NORM, err_est_filename, mesh_ref, y_prev, y_prev_ref);

  // Plot error wrt. exact solution (if available)
  if (EXACT_SOL_PROVIDED) {   
    const char *err_exact_filename = "error_exact.gp";
    mesh->plot_error_exact(NORM, err_exact_filename, y_prev, exact_sol);
  }
}

// ********************************************************************

double jacobian_1_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_1_3(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -R*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_1_4(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += +omega*L*u[i]*v[i]*weights[i];
    }
    return val;
};


double jacobian_2_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_2_3(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -omega*L*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_2_4(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -R*u[i]*v[i]*weights[i];
    }
    return val;
};


double jacobian_3_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -G*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_3_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += +omega*C*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_3_3(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};



double jacobian_4_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -omega*C*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_4_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -G*weights[i];
    }
    return val;
};

double jacobian_4_4(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double residual_1(int num, double *x, double *weights,
                  double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                  double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[0][i] - R*u_prev[2][i] +L*omega*u_prev[3][i])*v[i]*weights[i];
    }
    return val;
};

double residual_2(int num, double *x, double *weights,
                  double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                  double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[1][i] - R*u_prev[3][i]-omega*L*u_prev[2][i])*v[i]*weights[i];
    }
    return val;
};

double residual_3(int num, double *x, double *weights,
                  double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                  double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[2][i] - G*u_prev[0][i] + omega*C*u_prev[1][i])*v[i]*weights[i];
    }
    return val;
};


double residual_4(int num, double *x, double *weights,
                  double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                  double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[3][i] - G*u_prev[1][i]-omega*C*u_prev[0][i])*v[i]*weights[i];
    }
    return val;
};


double jacobian_surf_right_U_Re(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_EQN_NUM],
                                double du_prevdx[MAX_EQN_NUM], void *user_data)
{     
    return 1*u*v;
}


double jacobian_surf_right_U_Im(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_EQN_NUM],
                                double du_prevdx[MAX_EQN_NUM], void *user_data)
{
    return -Zl*u*v;
}

double jacobian_surf_right_I_Re(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_EQN_NUM],
                                double du_prevdx[MAX_EQN_NUM], void *user_data)
{
    return 1*u*v;
}


double jacobian_surf_right_I_Im(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_EQN_NUM],
                                double du_prevdx[MAX_EQN_NUM], void *user_data)
{
    return -Zl*u*v;
}

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
  graph.add_row("error estimate", "k", "--");

  // Create coarse mesh, set Dirichlet BC, enumerate basis functions
  mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left_1);
  mesh->set_bc_left_dirichlet(1, Val_dir_left_2);
  mesh->set_bc_left_dirichlet(2, Val_dir_left_3);
  mesh->set_bc_left_dirichlet(3, Val_dir_left_4);
  int N_dof = mesh->assign_dofs();
  printf("N_dof = %d\n", N_dof);

  // Create discrete problem on coarse mesh
  dp = new DiscreteProblem();
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

    // Obtain coarse mesh solution via Newton's iteration
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
      printf("Residual norm: %.15f\n", res_norm);
      if(res_norm < TOL_NEWTON_COARSE) break;

      // Change sign of vector res
      for(int i=0; i<N_dof; i++) res[i]*= -1;

      // Solve the matrix system
      solve_linear_system_umfpack((CooMatrix*)mat, res);

      // Update y_prev by new solution which is in res
      for(int i=0; i<N_dof; i++) y_prev[i] += res[i];

      newton_iterations++;
      printf("Finished coarse Newton iteration: %d\n", newton_iterations);
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
    
    // (Re)allocate Jacobi matrix mat_ref and vectors 
    // y_prev_ref and res_ref on reference mesh
    Matrix *mat_ref = new CooMatrix(N_dof_ref);
    y_prev_ref = new double[N_dof_ref];
    double *res_ref = new double[N_dof_ref];

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
      dp->assemble_matrix_and_vector(mesh_ref, mat_ref, res_ref, y_prev_ref); 

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
    // Update y_prev_ref by the increment stored in res
    for(int i=0; i<N_dof_ref; i++) y_prev_ref[i] += res_ref[i];

    // Estimate element errors (squared)
    double err_est_squared_array[MAX_ELEM_NUM]; 
    double err_est_total = calc_elem_est_errors_squared(NORM, 
                                 mesh, mesh_ref, y_prev, 
                                 y_prev_ref, err_est_squared_array);

    // Calculate the norm of the reference solution
    double ref_sol_norm = calc_approx_sol_norm(NORM, mesh_ref, y_prev_ref);

    // Calculate an estimate of the global relative error
    double err_est_rel = err_est_total/ref_sol_norm;
    printf("Relative error (est) = %g %%\n", 100.*err_est_rel);

    // add entry to DOF convergence graph
    graph.add_values(0, N_dof, 100 * err_est_rel);

    // Decide whether the relative error is sufficiently small
    if(err_est_rel*100 < TOL_ERR_REL) break;

    // debug
    if (adapt_iterations == 1) break;

    // Refine elements in the id_array list whose id_array >= 0
    mesh->adapt(NORM, ADAPT_TYPE, THRESHOLD, mesh_ref, y_prev, 
                y_prev_ref, err_est_squared_array);
    N_dof = mesh->assign_dofs();

    adapt_iterations++;
  };

  // Plot meshes, results, and errors
  plotting(mesh, mesh_ref, y_prev, y_prev_ref);

  // Save convergence graph
  graph.save("conv_dof.gp");

  printf("Done.\n");
  return 1;
}
