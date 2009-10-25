#include "hermes1d.h"
#include "solver_umfpack.h"

// ********************************************************************
// This example solves a simplified car control model consisting of
// five first-order equations
// x' - v cos(phi) cos(theta) = 0
// y' - v cos(phi) sin(theta) = 0
// v' - alpha                 = 0
// theta' - v sin(phi)        = 0
// phi' - zeta                = 0
// equipped with the constraints
// |alpha| <= alpha_max
// |zeta| <= zeta_max
// |v| <= v_max
// |phi| <= phi_max

// The problem is considered in a time interval (0, T), and Dirichlet 
// conditions are given for all quantities at both interval endpoints.

// Goal: Calculate the time-dependent controls alpha(t) and zeta(t)
// in such a way that the integral of v^2(t) over the interval (0, T)
// is minimal.

// General input:
static int N_eq = 5;
int N_elem = 3;             // number of elements
double A = 0, B = 10;       // domain end points
int P_init = 3;             // initial polynomal degree 

// Error tolerance
double TOL_NEWTON = 1e-5;        // tolerance for the Newton's method on basic mesh

// Boundary conditions
double X0_left = 0;
double Y0_left = 1;
double V_left = 0;
double THETA_left = 0;
double PHI_left = 0;

// ********************************************************************

double jacobian_0_0(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dudx[i] - DAMPING*(1 - u_prev[1][i]*u_prev[1][i])*u[i])*v[i]*weights[i];
  }
  return val;
};

double jacobian_1_2(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (DAMPING*2*u_prev[0][i]*u_prev[1][i] + 1)*u[i]*v[i]*weights[i];
  }
  return val;
};

double jacobian_2_1(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -u[i]*v[i]*weights[i];
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
    val += u[i]*v[i]*weights[i];
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
    val += -u[i]*v[i]*weights[i];
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

double jacobian_3_4(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += u[i]*v[i]*weights[i];
  }
  return val;
};

double jacobian_4_3(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -u[i]*v[i]*weights[i];
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
    val += (du_prevdx[0][i] - DAMPING*(1 - u_prev[1][i]*u_prev[1][i])*u_prev[0][i]
	    + u_prev[1][i])*v[i]*weights[i];
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
    val += (du_prevdx[1][i] - u_prev[0][i] + u_prev[2][i])*v[i]*weights[i];
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
    val += (du_prevdx[2][i] - u_prev[1][i] + u_prev[3][i])*v[i]*weights[i];
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
    val += (du_prevdx[3][i] - u_prev[2][i])*v[i]*weights[i];
  }
  return val;
};

/******************************************************************************/
int main() {
  // create mesh
  Mesh mesh(A, B, N_elem, P_init, N_eq);
  mesh.set_bc_left_dirichlet(0, Val_dir_left_1);
  mesh.set_bc_left_dirichlet(1, Val_dir_left_2);
  mesh.set_bc_left_dirichlet(2, Val_dir_left_3);
  mesh.set_bc_left_dirichlet(3, Val_dir_left_4);
  int N_dof_basic = mesh.assign_dofs();
  printf("N_dof_basic = %d\n", N_dof_basic);

  // register weak forms
  DiscreteProblem dp(&mesh);
  dp.add_matrix_form(0, 0, jacobian_1_1);
  dp.add_matrix_form(0, 1, jacobian_1_2);
  dp.add_matrix_form(1, 0, jacobian_2_1);
  dp.add_matrix_form(1, 1, jacobian_2_2);
  dp.add_matrix_form(1, 2, jacobian_2_3);
  dp.add_matrix_form(2, 1, jacobian_3_2);
  dp.add_matrix_form(2, 2, jacobian_3_3);
  dp.add_matrix_form(2, 3, jacobian_3_4);
  dp.add_matrix_form(3, 2, jacobian_4_3);
  dp.add_matrix_form(3, 3, jacobian_4_4);
  dp.add_vector_form(0, residual_1);
  dp.add_vector_form(1, residual_2);
  dp.add_vector_form(2, residual_3);
  dp.add_vector_form(3, residual_4);

  // Allocate Jacobi matrix and residual
  Matrix *mat = new CooMatrix(N_dof_basic);
  double *y_prev = new double[N_dof_basic];
  double *res = new double[N_dof_basic];

  // Set zero initial condition for the Newton's method
  // on the basic mesh
  for(int i=0; i<N_dof_basic; i++) y_prev[i] = 0; 

  // Damping loop
  Mesh *mesh_ref;
  double *y_prev_ref;
  printf("Damping: %g\n", DAMPING);
  printf("------------- Newton's iterations on basic mesh -------------- \n"); 

  // Newton's loop on coarse mesh
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

    // Update y_prev by new solution which is in res
    for(int i=0; i<N_dof_basic; i++) y_prev[i] += res[i];
  } // end of adaptivity loop

  // plotting the basic solution
  Linearizer l(&mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  printf("Done.\n");
  if (y_prev != NULL) delete[] y_prev;
  if (res != NULL) delete[] res;
  if (mat != NULL) delete mat;
  return 1;
}
