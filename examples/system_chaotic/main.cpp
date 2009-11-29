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
int DAMPING_STEPS = 20;          // This is the number of damping steps. The entire problem
                                 // will be run repeatedly, with the DAMPING parameter 
                                 // increased from 0 to 1 in DAMPING_STEPS. Every time, 
                                 // the last result will be used as an initial condition for the 
                                 // new computation.   
double DAMPING;                  // DAMPING is an artificial parameter used to reduce the strength 
                                 // of the nonlinearity. (The nonlinearity is multiplied with it.)

// Error tolerance
double TOL_NEWTON = 1e-5;        // tolerance for the Newton's method 

// Boundary conditions
double Val_dir_left_1 = 1;
double Val_dir_left_2 = 0;
double Val_dir_left_3 = 0;
double Val_dir_left_4 = 0;

// ********************************************************************

double jacobian_1_1(int num, double *x, double *weights,
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

  // Allocate vectors res and y_prev
  double *res = new double[N_dof];
  double *y_prev = new double[N_dof];
  if (res == NULL || y_prev == NULL)
    error("res or y_prev could not be allocated in main().");

  // Set y_prev zero
  for(int i=0; i<N_dof; i++) y_prev[i] = 0; 

  // Damping loop
  Mesh *mesh_ref;
  double *y_prev_ref;
  for(int damp_step = 1; damp_step < DAMPING_STEPS+1; damp_step++) {
    DAMPING = sin(damp_step*(1./DAMPING_STEPS)*M_PI/2.);

    printf("Damping: %g\n", DAMPING);

    // Newton's loop
    int newton_iterations = 1;
    CooMatrix *mat = NULL;
    while (1) {
      // Reset the matrix:
      if (mat != NULL) delete mat;
      mat = new CooMatrix();

      // Construct residual vector
      dp->assemble_matrix_and_vector(mesh, mat, res, y_prev); 

      // Calculate L2 norm of residual vector
      double res_norm = 0;
      for(int i=0; i<N_dof; i++) res_norm += res[i]*res[i];
      res_norm = sqrt(res_norm);

      // If residual norm less than TOL_NEWTON, quit
      // latest solution is in y_prev
      printf("Residual L2 norm: %.15f\n", res_norm);
      if(res_norm < TOL_NEWTON) break;

      // Change sign of vector res
      for(int i=0; i<N_dof; i++) res[i]*= -1;

      // Solve the matrix system
      solve_linear_system_umfpack((CooMatrix*)mat, res);

      // Update y_prev by new solution which is in res
      for(int i=0; i<N_dof; i++) y_prev[i] += res[i];

      newton_iterations++;
    }
    printf("Finished Newton loop (%d iter).\n", newton_iterations);
  } // end of the damping loop

  // plotting the solution
  Linearizer l(mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  printf("Done.\n");
  if (y_prev != NULL) delete[] y_prev;
  if (res != NULL) delete[] res;
  return 1;
}
