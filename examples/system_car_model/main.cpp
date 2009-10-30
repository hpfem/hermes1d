#include "hermes1d.h"
#include "solver_umfpack.h"

// ********************************************************************
// This example solves a simplified car control model consisting of
// five first-order equations
//
// x' - v cos(phi) cos(theta) = 0
// y' - v cos(phi) sin(theta) = 0
// v' - alpha                 = 0
// phi' - zeta                = 0
// theta' - v sin(phi)        = 0
//
// equipped with the constraints
//
// |alpha| <= alpha_max
// |zeta| <= zeta_max
// |v| <= v_max
// |phi| <= phi_max

// The problem is considered in a time interval (0, T), and Dirichlet 
// conditions are given for all quantities at the beginning.

// Goal: Calculate all possible trajectories of the car given the 
// initial condition and intervals for alpha and zeta.

// Print data ?
const int PRINT = 0;

// General input:
const int N_eq = 5;
const int N_elem = 3;              // number of elements
const double A = 0, B = 1;         // domain end points
const int P_init = 3;              // initial polynomal degree 

// Parameters
const double Alpha_max = 1.;
const double Zeta_max = M_PI/6;
const int Num_rays = 100;
const int N_steps_per_ray = 6;    // NOTE: be careful with this number,
                                  // there are four embedded loops of
                                  // this length!

// Boundary conditions
const double X0_left = 0;
const double Y0_left = 0;
const double Vel_left = 0;
const double Phi_left = 0;
const double Theta_left = 0;

// Controls
const int N_ctrl = 4;
double time_ctrl[N_ctrl] = 
  {A, A + (A+B)/(N_ctrl-1.), A + 2.*(A+B)/(N_ctrl-1.), B};
double alpha_ctrl[N_ctrl] = {0, 0, 0, 0};
double zeta_ctrl[N_ctrl] = {0, 0, 0, 0};

// Error tolerance
const double TOL_newton = 1e-5;        // tolerance for the Newton's method

// ********************************************************************

double get_ctrl_alpha(double t) 
{
  // FIXME: this implementation is not efficient
  for (int i = 0; i<N_ctrl-1; i++) {
    if (time_ctrl[i] <= t && t <= time_ctrl[i+1]) {
      // return linear interpolant
      return alpha_ctrl[i] + 
	(t - time_ctrl[i])*(alpha_ctrl[i+1] - alpha_ctrl[i])/(time_ctrl[i+1] - time_ctrl[i]);
    }
  }
  error("Internal: time interval not found in get_ctrl_alpha().");
}

double get_ctrl_zeta(double t) 
{
  // FIXME: this implementation is not efficient
  for (int i = 0; i<N_ctrl-1; i++) {
    if (time_ctrl[i] <= t && t <= time_ctrl[i+1]) {
      // return linear interpolant
      return zeta_ctrl[i] + 
	(t - time_ctrl[i])*(zeta_ctrl[i+1] - zeta_ctrl[i])/(time_ctrl[i+1] - time_ctrl[i]);
    }
  }
  error("Internal: time interval not found in get_ctrl_zeta().");
}


// ********************************************************************

double jacobian_0_0(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_0_2(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  // renaming for better readability
  double* phi = u_prev[3];
  double* theta = u_prev[4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -u[i] * cos(phi[i]) * cos(theta[i]) * v[i] * weights[i];
  }
  return val;
};

double jacobian_0_3(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  // renaming for better readability
  double* vel = u_prev[2];
  double* phi = u_prev[3];
  double* theta = u_prev[4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += vel[i] * sin(phi[i]) * u[i] * cos(theta[i]) * v[i] * weights[i];
  }
  return val;
};

double jacobian_0_4(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  // renaming for better readability
  double* vel = u_prev[2];
  double* phi = u_prev[3];
  double* theta = u_prev[4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += vel[i] * cos(phi[i]) * sin(theta[i]) * u[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_1_1(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_1_2(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  // renaming for better readability
  double* phi = u_prev[3];
  double* theta = u_prev[4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -u[i] * cos(phi[i]) * sin(theta[i]) * v[i] * weights[i];
  }
  return val;
};

double jacobian_1_3(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  // renaming for better readability
  double* vel = u_prev[2];
  double* phi = u_prev[3];
  double* theta = u_prev[4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += vel[i] * sin(phi[i]) * u[i] * sin(theta[i]) * v[i] * weights[i];
  }
  return val;
};

double jacobian_1_4(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  // renaming for better readability
  double* vel = u_prev[2];
  double* phi = u_prev[3];
  double* theta = u_prev[4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -vel[i] * cos(phi[i]) * cos(theta[i]) * u[i] * v[i] * weights[i];
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
    val += dudx[i] * v[i] * weights[i];
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
    val += dudx[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_4_2(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  // renaming for better readability
  double* phi = u_prev[3];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -u[i] * sin(phi[i]) * v[i] * weights[i];
  }
  return val;
};

double jacobian_4_3(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                void *user_data)
{
  // renaming for better readability
  double* vel = u_prev[2];
  double* phi = u_prev[3];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -vel[i] * cos(phi[i]) * u[i] * v[i] * weights[i];
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
    val += dudx[i] * v[i] * weights[i];
  }
  return val;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double residual_0(int num, double *x, double *weights,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
  // renaming for better readability
  double* dxposdt = du_prevdx[0];
  double* vel = u_prev[2];
  double* phi = u_prev[3];
  double* theta = u_prev[4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dxposdt[i] - vel[i] * cos(phi[i]) * cos(theta[i])) * v[i] * weights[i];
  }
  return val;
};

double residual_1(int num, double *x, double *weights,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
  // renaming for better readability
  double* dyposdt = du_prevdx[1];
  double* vel = u_prev[2];
  double* phi = u_prev[3];
  double* theta = u_prev[4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dyposdt[i] - vel[i] * cos(phi[i]) * sin(theta[i])) * v[i] * weights[i];
  }
  return val;
};

double residual_2(int num, double *x, double *weights,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
  // renaming for better readability
  double* dveldt = du_prevdx[2];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dveldt[i] - get_ctrl_alpha(x[i])) * v[i] * weights[i];
  }
  return val;
};


double residual_3(int num, double *x, double *weights,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
  // renaming for better readability
  double* dphidt = du_prevdx[3];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dphidt[i] - get_ctrl_zeta(x[i])) * v[i] * weights[i];
  }
  return val;
};

double residual_4(int num, double *x, double *weights,
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
  // renaming for better readability
  double* vel = u_prev[2];
  double* phi = u_prev[3];
  double* dthetadt = du_prevdx[4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dthetadt[i] - vel[i] * sin(phi[i])) * v[i] * weights[i];
  }
  return val;
};

/******************************************************************************/
int main() {
  // create mesh
  Mesh mesh(A, B, N_elem, P_init, N_eq);
  mesh.set_bc_left_dirichlet(0, X0_left);
  mesh.set_bc_left_dirichlet(1, Y0_left);
  mesh.set_bc_left_dirichlet(2, Vel_left);
  mesh.set_bc_left_dirichlet(3, Phi_left);
  mesh.set_bc_left_dirichlet(4, Theta_left);
  int N_dof = mesh.assign_dofs();
  printf("N_dof = %d\n", N_dof);

  // register weak forms
  DiscreteProblem dp(&mesh);
  dp.add_matrix_form(0, 0, jacobian_0_0);
  dp.add_matrix_form(0, 2, jacobian_0_2);
  dp.add_matrix_form(0, 3, jacobian_0_3);
  dp.add_matrix_form(0, 4, jacobian_0_4);
  dp.add_matrix_form(1, 1, jacobian_1_1);
  dp.add_matrix_form(1, 2, jacobian_1_2);
  dp.add_matrix_form(1, 3, jacobian_1_3);
  dp.add_matrix_form(1, 4, jacobian_1_4);
  dp.add_matrix_form(2, 2, jacobian_2_2);
  dp.add_matrix_form(3, 3, jacobian_3_3);
  dp.add_matrix_form(4, 2, jacobian_4_2);
  dp.add_matrix_form(4, 3, jacobian_4_3);
  dp.add_matrix_form(4, 4, jacobian_4_4);
  dp.add_vector_form(0, residual_0);
  dp.add_vector_form(1, residual_1);
  dp.add_vector_form(2, residual_2);
  dp.add_vector_form(3, residual_3);
  dp.add_vector_form(4, residual_4);

  // Allocate Jacobi matrix and residual
  Matrix *mat = new CooMatrix(N_dof);
  double *y_prev = new double[N_dof];
  double *res = new double[N_dof];

  // Loop over rays
  double radius = sqrt(Alpha_max*Alpha_max + Zeta_max*Zeta_max)/2.;
  for (int a = 0; a < Num_rays; a++) {
    double ray_angle = 2*M_PI*a/Num_rays;
    printf("Ray %d, angle %g\n", a, ray_angle);
    double alpha_actual = radius * cos(ray_angle);
    if (alpha_actual > Alpha_max) alpha_actual = Alpha_max;
    if (alpha_actual < -Alpha_max) alpha_actual = -Alpha_max;
    double zeta_actual = radius * sin(ray_angle); 
    if (zeta_actual > Zeta_max) zeta_actual = Zeta_max;
    if (zeta_actual < -Zeta_max) zeta_actual = -Zeta_max;

    // At the beginning of every ray: Set zero initial condition 
    // for the Newton's method
    for(int i=0; i<N_dof; i++) y_prev[i] = 0; 

    // Loop over control parameters alpha and zeta lying on the ray. 
    // Start from the origin and move towards the boundary of
    // the rectangle
    for (int step0 = 0; step0 < N_steps_per_ray; step0++) {
      alpha_ctrl[0] = alpha_actual * step0/double(N_steps_per_ray);      
      zeta_ctrl[0] = zeta_actual * step0/double(N_steps_per_ray);      
      for (int step1 = 0; step1 < N_steps_per_ray; step1++) {
        alpha_ctrl[1] = alpha_actual * step1/double(N_steps_per_ray);      
        zeta_ctrl[1] = zeta_actual * step1/double(N_steps_per_ray);      
        for (int step2 = 0; step2 < N_steps_per_ray; step2++) {
          alpha_ctrl[2] = alpha_actual * step2/double(N_steps_per_ray);      
          zeta_ctrl[2] = zeta_actual * step2/double(N_steps_per_ray);      
          for (int step3 = 0; step3 < N_steps_per_ray; step3++) {
            alpha_ctrl[3] = alpha_actual * step3/double(N_steps_per_ray);
            zeta_ctrl[3] = zeta_actual * step3/double(N_steps_per_ray);      

            // Newton's iteration
            int newton_iterations = 0;
              
            if (PRINT) {
              printf("------------- Newton's iterations -------------- \n"); 
              printf("alpha = (%g, %g, %g, %g), zeta = (%g, %g, %g, %g)\n", 
                     alpha_ctrl[0], alpha_ctrl[1], 
                     alpha_ctrl[2], alpha_ctrl[3], zeta_ctrl[0], 
                     zeta_ctrl[1], zeta_ctrl[2], zeta_ctrl[3]); 
            }
            while (1) {
              // Erase the matrix:
              mat->zero();

              // Construct residual vector
              dp.assemble_matrix_and_vector(mat, res, y_prev); 

              // Calculate L2 norm of residual vector
              double res_norm = 0;
              for(int i=0; i<N_dof; i++) res_norm += res[i]*res[i];
              res_norm = sqrt(res_norm);

              // If residual norm less than TOL_newton, quit
              // latest solution is in y_prev
              if (PRINT) printf("Residual L2 norm: %.15f\n", res_norm);
              if(res_norm < TOL_newton) break;

              // Change sign of vector res
              for(int i=0; i<N_dof; i++) res[i]*= -1;

              // Solve the matrix system
              solve_linear_system_umfpack((CooMatrix*)mat, res);

              newton_iterations++;
              if (PRINT) printf("Finished Newton iteration: %d\n", newton_iterations);

              // Update y_prev by new solution which is in res
              for(int i=0; i<N_dof; i++) y_prev[i] += res[i];
            } // end of Newton's loop

            // plotting the solution
            //Linearizer l(&mesh);
            //const char *out_filename = "solution.gp";
            //l.plot_solution(out_filename, y_prev, 10);

            // plotting the trajectory
            static int first_traj = 1;
            const char *traj_filename = "trajectory.gp";
            FILE *f;
            if (first_traj) {
              f = fopen(traj_filename, "wb");
              first_traj = 0;
            }
            else f = fopen(traj_filename, "ab");
            if(f == NULL) error("Problem opening trajectory file.");
            Linearizer l_traj(&mesh);
            l_traj.plot_trajectory(f, y_prev, 0, 1, 10);
            static int traj_count = 0; 
            traj_count++;
            if (PRINT) printf("Trajectory %d written to %s.\n", 
              traj_count, traj_filename);
            fclose(f);
	  }
	}
      }
    }
  }

  printf("Done.\n");
  if (y_prev != NULL) delete[] y_prev;
  if (res != NULL) delete[] res;
  if (mat != NULL) delete mat;
  return 1;
}
