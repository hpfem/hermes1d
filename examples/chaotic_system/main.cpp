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

// Disable reference solution and adaptivity
int ADAPTIVITY_DISABLED = 1;     

// Damping parameter
int DAMPING_STEPS = 20;          // This is the number of damping steps. The entire problem
                                 // will be run repeatedly, with the DAMPING parameter 
                                 // increased from 0 to 1 in DAMPING_STEPS. Every time, 
                                 // the last result will be used as an initial condition for the 
                                 // new computation.   
double DAMPING;                  // DAMPING is an artificial parameter used to reduce the strength 
                                 // of the nonlinearity. (The nonlinearity is multiplied with it.)

// Error tolerance
double TOL_NEWTON_BASIC = 1e-5;  // tolerance for the Newton's method on basic mesh
double TOL_NEWTON_REF = 1e-4;    // tolerance for the Newton's method on reference mesh
double TOL_ADAPT = 1e-5;         // tolerance for the adaptivity loop      

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
  for(int damp_step = 1; damp_step < DAMPING_STEPS+1; damp_step++) {
    DAMPING = sin(damp_step*(1./DAMPING_STEPS)*M_PI/2.);

    // Adaptivity loop
    int adapt_iterations = 1;
    while(1) {

      printf("------------------ Adaptivity iteration %d --------------------\n", 
	     adapt_iterations); 
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
      }
      // Update y_prev by new solution which is in res
      for(int i=0; i<N_dof_basic; i++) y_prev[i] += res[i];

      if(ADAPTIVITY_DISABLED == 0) {
        printf("--------- Starting Newton iterations on reference mesh (damping %g) ---------\n", DAMPING); 

        // Replicate current mesh
        mesh_ref = mesh.replicate();

        // Refines 'num_to_ref' elements starting with element 'start_elem_id'
        // For now, refine entire mesh uniformly in 'h' and 'p'
        int start_elem_id = 0; 
        int num_to_ref = mesh.get_n_active_elem();
        //mesh_ref->refine_elems(0, 1);
        static int first_time_refine_elems = 1;
        if(first_time_refine_elems == 1) {
          mesh_ref->refine_elems(start_elem_id, num_to_ref);
          first_time_refine_elems = 0;
        }

        // Enumerate DOF in the reference mesh
        int N_dof_ref = mesh_ref->assign_dofs();
        printf("N_dof_ref = %d\n", N_dof_ref);

        // Register weak forms
        // FIXME: the main damping loop is very much messed up, needs thorough cleaning !!!
        static int first_time_dp_ref = 1;
        DiscreteProblem dp_ref(mesh_ref);
        if(first_time_dp_ref == 1) {
          dp_ref.add_matrix_form(0, 0, jacobian_1_1);
          dp_ref.add_matrix_form(0, 1, jacobian_1_2);
          dp_ref.add_matrix_form(1, 0, jacobian_2_1);
          dp_ref.add_matrix_form(1, 1, jacobian_2_2);
          dp_ref.add_matrix_form(1, 2, jacobian_2_3);
          dp_ref.add_matrix_form(2, 1, jacobian_3_2);
          dp_ref.add_matrix_form(2, 2, jacobian_3_3);
          dp_ref.add_matrix_form(2, 3, jacobian_3_4);
          dp_ref.add_matrix_form(3, 2, jacobian_4_3);
          dp_ref.add_matrix_form(3, 3, jacobian_4_4);
          dp_ref.add_vector_form(0, residual_1);
          dp_ref.add_vector_form(1, residual_2);
          dp_ref.add_vector_form(2, residual_3);
          dp_ref.add_vector_form(3, residual_4);
          first_time_dp_ref = 0;
        }

        // Allocate Jacobi matrix and residual
        // FIXME: this should not be called every time!
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

        // Newton's loop on fine mesh
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

        // feeing matrices and vectors (N_dof_ref will change)
        if (y_prev_ref != NULL) delete[] y_prev_ref;
        if (res_ref != NULL) delete[] res_ref;
        if (mat_ref != NULL) delete mat_ref;

        adapt_iterations++;
      } // end of condition ADAPTIVITY_DISABLED == 0
      else break;
    } // end of the adaptivity loop
  } // end of the damping loop

  // plotting the basic solution
  Linearizer l(&mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  // plotting the reference solution
  //Linearizer lxx(mesh_ref);
  //const char *out_filename2 = "solution_ref.gp";
  //lxx.plot_solution(out_filename2, y_prev_ref);
 
  printf("Done.\n");
  if (y_prev != NULL) delete[] y_prev;
  if (res != NULL) delete[] res;
  if (mat != NULL) delete mat;
  return 1;
}
