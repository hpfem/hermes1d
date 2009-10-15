#include "hermes1d.h"
#include "solver_umfpack.h"

// ********************************************************************

// This example solves the Poisson equation -u'' - f = 0 in
// an interval (A, B), equipped with a Dirichlet boundary
// condition on the left and a Newton BC "u + \alpha dudn = \beta" 
// on the right. 

// General input:
static int N_eq = 1;
int N_elem = 30;                       // number of elements
double A = 0, B = 2*M_PI;              // domain end points
int P_init = 2;                        // initial polynomal degree

// Boundary conditions
double Val_dir_left = 0;
double Val_newton_alpha = 1;   // must be nonzero
double Val_newton_beta = 0.5;

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
// u_prev...previous solution
double jacobian_vol(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],  
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i]*dvdx[i]*weights[i];
  }
  return val;
};

double residual_vol(int num, double *x, double *weights, 
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (du_prevdx[0][i]*dvdx[i] - f(x[i])*v[i])*weights[i];
  }
  return val;
};

double jacobian_surf_right(double x, double u, double dudx,
        double v, double dvdx, double u_prev[MAX_EQN_NUM], 
        double du_prevdx[MAX_EQN_NUM], void *user_data)
{
  return (1./Val_newton_alpha)*u*v;
}

double residual_surf_right(double x, double u_prev[MAX_EQN_NUM], 
        double du_prevdx[MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
  return -(Val_newton_beta/Val_newton_alpha) * v; 
}

/******************************************************************************/
int main() {
  // create mesh
  Mesh mesh(A, B, N_elem, P_init, N_eq);

  // boundary conditions
  mesh.set_bc_left_dirichlet(0, Val_dir_left);
  mesh.set_bc_right_natural(0);
  int N_dof = mesh.assign_dofs();
  printf("N_dof = %d\n", N_dof);

  // register weak forms
  DiscreteProblem dp(&mesh);
  dp.add_matrix_form(0, 0, jacobian_vol);
  dp.add_vector_form(0, residual_vol);
  dp.add_matrix_form_surf(0, 0, jacobian_surf_right, BOUNDARY_RIGHT);
  dp.add_vector_form_surf(0, residual_surf_right, BOUNDARY_RIGHT);

  // allocate Jacobi matrix and residual
  Matrix *mat;
  double *y_prev = new double[N_dof];
  double *res = new double[N_dof];

  // zero initial condition for the Newton's method
  for(int i=0; i<N_dof; i++) y_prev[i] = 0; 

  int newton_iterations = 1;
  // Newton's loop
  while (1) {
    // zero the matrix:
    mat = new DenseMatrix(N_dof);

    // construct residual vector
    dp.assemble_matrix_and_vector(mat, res, y_prev); 

    // calculate L2 norm of residual vector
    double res_norm = 0;
    for(int i=0; i<N_dof; i++) res_norm += res[i]*res[i];
    res_norm = sqrt(res_norm);

    // if residual norm less than TOL, quit
    // latest solution is in y_prev
    printf("Residual L2 norm: %.15f\n", res_norm);
    if(res_norm < TOL) break;

    // changing sign of vector res
    for(int i=0; i<N_dof; i++) res[i]*= -1;

    //mat->print();

    // solving the matrix system
    solve_linear_system(mat, res);

    // updating y_prev by new solution which is in res
    for(int i=0; i<N_dof; i++) y_prev[i] += res[i];
    printf("Finished Newton iteration: %d\n", newton_iterations);
    newton_iterations++;
  }

  Linearizer l(&mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  printf("Done.\n");
  return 1;
}
