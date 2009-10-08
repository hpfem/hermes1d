#include "hermes1d.h"
#include "solver_umfpack.h"

// ********************************************************************

// general input:
static int N_eq = 1;
int N_elem = 30;                         // number of elements
double A = 0, B = 2*M_PI;              // domain end points
int P_init = 1;                        // initial polynomal degree

// boundary conditions
double val_newton_alpha_left = 2;
double val_newton_beta_left = -2;
double val_newton_alpha_right = 1;
double val_newton_beta_right = 1;

// Tolerance for Newton's method
double TOL = 1e-5;

// right-hand side
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
double jacobian(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double *u_prev, double *du_prevdx, void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i]*dvdx[i]*weights[i];
  }
  return val;
};

double residual_vol(int num, double *x, double *weights, 
                double *u_prev, double *du_prevdx, double *v, double *dvdx,
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (du_prevdx[i]*dvdx[i] - f(x[i])*v[i])*weights[i];
  }
  if(DEBUG) {
    /*printf("u = ");
    for(int i=0; i<num; i++) printf("%g, ", u[i]);
    printf("\n");
    printf("dudx = ");
    for(int i=0; i<num; i++) printf("%g, ", dudx[i]);
    printf("\n");*/
    printf("v = ");
    for(int i=0; i<num; i++) printf("%g, ", v[i]);
    printf("\n");
    printf("dvdx = ");
    for(int i=0; i<num; i++) printf("%g, ", dvdx[i]);
    printf("\n");
    printf("u_prev = ");
    for(int i=0; i<num; i++) printf("%g, ", u_prev[i]);
    printf("\n");
    printf("du_prevdx = ");
    for(int i=0; i<num; i++) printf("%g, ", du_prevdx[i]);
    printf("\n");
    printf("f = ");
    for(int i=0; i<num; i++) printf("%g, ", f(x[i]));
    printf("\n");
    printf("val = %g\n", val);
  }
  return val;
};

double jacobian_surf_left(double x, double u, double dudx,
        double v, double dvdx, double u_prev, double du_prevdx,
        void *user_data)
{
  return (1/val_newton_alpha_left)*u*v;
}

double jacobian_surf_right(double x, double u, double dudx,
        double v, double dvdx, double u_prev, double du_prevdx,
        void *user_data)
{
  return (1/val_newton_alpha_right)*u*v;
}

double residual_surf_left(double x, double u_prev, double du_prevdx, double v,
        double dvdx, void *user_data)
{
  return -(val_newton_beta_left/val_newton_alpha_left) * v; 
}

double residual_surf_right(double x, double u_prev, double du_prevdx, double v,
        double dvdx, void *user_data)
{
  return -(val_newton_beta_right/val_newton_alpha_right) * v; 
}

/******************************************************************************/
int main() {
  // create mesh
  Mesh mesh(N_eq);
  mesh.create(A, B, N_elem);
  mesh.set_poly_orders(P_init);

  // boundary conditions
  mesh.set_bc_left_natural(0);
  mesh.set_bc_right_natural(0);
  mesh.assign_dofs();

  // register weak forms
  DiscreteProblem dp(&mesh);
  dp.add_matrix_form(0, 0, jacobian);
  dp.add_vector_form(0, residual_vol);
  dp.add_matrix_form_surf(0, 0, jacobian_surf_left, BOUNDARY_LEFT);
  dp.add_vector_form_surf(0, residual_surf_left, BOUNDARY_LEFT);
  dp.add_matrix_form_surf(0, 0, jacobian_surf_right, BOUNDARY_RIGHT);
  dp.add_vector_form_surf(0, residual_surf_right, BOUNDARY_RIGHT);

  // variable for the total number of DOF 
  int N_dof = mesh.get_n_dof();
  printf("N_dof = %d\n", N_dof);

  // allocate Jacobi matrix and residual
  Matrix *mat;
  double *y_prev = new double[N_dof];
  double *res = new double[N_dof];

  // zero initial condition for the Newton's method
  for(int i=0; i<N_dof; i++) y_prev[i] = 0; 

  // Newton's loop
  int newton_iterations = 0;
  while (1) {
    // zero the matrix:
    mat = new DenseMatrix(N_dof);

    // construct residual vector
    dp.assemble_matrix_and_vector(mat, res, y_prev); 

    if (DEBUG) {
        printf("RHS:");
        for(int i=0; i<N_dof; i++)
            printf("%f ", res[i]);
        printf("\n");
    }
  
    // calculate L2 norm of residual vector
    double res_norm = 0;
    for(int i=0; i<N_dof; i++) res_norm += res[i]*res[i];
    res_norm = sqrt(res_norm);

    // if residual norm less than TOL, quit
    // latest solution is in y_prev
    printf("Residual L2 norm: %.15f\n", res_norm);
    if (DEBUG)
        printf("TOL: %.15f\n", TOL);
    if(res_norm < TOL) break;

    // changing sign of vector res
    for(int i=0; i<N_dof; i++) res[i]*= -1;

    //mat->print();

    // solving the matrix system
    solve_linear_system(mat, res);

    // DEBUG: print solution
    if(DEBUG) {
      printf("New Y:\n");
      for(int i=0; i<N_dof; i++) {
        printf("%g ", res[i]);
      }
      printf("\n");
    }

    // updating y_prev by new solution which is in res
    for(int i=0; i<N_dof; i++) y_prev[i] += res[i];
    newton_iterations++;
    printf("Finished Newton iteration: %d\n", newton_iterations);
  }
  printf("Total number of Newton iterations: %d\n", newton_iterations);

  Linearizer l(&mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  printf("Done.\n");
  return 1;
}
