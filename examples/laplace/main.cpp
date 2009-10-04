#include "hermes1d.h"

int DEBUG = 1;

// ********************************************************************

// general input:
static int NUM_EQ = 1;
int Nelem = 2;                         // number of elements
double A = 0, B = 2*M_PI;                // domain end points
int P_INIT = 3;                        // initial polynomal degree

// Tolerance for Newton's method
double TOL = 1e-8;

// right-hand side
double f(double x) {
  return sin(x);
  //return 1;
}

// ********************************************************************

// bilinear form for the Jacobi matrix 
// pts_num...number of Gauss points in element
// pts[]...Gauss points
// weights[]...Gauss weights for points in pts[]
// u...basis function
// v...test function
// u_prev...previous solution
double jacobian(int pts_num, double *pts, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double *u_prev, double *du_prevdx)
{
  double val = 0;
  for(int i = 0; i<pts_num; i++) {
    val += dudx[i]*dvdx[i]*weights[i];
  }
  return val;
};

// (nonlinear) form for the residual vector
// pts_num...number of Gauss points in element
// pts[]...Gauss points
// weights[]...Gauss weights for points in pts[]
// u...approximate solution
// v...test function
// u_prev...previous solution
double residual(int pts_num, double *pts, double *weights, 
                double *u_prev, double *du_prevdx, double *v, double *dvdx)
{
  double val = 0;
  for(int i = 0; i<pts_num; i++) {
    val += (du_prevdx[i]*dvdx[i] - f(pts[i])*v[i])*weights[i];
  }
  if(DEBUG) {
    /*printf("u = ");
    for(int i=0; i<pts_num; i++) printf("%g, ", u[i]);
    printf("\n");
    printf("dudx = ");
    for(int i=0; i<pts_num; i++) printf("%g, ", dudx[i]);
    printf("\n");*/
    printf("v = ");
    for(int i=0; i<pts_num; i++) printf("%g, ", v[i]);
    printf("\n");
    printf("dvdx = ");
    for(int i=0; i<pts_num; i++) printf("%g, ", dvdx[i]);
    printf("\n");
    printf("u_prev = ");
    for(int i=0; i<pts_num; i++) printf("%g, ", u_prev[i]);
    printf("\n");
    printf("du_prevdx = ");
    for(int i=0; i<pts_num; i++) printf("%g, ", du_prevdx[i]);
    printf("\n");
    printf("f = ");
    for(int i=0; i<pts_num; i++) printf("%g, ", f(pts[i]));
    printf("\n");
    printf("val = %g\n", val);
  }
  return val;
};

/******************************************************************************/
int main() {
  // create mesh
  Mesh mesh(NUM_EQ);
  mesh.create(A, B, Nelem);
  mesh.set_poly_orders(P_INIT);
  mesh.set_dirichlet_bc_left(0, 1);
  mesh.set_dirichlet_bc_right(0, 5);
  mesh.assign_dofs();

  // register weak forms
  DiscreteProblem dp(NUM_EQ, &mesh);
  dp.add_matrix_form(0, 0, jacobian);
  dp.add_vector_form(0, residual);

  // variable for the total number of DOF 
  int Ndof = mesh.get_n_dof();

  // allocate Jacobi matrix and residual
  double **mat = new_matrix<double>(Ndof,Ndof);
  double *y_prev = new double[Ndof];
  double *res = new double[Ndof];

  // zero initial condition for the Newton's method
  for(int i=0; i<Ndof; i++) y_prev[i] = 0; 

  // Newton's loop
  while (1) {
    // construct residual vector
    dp.assemble_matrix_and_vector(mat, res, y_prev); 

    // calculate L2 norm of residual vector
    double res_norm = 0;
    for(int i=0; i<Ndof; i++) res_norm += res[i]*res[i];
    res_norm = sqrt(res_norm);

    // if residual norm less than TOL, quit
    // latest solution is in y_prev
    if(res_norm < TOL) break;

    // changing sign of vector res
    for(int i=0; i<Ndof; i++) res[i]*= -1;

    // solve linear system
    int *indx = new int[Ndof];
    double d;
    ludcmp(mat, Ndof, indx, &d);
    lubksb(mat, Ndof, indx, res);

    // DEBUG: print solution
    if(DEBUG) {
      printf("New Y:\n");
      for(int i=0; i<Ndof; i++) {
        printf("%g ", res[i]);
      }
      printf("\n");
    }

    // updating y_prev by new solution which is in res
    for(int i=0; i<Ndof; i++) y_prev[i] += res[i];
  }

  Linearizer l(&mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  printf("Output written to %s.\n", out_filename);
  printf("Done.\n");
  return 1;
}
