#include "hermes1d.h"

static int NUM_EQ = 1;
int Nelem = 6;                         // number of elements
double A = 0, B = 2*M_PI;                // domain end points
int P_INIT = 1;                        // initial polynomal degree

double l = 1;
// bilinear form for the Jacobi matrix 
// pts_num...number of Gauss points in element
// pts[]...Gauss points
// weights[]...Gauss weights for points in pts[]
// u...basis function
// v...test function
// u_prev...previous solution
double lhs(int pts_num, double *pts, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double *u_prev, double *du_prevdx)
{
    double val = 0;
    for(int i = 0; i<pts_num; i++) {
        double coeff;
        coeff = 0.5*pts[i]*pts[i]*dudx[i]*dvdx[i] -u[i]*v[i]*pts[i]
            + 0.5 * (l + 1)*l *u[i]*v[i];
        val += coeff*weights[i];
    }
    return val;
}

double rhs(int pts_num, double *pts, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double *u_prev, double *du_prevdx)
{
  double val = 0;
  for(int i = 0; i<pts_num; i++) {
    val += u[i]*v[i]*pts[i]*pts[i]*weights[i];
  }
  return val;
}

/******************************************************************************/
int main() {
  // create mesh
  Mesh mesh(NUM_EQ);
  mesh.create(A, B, Nelem);
  mesh.set_poly_orders(P_INIT);
  mesh.set_dirichlet_bc_left(0, 0);
  mesh.assign_dofs();

  // register weak forms
  DiscreteProblem dp1(NUM_EQ, &mesh);
  dp1.add_matrix_form(0, 0, lhs);
  DiscreteProblem dp2(NUM_EQ, &mesh);
  dp2.add_matrix_form(0, 0, rhs);

  // variable for the total number of DOF 
  int Ndof = mesh.get_n_dof();

  // allocate Jacobi matrix and residual
  double **mat1 = new_matrix<double>(Ndof, Ndof);
  double **mat2 = new_matrix<double>(Ndof, Ndof);
  double *y_prev = new double[Ndof];

  // zero initial condition for the Newton's method
  for(int i=0; i<Ndof; i++) y_prev[i] = 0; 

  dp1.assemble_matrix(mat1, y_prev);
  dp2.assemble_matrix(mat2, y_prev);

  /*
  for(int i=0; i<Ndof; i++) {
      for(int j=0; j<Ndof; j++)
          printf("%f ", mat2[i][j]);
      printf("\n");
  }
  */


  /*
  Linearizer l(&mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  printf("Output written to %s.\n", out_filename);
  printf("Done.\n");
  */
  return 0;
}
