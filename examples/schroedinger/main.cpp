#include "hermes1d.h"

#include "_hermes1d_api.h"

static int NUM_EQ = 1;
int Nelem = 100;                         // number of elements
double A = 0, B = 30;                // domain end points
int P_INIT = 2;                        // initial polynomal degree

double l = 1;

double lhs(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[10][100], double du_prevdx[10][100], void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        double coeff;
        coeff = 0.5*x[i]*x[i]*dudx[i]*dvdx[i] -u[i]*v[i]*x[i]
            + 0.5 * (l + 1)*l *u[i]*v[i];
        val += coeff*weights[i];
    }
    return val;
}

double rhs(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[10][100], double du_prevdx[10][100], void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += u[i]*v[i]*x[i]*x[i]*weights[i];
  }
  return val;
}

void insert_matrix(DenseMatrix *mat, int len)
{
  double _mat[len*len];
  for(int i=0; i<len; i++)
      for(int j=0; j<len; j++)
          _mat[i*len+j] = mat->get(i, j);
  insert_object("len", c2py_int(len));
  insert_object("mat", c2numpy_double(_mat, len*len));
  cmd("_ = mat.reshape((len, len))");
  cmd("del len");
  cmd("del mat");
}

/******************************************************************************/
int main(int argc, char* argv[]) {
  // create mesh
  Mesh mesh(NUM_EQ);
  mesh.create(A, B, Nelem);
  mesh.set_uniform_poly_order(P_INIT);
  mesh.set_bc_left_dirichlet(0, 0);
  mesh.assign_dofs();

  // register weak forms
  DiscreteProblem dp1(&mesh);
  dp1.add_matrix_form(0, 0, lhs);
  DiscreteProblem dp2(&mesh);
  dp2.add_matrix_form(0, 0, rhs);

  // variable for the total number of DOF 
  int Ndof = mesh.get_n_dof();

  // allocate Jacobi matrix and residual
  DenseMatrix *mat1 = new DenseMatrix(Ndof);
  DenseMatrix *mat2 = new DenseMatrix(Ndof);
  double *y_prev = new double[Ndof];

  // zero initial condition for the Newton's method
  for(int i=0; i<Ndof; i++) y_prev[i] = 0; 

  dp1.assemble_matrix(mat1, y_prev);
  dp2.assemble_matrix(mat2, y_prev);

  printf("Importing hermes1d\n");
  // Initialize Python
  Py_Initialize();
  PySys_SetArgv(argc, argv);
  if (import_hermes1d___hermes1d())
      throw std::runtime_error("hermes1d failed to import.");

  cmd("print 'Python initialized'");
  insert_matrix(mat1, Ndof); cmd("A = _");
  insert_matrix(mat2, Ndof); cmd("B = _");
  cmd("from utils import solve");
  cmd("v = solve(A, B)");
  double *v;
  int n;
  numpy2c_double_inplace(get_object("v"), &v, &n);

  Linearizer l(&mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, v);

  printf("Output written to %s.\n", out_filename);
  cmd("import plot");
  printf("Done.\n");
  return 0;
}
