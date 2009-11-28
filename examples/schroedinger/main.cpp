#include "hermes1d.h"

#include "_hermes1d_api.h"

static int N_eq = 1;
int N_elem = 100;                         // number of elements
double A = 0, B = 20;                // domain end points
int P_init = 2;                        // initial polynomal degree

double l = 0;

double lhs(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM], double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM], 
                void *user_data)
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
                double u_prev[MAX_EQN_NUM][MAX_PTS_NUM], double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += u[i]*v[i]*x[i]*x[i]*weights[i];
  }
  return val;
}

double E;

double residual(int num, double *x, double *weights, double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
        double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM], double *v, double *dvdx, void *user_data)
{
    double *u = &u_prev[0][0];
    double *dudx = &du_prevdx[0][0];
    double val = 0;
    for(int i = 0; i<num; i++) {
        double coeff;
        coeff = 0.5*x[i]*x[i]*dudx[i]*dvdx[i] -u[i]*v[i]*x[i]
            + 0.5 * (l + 1)*l *u[i]*v[i];
        coeff -= E*u[i]*v[i]*x[i]*x[i];
        val += coeff*weights[i];
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
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  // you can set the zero dirichlet at the right hand side
  //mesh.set_bc_right_dirichlet(0, 0);

  // variable for the total number of DOF 
  int N_dof = mesh->assign_dofs();
  printf("ndofs: %d", N_dof);

  // register weak forms
  DiscreteProblem *dp1 = new DiscreteProblem();
  dp1->add_matrix_form(0, 0, lhs);
  DiscreteProblem *dp2 = new DiscreteProblem();
  dp2->add_matrix_form(0, 0, rhs);

  DiscreteProblem *dp3 = new DiscreteProblem();
  dp3->add_vector_form(0, residual);

  // allocate Jacobi matrix and residual
  DenseMatrix *mat1 = new DenseMatrix(N_dof);
  DenseMatrix *mat2 = new DenseMatrix(N_dof);
  double *y_prev = new double[N_dof];

  // zero initial condition for the Newton's method
  for(int i=0; i<N_dof; i++) y_prev[i] = 0; 

  dp1->assemble_matrix(mesh, mat1, y_prev);
  dp2->assemble_matrix(mesh, mat2, y_prev);

  printf("Importing hermes1d\n");
  // Initialize Python
  Py_Initialize();
  PySys_SetArgv(argc, argv);
  if (import_hermes1d___hermes1d())
      throw std::runtime_error("hermes1d failed to import.");

  cmd("print 'Python initialized'");
  insert_matrix(mat1, N_dof); cmd("A = _");
  insert_matrix(mat2, N_dof); cmd("B = _");
  cmd("from utils import solve");
  cmd("eigs = solve(A, B)");
  cmd("E, v = eigs[0]");
  //cmd("print 'v[0]', v[0]");
  //cmd("print E");
  double *v;
  int n;
  numpy2c_double_inplace(get_object("v"), &v, &n);

  double *res = new double[N_dof];
  E = py2c_double(get_object("E"));
  printf("E=%.10f\n", E);
  E = -0.5;
  dp3->assemble_vector(mesh, res, v);
  // calculate L2 norm of residual vector
  double res_norm = 0;
  for(int i=0; i<N_dof; i++) res_norm += res[i]*res[i];
  res_norm = sqrt(res_norm);
  printf("L2 norm of the residual: %f\n", res_norm);


  Linearizer l(mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, v);

  printf("still ok\n");
  insert_object("mesh", c2py_mesh(&mesh));
  printf("2\n");
  cmd("from plot import plot_eigs, plot_file");
  cmd("plot_eigs(mesh, eigs)");
  //cmd("plot_eigs(eigs)");
  printf("Done.\n");
  return 0;
}
