// including higher-order shape functions
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

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

// Dirichlet boundary conditions at both endpoints
// (first integer in pair indicates whether there is 
// a Dirichlet condition for that equation, the other
// one tells the value)
int2 DIR_BC_LEFT[] = { {1, 1} };
int2 DIR_BC_RIGHT[] = { {1, 2} };

// right-hand side
double f(double x) {
  return sin(x);
  //return 1;
}

int plotting_elem_subdivision = 100; 

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

// transformation of quadrature to physical element
void element_quadrature(double a, double b, 
                        int order, double *pts, double *weights, int &num) {
  double2 *ref_tab = g_quad_1d_std.get_points(order);
  num = g_quad_1d_std.get_num_points(order);
  for (int i=0;i<num;i++) {
    //change points and weights to interval (a, b)
    pts[i] = (b-a)/2.*ref_tab[i][0]+(b+a)/2.; 
    weights[i] = ref_tab[i][1]*(b-a)/2.;
  }
};

// transformation of k-th shape function defined on Gauss points 
// corresponding to 'order' to physical element.
void element_shapefn(double a, double b, 
		     int k, int order, double *val, double *der) {
  double2 *ref_tab = g_quad_1d_std.get_points(order);
  int pts_num = g_quad_1d_std.get_num_points(order);
  for (int i=0 ; i<pts_num; i++) {
    // change function values and derivatives to interval (a, b)
    val[i] = lobatto_fn_tab_1d[k](ref_tab[i][0]);
    double jac = (b-a)/2.; 
    der[i] = lobatto_der_tab_1d[k](ref_tab[i][0]) / jac; 
  }
};

// evaluate previous solution and its derivative 
// in the Gauss quadrature points of order 'order'
void element_solution(Element *e, double *y_prev, int order, 
                      double *val, double *der) {
  double a = e->v1->x;
  double b = e->v2->x;
  double jac = (b-a)/2.; 
  int p = e->p;
  double2 *ref_tab = g_quad_1d_std.get_points(order);
  int pts_num = g_quad_1d_std.get_num_points(order);
  for (int i=0 ; i<pts_num; i++) {
    der[i] = val[i] = 0;
    for(int j=0; j<=p; j++) {
      val[i] += y_prev[e->dof[j]]*lobatto_fn_tab_1d[j](ref_tab[i][0]);
      der[i] += y_prev[e->dof[j]]*lobatto_der_tab_1d[j](ref_tab[i][0]);
    }
    der[i] /= jac;
  }
  return;
} 

// evaluate approximate solution at element 'm' at reference 
// point 'x_ref'. Here 'y' is the global vector of coefficients
void eval_approx(Element *e, double x_ref, double *y, double &x_phys, double &val) {
  val = 0;
  for(int i=0; i <= e->p; i++) {
    if(e->dof[i] >= 0) val += y[e->dof[i]]*lobatto_fn_tab_1d[i](x_ref);
  }
  double a = e->v1->x;
  double b = e->v2->x;
  x_phys = (a+b)/2 + x_ref*(b-a)/2; 
  return;
} 

// construct Jacobi matrix or residual vector
// matrix_flag == 0... assembling Jacobi matrix and residual vector together
// matrix_flag == 1... assembling Jacobi matrix only
// matrix_flag == 2... assembling residual vector only
// NOTE: Simultaneous assembling of the Jacobi matrix and residual
// vector is more efficient than if they are assembled separately
void assemble(int ndof, Element *elems, double **mat, double *res, 
              double *y_prev, int matrix_flag) {
  // erase matrix
  if(matrix_flag == 0 || matrix_flag == 1) {
    for(int i=0; i<ndof; i++) {
      for(int j=0; j<ndof; j++) mat[i][j] = 0;
    }
  }

  // erase residual vector
  if(matrix_flag == 0 || matrix_flag == 2) 
    for(int i=0; i<ndof; i++) res[i] = 0;

  // memory for quadrature data
  // FIXME - this is a memory leak
  int     pts_num = 0;                    // num of quad points
  double *phys_pts =  new double[100];    // quad points
  double *phys_weights = new double[100]; // quad weights
  double *phys_u =    new double[100];
  double *phys_dudx = new double[100];
  double *phys_v =    new double[100];
  double *phys_dvdx = new double[100];
  double *phys_u_prev =    new double[100];
  double *phys_du_prevdx = new double[100];

  // element loop
  for(int m=0; m<Nelem; m++) {
    // decide quadrature order and set up 
    // quadrature weights and points in element m
    int order = 2*elems[m].p; // Now the order is the same for the matrix 
                              // and residual. FIXME - this needs to be improved.
    element_quadrature(elems[m].v1->x, elems[m].v2->x,  
	               order, phys_pts, phys_weights, pts_num); 
    if(DEBUG) {
      printf("--------------------------------------------\n");
      printf("Element %d     interval (%g, %g)\n", m, elems[m].v1->x, elems[m].v2->x);
      printf("p = %d      order = %d     pts_num = %d\n", elems[m].p, order, pts_num);
      printf("Pts and weights: \n");
      for(int s=0; s < pts_num; s++) printf("[%g, %g]\n", 
             phys_pts[s], phys_weights[s]);
    }
    // evaluate previous solution and its derivative in the quadrature points
    element_solution(elems + m, y_prev, order, phys_u_prev, phys_du_prevdx); 
    // loop over test functions (rows)
    for(int i=0; i<elems[m].p + 1; i++) {
      // if i-th test function is active
      int pos_i = elems[m].dof[i]; // row index in matrix or residual vector
      if(pos_i != -1) {
        // transform i-th test function to element 'm'
        element_shapefn(elems[m].v1->x, elems[m].v2->x,  
                        i, order, phys_v, phys_dvdx); 
        // if we are constructing the matrix
        if(matrix_flag == 0 || matrix_flag == 1) {
          // loop over basis functions (columns)
          for(int j=0; j < elems[m].p + 1; j++) {
            int pos_j = elems[m].dof[j]; // matrix column index
            // if j-th basis function is active
	    if(pos_j != -1) {
              // transform j-th basis function to element 'm'
              element_shapefn(elems[m].v1->x, elems[m].v2->x,  
			      j, order, phys_u, phys_dudx); 
              // evaluate the bilinear form
              double val_ji = jacobian(pts_num, phys_pts, phys_weights, 
				       phys_u, phys_dudx, phys_v, phys_dvdx, 
                                       phys_u_prev, phys_du_prevdx); 
              // add the result to the matrix
              mat[pos_j][pos_i] += val_ji; 
	    }
	  }
        }
        // contribute to residual vector
        if(matrix_flag == 0 || matrix_flag == 2) {
     	  double val_i = residual(pts_num, phys_pts, phys_weights, 
                                  phys_u_prev, phys_du_prevdx, phys_v,
                                  phys_dvdx);
          // add the contribution to the residual vector
          printf("Adding to residual pos %d value %g\n", pos_i, val_i);
          res[pos_i] += val_i;
        }
      }
    }
  }

  // DEBUG: print Jacobi matrix
  if(DEBUG && (matrix_flag == 0 || matrix_flag == 1)) {
    printf("Jacobi matrix:\n");
    for(int i=0; i<ndof; i++) {
      for(int j=0; j<ndof; j++) {
        printf("%g ", mat[i][j]);
      }
      printf("\n");
    }
  }

  // DEBUG: print residual vector
  if(DEBUG && (matrix_flag == 0 || matrix_flag == 2)) {
    printf("Residual:\n");
    for(int i=0; i<ndof; i++) {
      printf("%g ", res[i]);
    }
    printf("\n");
  }

  return;
} 

// construct both the Jacobi matrix and the residual vector
void assemble_matrix_and_vector(int ndof, Element *elems, double **mat, double *res, double *y_prev) {
  assemble(ndof, elems, mat, res, y_prev, 0);
  return;
} 

// construct Jacobi matrix only
void assemble_matrix(int ndof, Element *elems, double **mat, double *y_prev) {
  double *void_res = NULL;
  assemble(ndof, elems, mat, void_res, y_prev, 1);
  return;
} 

// construct residual vector only
void assemble_vector(int ndof, Element *elems, double *res, double *y_prev) {
  double **void_mat = NULL;
  assemble(ndof, elems, void_mat, res, y_prev, 2);
  return;
} 

void intro() {
  printf("\n-------------------------------------------\n");
  printf("   This is Hermes1D - a free ODE solver\n");
  printf(" based on the hp-FEM and Newton's method,\n");
  printf("   developed by the hp-FEM group at UNR\n");
  printf("  and distributed under the GPL license.\n");
  printf(" For more details visit http://hpfem.org/.\n");
  printf("-------------------------------------------\n");
}

/******************************************************************************/
int main() {
  // introductory text
  intro();

  // variable for the total number of DOF 
  int Ndof;				       

  Mesh mesh;
  mesh.create(A, B, Nelem);
  mesh.set_poly_orders(P_INIT);

  // create mesh
  Vertex *Vertices = new Vertex[Nelem+1];    // allocate array of vertices
  double h = (B - A)/Nelem;
  for(int i = 0; i < Nelem+1; i++) {
    Vertices[i].x = A + i*h;             // so far equidistant division only
  }
  Element *Elems = new Element[Nelem];                // allocate array of elements
  for(int i=0; i<Nelem; i++) {
    Elems[i].p = P_INIT; 
    Elems[i].v1 = Vertices + i;
    Elems[i].v2 = Vertices + i + 1;
    Elems[i].dof = new int[P_INIT+1];
  }

  // define element connectivities
  // (so far only for zero Dirichlet conditions)
  // (a) enumerate vertex dofs
  Elems[0].dof[0] = -1;        // Dirichlet left end of domain
  Elems[0].dof[1] = 0;         // first vertex dof
  for(int i=1; i<Nelem-1; i++) {
    Elems[i].dof[0] = i-1;     // Dirichlet left end of domain
    Elems[i].dof[1] = i;       // first dof
  }
  Elems[Nelem-1].dof[0] = Nelem-2;     // last vertex dof
  Elems[Nelem-1].dof[1] = -1;      // Dirichlet left end of domain
  // (b) enumerate bubble dofs
  Ndof = Nelem-1; 
  for(int i=0; i<Nelem; i++) {
    for(int j=2; j<=P_INIT; j++) {
      Elems[i].dof[j] = Ndof++;     // enumerating higher-order dofs
    }
  }

  // test (print DOF in elements)
  int n_test = 0;
  for(int i=0; i<Nelem; i++) n_test += Elems[i].p;
  n_test -= 1;
  if(Ndof != n_test) {
    printf("Ndof = %d, n_test = %d\n", Ndof, n_test);
    error("Internal: Test of total DOF number failed."); 
  }

  // test (print element connectivities)
  if(DEBUG) {
    printf("Printing element DOF arrays:\n");
    printf("Elements = %d\n", Nelem);
    printf("DOF = %d", Ndof);
    for (int i = 0; i < Nelem; i++) {
      printf("\nElement[%d]: ", i); 
      for(int j = 0; j<Elems[i].p+1; j++) {
        printf("%d, ", Elems[i].dof[j]);
      }
    }
    printf("\n"); 
  }

  DiscreteProblem dp(1);
  dp.add_matrix_form(0, 0, jacobian);
  dp.add_vector_form(0, 0, residual);

  // allocate Jacobi matrix and residual
  double **mat = new_matrix<double>(Ndof,Ndof);
  double *y_prev = new double[Ndof];
  double *res = new double[Ndof];

  // zero initial condition for the Newton's method
  for(int i=0; i<Ndof; i++) y_prev[i] = 0; 

  // Newton's loop
  while (1) {
    // construct residual vector
    assemble_matrix_and_vector(Ndof, Elems, mat, res, y_prev); 

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

  // Plot solution in Gnuplot format
  char out_filename[100];
  strcpy(out_filename, "solution.gp");
  FILE *f = fopen(out_filename, "wb"); 
  for(int m=0; m<Nelem; m++) {
    double h = 2./plotting_elem_subdivision; 
    double x_phys, val;
    for(int i=0; i<=plotting_elem_subdivision; i++) {
      double x_ref = -1. + i*h;  
      eval_approx(Elems + m, x_ref, y_prev, x_phys, val);
      fprintf(f, "%g %g\n", x_phys, val);
    }
  }
  fclose(f);  

  printf("Output written to %s.\n", out_filename);
  printf("Done.\n");
  return 1;
  };
