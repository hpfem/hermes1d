#include "discrete.h"

DiscreteProblem::DiscreteProblem(int neq, Mesh *mesh)
{
    this->neq = neq;
    this->mesh = mesh;
}

void DiscreteProblem::add_matrix_form(int i, int j, matrix_form fn)
{
    this->_matrix_form = fn;
}

void DiscreteProblem::add_vector_form(int i, int j, vector_form fn)
{
    this->_vector_form = fn;
}

// construct Jacobi matrix or residual vector
// matrix_flag == 0... assembling Jacobi matrix and residual vector together
// matrix_flag == 1... assembling Jacobi matrix only
// matrix_flag == 2... assembling residual vector only
// NOTE: Simultaneous assembling of the Jacobi matrix and residual
// vector is more efficient than if they are assembled separately
void DiscreteProblem::assemble(int ndof, Element *elems, double **mat, double *res, 
              double *y_prev, int matrix_flag) {

  int DEBUG = 1;
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
  for(int m=0; m < this->mesh->get_n_elems(); m++) {
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
              double val_ji = this->_matrix_form(pts_num, phys_pts, phys_weights, 
				       phys_u, phys_dudx, phys_v, phys_dvdx, 
                                       phys_u_prev, phys_du_prevdx); 
              // add the result to the matrix
              mat[pos_j][pos_i] += val_ji; 
	    }
	  }
        }
        // contribute to residual vector
        if(matrix_flag == 0 || matrix_flag == 2) {
     	  double val_i = this->_vector_form(pts_num, phys_pts, phys_weights, 
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
void DiscreteProblem::assemble_matrix_and_vector(int ndof, Element *elems, double **mat, double *res, double *y_prev) {
  assemble(ndof, elems, mat, res, y_prev, 0);
  return;
} 

// construct Jacobi matrix only
void DiscreteProblem::assemble_matrix(int ndof, Element *elems, double **mat, double *y_prev) {
  double *void_res = NULL;
  assemble(ndof, elems, mat, void_res, y_prev, 1);
  return;
} 

// construct residual vector only
void DiscreteProblem::assemble_vector(int ndof, Element *elems, double *res, double *y_prev) {
  double **void_mat = NULL;
  assemble(ndof, elems, void_mat, res, y_prev, 2);
  return;
} 

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
