#include "discrete.h"

DiscreteProblem::DiscreteProblem(int neq, Mesh *mesh)
{
    this->neq = neq;
    this->mesh = mesh;
}

void DiscreteProblem::add_matrix_form(int i, int j, matrix_form fn)
{
    MatrixFormVol form = {i, j, fn};
    this->matrix_forms_vol.push_back(form);
}

void DiscreteProblem::add_vector_form(int i, vector_form fn)
{
    VectorFormVol form = {i, fn};
    this->vector_forms_vol.push_back(form);
}

void DiscreteProblem::add_matrix_form_surf(int i, int j, matrix_form_surf fn, int bdy_index)
{
    MatrixFormSurf form = {i, j, bdy_index, fn};
    this->matrix_forms_surf.push_back(form);
}

void DiscreteProblem::add_vector_form_surf(int i, vector_form_surf fn, int bdy_index)
{
    VectorFormSurf form = {i, bdy_index, fn};
    this->vector_forms_surf.push_back(form);
}

// c is solution component
void calculate_elem_coeffs(Mesh *mesh, int m, double *y_prev, double *coeffs, int c)
{ 
    Element *elems = mesh->get_elems();
    if (m == 0 && elems[m].dof[c][0] == -1) {
        coeffs[0] = mesh->bc_left_dir_values[c][0];
    }
    else {
        coeffs[0] = y_prev[elems[m].dof[c][0]];
    }
    if (m == mesh->get_n_elems()-1 && elems[m].dof[c][1] == -1) {
        coeffs[1] = mesh->bc_right_dir_values[c][0];
    }
    else {
        coeffs[1] = y_prev[elems[m].dof[c][1]];
    }
    for (int j=2; j<=elems[m].p; j++) {
        coeffs[j] = y_prev[elems[m].dof[c][j]];
    }
}

// process volumetric weak forms
// c is solution component
void DiscreteProblem::process_vol_forms(Matrix *mat, double *res, 
					double *y_prev, int matrix_flag, int c) {
  Element *elems = this->mesh->get_elems();
  for(int m=0; m < this->mesh->get_n_elems(); m++) {
    // to store quadrature data
    int    pts_num = 0;       // num of quad points
    double phys_pts[100];     // quad points
    double phys_weights[100]; // quad weights
    double phys_u[100];
    double phys_dudx[100];
    double phys_v[100];
    double phys_dvdx[100];
    double phys_u_prev[100];
    double phys_du_prevdx[100];
    // decide quadrature order and set up 
    // quadrature weights and points in element m
    int order = 2*elems[m].p; // FIXME - this needs to be improved.
    element_quadrature(elems[m].v1->x, elems[m].v2->x,  
	               order, phys_pts, phys_weights, &pts_num); 

    // evaluate previous solution and its derivative 
    // at all quadrature points in the element
    double coeffs[100];
    calculate_elem_coeffs(this->mesh, m, y_prev, coeffs, c); 
    double2 *ref_tab = g_quad_1d_std.get_points(order);
    double pts_array[pts_num];
    for (int j=0; j<pts_num; j++)
        pts_array[j] = ref_tab[j][0];
    element_solution(elems + m, coeffs, pts_num, 
                     pts_array, phys_u_prev, phys_du_prevdx, c); 

    // loop over test functions (rows)
    for(int i=0; i<elems[m].p + 1; i++) {
      // if i-th test function is active
      int pos_i = elems[m].dof[c][i]; // row index in matrix or residual vector
      if(pos_i != -1) {
        // transform i-th test function to element 'm'
        element_shapefn(elems[m].v1->x, elems[m].v2->x,  
                        i, order, phys_v, phys_dvdx); 
        // if we are constructing the matrix
        if(matrix_flag == 0 || matrix_flag == 1) {
          // loop over basis functions (columns)
          for(int j=0; j < elems[m].p + 1; j++) {
            int pos_j = elems[m].dof[c][j]; // matrix column index
            // if j-th basis function is active
	    if(pos_j != -1) {
              // transform j-th basis function to element 'm'
              element_shapefn(elems[m].v1->x, elems[m].v2->x,  
		 	      j, order, phys_u, phys_dudx); 
              // evaluate the bilinear form
              double val_ji = this->matrix_forms_vol[c].fn(pts_num, phys_pts,
                        phys_weights, phys_u, phys_dudx, phys_v, phys_dvdx,
                        phys_u_prev, phys_du_prevdx, NULL); 
              // add the result to the matrix
              mat->add(pos_j, pos_i, val_ji);
	    }
	  }
        }
        // contribute to residual vector
        if(matrix_flag == 0 || matrix_flag == 2) {
     	  double val_i = this->vector_forms_vol[c].fn(pts_num, phys_pts, phys_weights, 
                                  phys_u_prev, phys_du_prevdx, phys_v,
                                  phys_dvdx, NULL);
          // add the contribution to the residual vector
          if (DEBUG)
              printf("Adding to residual pos %d value %g\n", pos_i, val_i);
          res[pos_i] += val_i;
        }
      }
    }
  } 
}

// process left boundary weak forms
// c is solution component
void DiscreteProblem::process_surf_forms(Matrix *mat, double *res, 
					 double *y_prev, int matrix_flag, int bdy_index, int c) {
  Element *elems = this->mesh->get_elems();
  // evaluate previous solution and its derivative at the end point
  double phys_u_prev, phys_du_prevdx; // at the end point
  int m;
  if(bdy_index == BOUNDARY_LEFT) m = 0; // first element
  else m = this->mesh->get_n_elems()-1; // last element
  double coeffs[100];
  calculate_elem_coeffs(this->mesh, m, y_prev, coeffs, c); 
  double x_ref; 
  if(bdy_index == BOUNDARY_LEFT) x_ref = -1; // left end of reference element
  else x_ref = 1;                            // right end of reference element
  // getting solution value and derivative at the boundary point
  element_solution_point(x_ref, elems + m, coeffs,
                         &phys_u_prev, &phys_du_prevdx, c); 

  // surface integrals at the end point
  double phys_v, phys_dvdx; 
  double phys_u, phys_dudx;
  // loop over test functions on the boundary element
  for(int i=0; i<elems[m].p + 1; i++) {
    // if i-th test function is active
    int pos_i = elems[m].dof[c][i]; // row index in matrix or residual vector
    if(pos_i != -1) {
      // transform i-th test function to the boundary element
      element_shapefn_point(x_ref, elems[m].v1->x, elems[m].v2->x,  
                        i, &phys_v, &phys_dvdx); 
      // contribute to the matrix
      if(matrix_flag == 0 || matrix_flag == 1) {
        DiscreteProblem::MatrixFormSurf *matrix_form_surf=NULL;
        for (int ww = 0; ww < this->matrix_forms_surf.size(); ww++) {
          DiscreteProblem::MatrixFormSurf *vfs = &(this->matrix_forms_surf[ww]);
          if (vfs->bdy_index == bdy_index) {
            matrix_form_surf = vfs;
              break;
          }
	}
        if (matrix_form_surf != NULL) {
          // loop over basis functions on the boundary element
          for(int j=0; j < elems[m].p + 1; j++) {
            int pos_j = elems[m].dof[c][j]; // matrix column index
            // if j-th basis function is active
	    if(pos_j != -1) {
              // transform j-th basis function to the boundary element
              element_shapefn_point(x_ref, elems[m].v1->x, 
                                    elems[m].v2->x, j, &phys_u, 
                                    &phys_dudx); 
              // evaluate the surface bilinear form
              double val_ji_surf = matrix_form_surf->fn(elems[m].v1->x,
                      phys_u, phys_dudx, phys_v, 
                      phys_dvdx, phys_u_prev, phys_du_prevdx, 
                      NULL); 
                // add the result to the matrix
                mat->add(pos_j, pos_i, val_ji_surf);
	    }
	  }
	}
      }
      // contribute to residual vector
      if(matrix_flag == 0 || matrix_flag == 2) {
        DiscreteProblem::VectorFormSurf *vector_form_surf=NULL;
        for (int ww = 0; ww < this->vector_forms_surf.size(); ww++) {
          DiscreteProblem::VectorFormSurf *vfs = &(this->vector_forms_surf[ww]);
          if (vfs->bdy_index == bdy_index) {
            vector_form_surf = vfs;
            break;
          }
        }
        //if (vector_form_surf == NULL) error("Surface form not found.");
        if (vector_form_surf != NULL) {
       	  double val_i_surf = vector_form_surf->fn(elems[m].v1->x,  
                                  phys_u_prev, phys_du_prevdx, 
                                  phys_v, phys_dvdx, NULL);
          // add the contribution to the residual vector
          if (DEBUG)
              printf("Adding to residual pos %d value %g\n", pos_i, val_i_surf);
          res[pos_i] += val_i_surf;
        }
      }
    }
  }
}

// construct Jacobi matrix or residual vector
// matrix_flag == 0... assembling Jacobi matrix and residual vector together
// matrix_flag == 1... assembling Jacobi matrix only
// matrix_flag == 2... assembling residual vector only
// NOTE: Simultaneous assembling of the Jacobi matrix and residual
// vector is more efficient than if they are assembled separately
void DiscreteProblem::assemble(Matrix *mat, double *res, 
              double *y_prev, int matrix_flag) {
  // number of equations in the system
  int n_eq = this->mesh->n_eq;

  // total number of unknowns
  int n_dof = this->mesh->get_n_dof()

  // erase residual vector
  if(matrix_flag == 0 || matrix_flag == 2) 
    for(int i=0; i<n_dof; i++) res[i] = 0;

  // process volumetric weak forms via an element loop
  for(int c=0; c<n_eq; c++) process_vol_forms(mat, res, y_prev, matrix_flag, c);

  // process surface weak forms for the left boundary
  for(int c=0; c<n_eq; c++) 
    if(this->mesh->bc_left_dir[c] != 1) 
      process_surf_forms(mat, res, y_prev, matrix_flag, BOUNDARY_LEFT, c);

  // process surface weak forms for the right boundary
  for(int c=0; c<n_eq; c++) 
    if(this->mesh->bc_right_dir[0] != 1) process_surf_forms(mat, res, y_prev, matrix_flag, BOUNDARY_RIGHT, c);

  // DEBUG: print Jacobi matrix
  if(DEBUG && (matrix_flag == 0 || matrix_flag == 1)) {
    printf("Jacobi matrix:\n");
    for(int i=0; i<n_dof; i++) {
      for(int j=0; j<n_dof; j++) {
        printf("%g ", mat->get(i, j));
      }
    }
  }

  // DEBUG: print residual vector
  if(DEBUG && (matrix_flag == 0 || matrix_flag == 2)) {
    printf("Residual:\n");
    for(int i=0; i<n_dof; i++) {
      printf("%g ", res[i]);
    }
    printf("\n");
  }
} 

// construct both the Jacobi matrix and the residual vector
void DiscreteProblem::assemble_matrix_and_vector(Matrix *mat, double *res, double *y_prev) {
  assemble(mat, res, y_prev, 0);
} 

// construct Jacobi matrix only
void DiscreteProblem::assemble_matrix(Matrix *mat, double *y_prev) {
  double *void_res = NULL;
  assemble(mat, void_res, y_prev, 1);
} 

// construct residual vector only
void DiscreteProblem::assemble_vector(double *res, double *y_prev) {
  Matrix *void_mat = NULL;
  assemble(void_mat, res, y_prev, 2);
} 

// transformation of quadrature to physical element
void element_quadrature(double a, double b, 
                        int order, double *pts, double *weights, int *num) {
  double2 *ref_tab = g_quad_1d_std.get_points(order);
  *num = g_quad_1d_std.get_num_points(order);
  for (int i=0;i<*num;i++) {
    //change points and weights to interval (a, b)
    pts[i] = (b-a)/2.*ref_tab[i][0]+(b+a)/2.; 
    weights[i] = ref_tab[i][1]*(b-a)/2.;
  }
};

// transformation of k-th shape function defined on Gauss points 
// corresponding to 'order' to physical interval (a,b)
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
// in the "pts_array" points
void element_solution(Element *e, double **coeff, int pts_num, 
		      double *pts_array, double **val, double **der, int num_eq)
{
  double a = e->v1->x;
  double b = e->v2->x;
  double jac = (b-a)/2.; 
  int p = e->p;
  for (int i=0 ; i<pts_num; i++) {
    for(int c=0; i<num_eq; c++) der[c][i] = val[c][i] = 0;
    for(int j=0; j<=p; j++) {
      for(int c=0; i<num_eq; c++) val[c][i] += coeff[c][j]*lobatto_fn_tab_1d[j](pts_array[i]);
      for(int c=0; i<num_eq; c++) der[c][i] += coeff[c][j]*lobatto_der_tab_1d[j](pts_array[i]);
    }
    der[i] /= jac;
  }
} 

// transformation of k-th shape function at the reference 
// point x_ref to physical interval (a,b).
void element_shapefn_point(double x_ref, double a, double b, 
		     int k, double *val, double *der) {
    // change function values and derivatives to interval (a, b)
    *val = lobatto_fn_tab_1d[k](x_ref);
    double jac = (b-a)/2.; 
    *der = lobatto_der_tab_1d[k](x_ref) / jac; 
}

// evaluate previous solution and its derivative 
// at the reference point x_ref to element 'e'.
void element_solution_point(double x_ref, Element *e, 
			    double *coeff, double *val, double *der, int num_eq)
{
  double a = e->v1->x;
  double b = e->v2->x;
  double jac = (b-a)/2.; 
  int p = e->p;
  for(int c=0; i<num_eq; c++) *der = *val = 0;
  for(int j=0; j<=p; j++) {
    for(int c=0; i<num_eq; c++) *val += coeff[c][j]*lobatto_fn_tab_1d[j](x_ref);
    for(int c=0; i<num_eq; c++) *der += coeff[c][j]*lobatto_der_tab_1d[j](x_ref);
  }
  for(int c=0; i<num_eq; c++) *der /= jac;
} 


