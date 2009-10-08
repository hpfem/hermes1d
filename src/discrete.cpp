#include "discrete.h"

DiscreteProblem::DiscreteProblem(Mesh *mesh)
{
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

// process volumetric weak forms
// c is solution component
void DiscreteProblem::process_vol_forms(Matrix *mat, double *res, 
					double *y_prev, int matrix_flag, int c) {
  int n_eq = this->mesh->get_n_eq();
  Element *elems = this->mesh->get_elems();
  int n_elem = this->mesh->get_n_elems();
  for(int m=0; m < n_elem; m++) {
    // variables to store quadrature data
    // FIXME: now maximum number of Gauss points is 100
    int    pts_num = 0;       // num of quad points
    double phys_pts[100];     // quad points
    double phys_weights[100]; // quad weights
    double phys_u[100];       // basis function 
    double phys_dudx[100];    // basis function x-derivative
    double phys_v[100];       // test function
    double phys_dvdx[100];    // test function x-derivative
    // FIXME: now maximum limit of equations is 10, 
    // and number of Gauss points is limited to 100
    if(n_eq > 10) error("number of equations too high in process_vol_forms().");
    double phys_u_prev[10][100];     // previous solution, all components
    double phys_du_prevdx[10][100];  // previous solution x-derivative, all components
    // decide quadrature order and set up 
    // quadrature weights and points in element m
    // FIXME: for some equations this may not be enough!
    int order = 2*elems[m].p;

    // prepare quadrature points and weights in physical element m
    create_element_quadrature(elems[m].v1->x, elems[m].v2->x,  
	               order, phys_pts, phys_weights, &pts_num); 

    // prepare quadrature points in reference interval (-1, 1)
    double ref_pts_array[100];
    double2 *ref_tab = g_quad_1d_std.get_points(order);
    for (int j=0; j<pts_num; j++) ref_pts_array[j] = ref_tab[j][0];

    // evaluate previous solution and its derivative 
    // at all quadrature points in the element, 
    // for every solution component
    double coeffs[10][100];
    this->mesh->calculate_elem_coeffs(m, y_prev, coeffs); 
    this->mesh->element_solution(elems + m, coeffs, pts_num, 
                     ref_pts_array, phys_u_prev, phys_du_prevdx); 

    // loop over test functions (rows)
    for(int i=0; i<elems[m].p + 1; i++) {
      // if i-th test function is active
      int pos_i = elems[m].dof[c][i]; // row index in matrix or residual vector
      if(pos_i != -1) {
        // transform i-th test function to element 'm'
        this->mesh->element_shapefn(elems[m].v1->x, elems[m].v2->x,  
                        i, order, phys_v, phys_dvdx); 
        // if we are constructing the matrix
        if(matrix_flag == 0 || matrix_flag == 1) {
          // loop over basis functions (columns)
          for(int j=0; j < elems[m].p + 1; j++) {
            int pos_j = elems[m].dof[c][j]; // matrix column index
            // if j-th basis function is active
	    if(pos_j != -1) {
              // transform j-th basis function to element 'm'
              this->mesh->element_shapefn(elems[m].v1->x, elems[m].v2->x,  
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
  // FIXME: maximum number of equations limited by 10
  double phys_u_prev[10], phys_du_prevdx[10]; // at the end point
  int m;
  if(bdy_index == BOUNDARY_LEFT) m = 0; // first element
  else m = this->mesh->get_n_elems()-1; // last element
  double coeffs[10][100];
  this->mesh->calculate_elem_coeffs(m, y_prev, coeffs); 
  double x_ref; 
  if(bdy_index == BOUNDARY_LEFT) x_ref = -1; // left end of reference element
  else x_ref = 1;                            // right end of reference element
  // getting solution value and derivative at the boundary point
  this->mesh->element_solution_point(x_ref, elems + m, coeffs,
                         phys_u_prev, phys_du_prevdx); 

  // surface integrals at the end point
  double phys_v, phys_dvdx; 
  double phys_u, phys_dudx;
  // loop over test functions on the boundary element
  for(int i=0; i<elems[m].p + 1; i++) {
    // if i-th test function is active
    int pos_i = elems[m].dof[c][i]; // row index in matrix or residual vector
    if(pos_i != -1) {
      // transform i-th test function to the boundary element
      this->mesh->element_shapefn_point(x_ref, elems[m].v1->x, elems[m].v2->x,  
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
              this->mesh->element_shapefn_point(x_ref, elems[m].v1->x, 
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
  int n_eq = this->mesh->get_n_eq();

  // total number of unknowns
  int n_dof = this->mesh->get_n_dof();

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






