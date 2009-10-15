// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

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
void DiscreteProblem::process_vol_forms(Matrix *mat, double *res, 
					double *y_prev, int matrix_flag) {
  int n_eq = this->mesh->get_n_eq();
  Element *elems = this->mesh->get_elems();
  int n_elem = this->mesh->get_n_elems();
  for(int m=0; m < n_elem; m++) {
    //printf("Processing elem %d\n", m);
    // variables to store quadrature data
    // FIXME: now maximum number of Gauss points is [MAX_EQN_NUM][MAX_PTS_NUM]0
    int    pts_num = 0;       // num of quad points
    double phys_pts[MAX_PTS_NUM];                  // quad points
    double phys_weights[MAX_PTS_NUM];              // quad weights
    double phys_u[MAX_PTS_NUM];                    // basis function 
    double phys_dudx[MAX_PTS_NUM];                 // basis function x-derivative
    double phys_v[MAX_PTS_NUM];                    // test function
    double phys_dvdx[MAX_PTS_NUM];                 // test function x-derivative
    // FIXME: now maximum limit of equations is [MAX_EQN_NUM][MAX_PTS_NUM], 
    // and number of Gauss points is limited to [MAX_EQN_NUM][MAX_PTS_NUM]0
    if(n_eq > MAX_EQN_NUM) error("number of equations too high in process_vol_forms().");
    double phys_u_prev[MAX_EQN_NUM][MAX_PTS_NUM];     // previous solution, all components
    double phys_du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM];  // previous solution x-derivative, all components
    // decide quadrature order and set up 
    // quadrature weights and points in element m
    // FIXME: for some equations this may not be enough!
    int order = 20; // 2*elems[m].p;

    // prepare quadrature points and weights in physical element m
    create_element_quadrature(elems[m].x1, elems[m].x2,  
	               order, phys_pts, phys_weights, &pts_num); 

    // prepare quadrature points in reference interval (-1, 1)
    double ref_pts_array[MAX_PTS_NUM];
    double2 *ref_tab = g_quad_1d_std.get_points(order);
    for (int j=0; j<pts_num; j++) ref_pts_array[j] = ref_tab[j][0];

    // evaluate previous solution and its derivative 
    // at all quadrature points in the element, 
    // for every solution component
    double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
    this->mesh->calculate_elem_coeffs(m, y_prev, coeffs); 
    this->mesh->element_solution(elems + m, coeffs, pts_num, 
                     ref_pts_array, phys_u_prev, phys_du_prevdx); 

    // DEBUG - SETTING THE PREVIOUS SOLUTION 
    // TO BE THE EXACT SOLUTION  
    /*
    for (int i=0 ; i<pts_num; i++) {
      phys_u_prev[0][i] = exp(phys_pts[i]);
      phys_du_prevdx[0][i] = exp(phys_pts[i]);
      phys_u_prev[1][i] = exp(-phys_pts[i]);
      phys_du_prevdx[1][i] = -exp(-phys_pts[i]);
    }
    */    

    // volumetric bilinear forms
    if(matrix_flag == 0 || matrix_flag == 1) 
    {
      for (int ww = 0; ww < this->matrix_forms_vol.size(); ww++)
      {
	MatrixFormVol *mfv = &this->matrix_forms_vol[ww];
	int c_i = mfv->i;  
	int c_j = mfv->j;  

	// loop over test functions (rows)
	for(int i=0; i<elems[m].p + 1; i++) {
	  // if i-th test function is active
	  int pos_i = elems[m].dof[c_i][i]; // row in matrix
	  if(pos_i != -1) {
	    // transform i-th test function to element 'm'
	    this->mesh->element_shapefn(elems[m].x1, elems[m].x2,  
			    i, order, phys_v, phys_dvdx); 
	    // if we are constructing the matrix
	    if(matrix_flag == 0 || matrix_flag == 1) {
	      // loop over basis functions (columns)
	      for(int j=0; j < elems[m].p + 1; j++) {
		int pos_j = elems[m].dof[c_j][j]; // matrix column
		// if j-th basis function is active
		if(pos_j != -1) {
		  // transform j-th basis function to element 'm'
		  this->mesh->element_shapefn(elems[m].x1, elems[m].x2,  
				  j, order, phys_u, phys_dudx); 
		  // evaluate the bilinear form
		  double val_ji = mfv->fn(pts_num, phys_pts,
			    phys_weights, phys_u, phys_dudx, phys_v, phys_dvdx,
			    phys_u_prev, phys_du_prevdx, NULL); 
		  //truncating
		  if (fabs(val_ji) < 1e-12) val_ji = 0.0; 
		  // add the result to the matrix
		  if (val_ji != 0) mat->add(pos_j, pos_i, val_ji);
		  if (DEBUG) {
		    printf("Elem %d: add to matrix pos %d, %d value %g (comp %d, %d)\n", 
		    m, pos_i, pos_j, val_ji, c_i, c_j);
		}
	      }
	    }
	  }
	}
      }
      }}

    // volumetric part of residual
    if(matrix_flag == 0 || matrix_flag == 2) {
      for (int ww = 0; ww < this->vector_forms_vol.size(); ww++)
      {
        VectorFormVol *vfv = &this->vector_forms_vol[ww];
        int c_i = vfv->i;  

        // loop over test functions (rows)
        for(int i=0; i<elems[m].p + 1; i++) {
	  // if i-th test function is active
	  int pos_i = elems[m].dof[c_i][i]; // row in residual vector
	  if(pos_i != -1) {
	    // transform i-th test function to element 'm'
	    this->mesh->element_shapefn(elems[m].x1, elems[m].x2,  
			  i, order, phys_v, phys_dvdx); 

	    // contribute to residual vector
	    if(matrix_flag == 0 || matrix_flag == 2) {
	      double val_i = vfv->fn(pts_num, phys_pts, phys_weights, 
				   phys_u_prev, phys_du_prevdx, phys_v,
				   phys_dvdx, NULL);
	      // truncating
	      if(fabs(val_i) < 1e-12) val_i = 0.0; 
	      // add the contribution to the residual vector
 	      if (val_i != 0) res[pos_i] += val_i;
	      if (DEBUG) {
		if (val_i != 0) {
	          printf("Elem %d: add to residual pos %d value %g (comp %d)\n", 
                  m, pos_i, val_i, c_i);
                }
              }
            }
	  }
	}
      }
    } 
  }
}

// process boundary weak forms
void DiscreteProblem::process_surf_forms(Matrix *mat, double *res, 
					 double *y_prev, int matrix_flag, 
                                         int bdy_index) {
  Element *elems = this->mesh->get_elems();
  // evaluate previous solution and its derivative at the end point
  // FIXME: maximum number of equations limited by [MAX_EQN_NUM][MAX_PTS_NUM]
  double phys_u_prev[MAX_EQN_NUM], 
         phys_du_prevdx[MAX_EQN_NUM]; // at the end point

  // decide whether we are on the left-most or right-most one
  int m;
  if(bdy_index == BOUNDARY_LEFT) m = 0; // first element
  else m = this->mesh->get_n_elems() - 1; // last element

  // calculate coefficients of shape functions on element m
  double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
  this->mesh->calculate_elem_coeffs(m, y_prev, coeffs); 
  double x_ref, x_phys; 
  if(bdy_index == BOUNDARY_LEFT) 
  {
    x_ref = -1; // left end of reference element
    x_phys = this->mesh->left_endpoint;
  } 
  else 
  {
    x_ref = 1;                            // right end of reference element
    x_phys = this->mesh->right_endpoint;
  }

  // get solution value and derivative at the boundary point
  this->mesh->element_solution_point(x_ref, elems + m, coeffs,
                         phys_u_prev, phys_du_prevdx); 

  // surface bilinear forms
  if(matrix_flag == 0 || matrix_flag == 1) {
    for (int ww = 0; ww < this->matrix_forms_surf.size(); ww++)
    {
      MatrixFormSurf *mfs = &this->matrix_forms_surf[ww];
      if (mfs->bdy_index != bdy_index) continue;
      int c_i = mfs->i;  
      int c_j = mfs->j;  

      // loop over test functions on the boundary element
      for(int i=0; i<elems[m].p + 1; i++) {
        double phys_v, phys_dvdx; 
        int pos_i = elems[m].dof[c_i][i]; // matrix row
        if(pos_i != -1) {
          // transform j-th basis function to the boundary element
          this->mesh->element_shapefn_point(x_ref, elems[m].x1, 
                             elems[m].x2, i, &phys_v, 
                             &phys_dvdx); 
          // loop over basis functions on the boundary element
          for(int j=0; j < elems[m].p + 1; j++) {
            double phys_u, phys_dudx;
            int pos_j = elems[m].dof[c_j][j]; // matrix column
            // if j-th basis function is active
            if(pos_j != -1) {
              // transform j-th basis function to the boundary element
              this->mesh->element_shapefn_point(x_ref, elems[m].x1, 
                             elems[m].x2, j, &phys_u, 
                             &phys_dudx); 
              // evaluate the surface bilinear form
              double val_ji_surf = mfs->fn(x_phys,
                               phys_u, phys_dudx, phys_v, 
                               phys_dvdx, phys_u_prev, phys_du_prevdx, 
                               NULL); 
  	      // truncating
	      if(fabs(val_ji_surf) < 1e-12) val_ji_surf = 0.0; 
              // add the result to the matrix
              if (val_ji_surf != 0) mat->add(pos_j, pos_i, val_ji_surf);
            }
          }
	}
      }
    }
  }

  // surface part of residual
  if(matrix_flag == 0 || matrix_flag == 2) {
    for (int ww = 0; ww < this->vector_forms_surf.size(); ww++)
    {
      VectorFormSurf *vfs = &this->vector_forms_surf[ww];
      if (vfs->bdy_index != bdy_index) continue;
      int c_i = vfs->i;  

      // loop over test functions on the boundary element
      for(int i=0; i<elems[m].p + 1; i++) {
        double phys_v, phys_dvdx; 
        int pos_i = elems[m].dof[c_i][i]; // matrix row
        if(pos_i != -1) {
          // transform j-th basis function to the boundary element
          this->mesh->element_shapefn_point(x_ref, elems[m].x1, 
                             elems[m].x2, i, &phys_v, 
                             &phys_dvdx); 
          // evaluate the surface bilinear form
          double val_i_surf = vfs->fn(x_phys,
                          phys_u_prev, phys_du_prevdx, phys_v, phys_dvdx, 
                          NULL); 
          // truncating
          if(fabs(val_i_surf) < 1e-12) val_i_surf = 0.0; 
          // add the result to the matrix
          if (val_i_surf != 0) res[pos_i] += val_i_surf;
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
  process_vol_forms(mat, res, y_prev, matrix_flag);

  // process surface weak forms for the left boundary
  process_surf_forms(mat, res, y_prev, matrix_flag, BOUNDARY_LEFT);

  // process surface weak forms for the right boundary
  process_surf_forms(mat, res, y_prev, matrix_flag, BOUNDARY_RIGHT);

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






