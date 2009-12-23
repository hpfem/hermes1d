// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "discrete.h"
#include "solver_umfpack.h"

DiscreteProblem::DiscreteProblem() {
  // precalculating values and derivatives 
  // of all polynomials at all possible 
  // integration points
  fprintf(stderr, "Precalculating Legendre polynomials...");
  fflush(stderr);
  precalculate_legendre_1d();
  precalculate_legendre_1d_left();
  precalculate_legendre_1d_right();
  fprintf(stderr, "done.\n");

  fprintf(stderr, "Precalculating Lobatto shape functions...");
  fflush(stderr);
  precalculate_lobatto_1d();
  precalculate_lobatto_1d_left();
  precalculate_lobatto_1d_right();
  fprintf(stderr, "done.\n");
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
void DiscreteProblem::process_vol_forms(Mesh *mesh, Matrix *mat, double *res, 
					int matrix_flag) {
  int n_eq = mesh->get_n_eq();
  Element *elems = mesh->get_base_elems();
  int n_elem = mesh->get_n_base_elem();
  Iterator *I = new Iterator(mesh);

  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    //printf("Processing elem %d\n", m);
    int    pts_num;                                     // num of quad points
    double phys_pts[MAX_QUAD_PTS_NUM];                  // quad points
    double phys_weights[MAX_QUAD_PTS_NUM];              // quad weights
    double phys_u[MAX_QUAD_PTS_NUM];                    // basis function 
    double phys_dudx[MAX_QUAD_PTS_NUM];                 // basis function x-derivative
    double phys_v[MAX_QUAD_PTS_NUM];                    // test function
    double phys_dvdx[MAX_QUAD_PTS_NUM];                 // test function x-derivative
    if (n_eq > MAX_EQN_NUM) error("number of equations exceeded in process_vol_forms().");
    double phys_u_prev[MAX_EQN_NUM][MAX_QUAD_PTS_NUM];     // previous solution, all components
    double phys_du_prevdx[MAX_EQN_NUM][MAX_QUAD_PTS_NUM];  // previous solution x-derivative, all components
    // decide quadrature order and set up 
    // quadrature weights and points in element m
    // CAUTION: This is heuristic
    int order = 4*e->p;

    // prepare quadrature points and weights in element 'e'
    create_phys_element_quadrature(e->x1, e->x2,  
                               order, phys_pts, phys_weights, &pts_num); 

    // evaluate previous solution and its derivative 
    // at all quadrature points in the element, 
    // for every solution component
    // 0... in the entire element
    e->get_solution_quad(0, order,
                         phys_u_prev, phys_du_prevdx); 

    // volumetric bilinear forms
    if(matrix_flag == 0 || matrix_flag == 1) 
    {
      for (int ww = 0; ww < this->matrix_forms_vol.size(); ww++)
      {
	MatrixFormVol *mfv = &this->matrix_forms_vol[ww];
	int c_i = mfv->i;  
	int c_j = mfv->j;  

	// loop over test functions (rows)
	for(int i=0; i<e->p + 1; i++) {
	  // if i-th test function is active
	  int pos_i = e->dof[c_i][i]; // row in matrix
          //printf("elem (%g, %g): pos_i = %d\n", e->x1, e->x2, pos_i);
	  if(pos_i != -1) {
	    // transform i-th test function to element 'm'
            //printf("Elem (%g, %g): i = %d, order = %d\n", e->x1, e->x2, i, order);
	    element_shapefn(e->x1, e->x2,  
			    i, order, phys_v, phys_dvdx); 
	    // if we are constructing the matrix
	    if(matrix_flag == 0 || matrix_flag == 1) {
	      // loop over basis functions (columns)
	      for(int j=0; j < e->p + 1; j++) {
		int pos_j = e->dof[c_j][j]; // matrix column
                //printf("elem (%g, %g): pos_j = %d\n", e->x1, e->x2, pos_j);
		// if j-th basis function is active
		if(pos_j != -1) {
		  // transform j-th basis function to element 'm'
		  element_shapefn(e->x1, e->x2,  
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
		    printf("Adding to matrix pos %d, %d value %g (comp %d, %d)\n", 
		    pos_i, pos_j, val_ji, c_i, c_j);
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
        for(int i=0; i<e->p + 1; i++) {
	  // if i-th test function is active
	  int pos_i = e->dof[c_i][i]; // row in residual vector
	  if(pos_i != -1) {
	    // transform i-th test function to element 'm'
	    element_shapefn(e->x1, e->x2,  
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
	          printf("Adding to residual pos %d value %g (comp %d)\n", 
                  pos_i, val_i, c_i);
                }
              }
            }
	  }
	}
      }
    }
  } // end while
  delete I;
}

// process boundary weak forms
void DiscreteProblem::process_surf_forms(Mesh *mesh, Matrix *mat, double *res, 
					 int matrix_flag, int bdy_index) {
  Iterator *I = new Iterator(mesh);
  Element *e; 

  // evaluate previous solution and its derivative at the end point
  // FIXME: maximum number of equations limited by [MAX_EQN_NUM][MAX_QUAD_PTS_NUM]
  double phys_u_prev[MAX_EQN_NUM], 
         phys_du_prevdx[MAX_EQN_NUM]; // at the end point

  // decide whether we are on the left-most or right-most one
  double x_ref, x_phys; 
  if(bdy_index == BOUNDARY_LEFT) {
    e = I->first_active_element(); 
    x_ref = -1; // left end of reference element
    x_phys = mesh->get_left_endpoint();
  }
  else {
    e = I->last_active_element(); 
    x_ref = 1;  // right end of reference element
    x_phys = mesh->get_right_endpoint();
  }

  // get solution value and derivative at the boundary point
  e->get_solution_point(x_phys, phys_u_prev, phys_du_prevdx); 

  // surface bilinear forms
  if(matrix_flag == 0 || matrix_flag == 1) {
    for (int ww = 0; ww < this->matrix_forms_surf.size(); ww++)
    {
      MatrixFormSurf *mfs = &this->matrix_forms_surf[ww];
      if (mfs->bdy_index != bdy_index) continue;
      int c_i = mfs->i;  
      int c_j = mfs->j;  

      // loop over test functions on the boundary element
      for(int i=0; i<e->p + 1; i++) {
        double phys_v, phys_dvdx; 
        int pos_i = e->dof[c_i][i]; // matrix row
        if(pos_i != -1) {
          // transform j-th basis function to the boundary element
          element_shapefn_point(x_ref, e->x1, e->x2, i, phys_v, 
                                phys_dvdx); 
          // loop over basis functions on the boundary element
          for(int j=0; j < e->p + 1; j++) {
            double phys_u, phys_dudx;
            int pos_j = e->dof[c_j][j]; // matrix column
            // if j-th basis function is active
            if(pos_j != -1) {
              // transform j-th basis function to the boundary element
              element_shapefn_point(x_ref, e->x1, e->x2, j, phys_u, 
                                    phys_dudx); 
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
      for(int i=0; i<e->p + 1; i++) {
        double phys_v, phys_dvdx; 
        int pos_i = e->dof[c_i][i]; // matrix row
        if(pos_i != -1) {
          // transform j-th basis function to the boundary element
          element_shapefn_point(x_ref, e->x1, e->x2, i, phys_v, 
                                phys_dvdx); 
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
  delete I;
}

// construct Jacobi matrix or residual vector
// matrix_flag == 0... assembling Jacobi matrix and residual vector together
// matrix_flag == 1... assembling Jacobi matrix only
// matrix_flag == 2... assembling residual vector only
// NOTE: Simultaneous assembling of the Jacobi matrix and residual
// vector is more efficient than if they are assembled separately
void DiscreteProblem::assemble(Mesh *mesh, Matrix *mat, double *res, 
                               int matrix_flag) {
  // number of equations in the system
  int n_eq = mesh->get_n_eq();

  // total number of unknowns
  int n_dof = mesh->get_n_dof();

  // erase residual vector
  if(matrix_flag == 0 || matrix_flag == 2) 
    for(int i=0; i<n_dof; i++) res[i] = 0;

  // process volumetric weak forms via an element loop
  process_vol_forms(mesh, mat, res, matrix_flag);

  // process surface weak forms for the left boundary
  process_surf_forms(mesh, mat, res, matrix_flag, BOUNDARY_LEFT);

  // process surface weak forms for the right boundary
  process_surf_forms(mesh, mat, res, matrix_flag, BOUNDARY_RIGHT);

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
void DiscreteProblem::assemble_matrix_and_vector(Mesh *mesh, 
                      Matrix *mat, double *res) {
  assemble(mesh, mat, res, 0);
} 

// construct Jacobi matrix only
void DiscreteProblem::assemble_matrix(Mesh *mesh, Matrix *mat) {
  double *void_res = NULL;
  assemble(mesh, mat, void_res, 1);
} 

// construct residual vector only
void DiscreteProblem::assemble_vector(Mesh *mesh, double *res) {
  Matrix *void_mat = NULL;
  assemble(mesh, void_mat, res, 2);
} 

void solve_linear_system_cg(Matrix* mat, double *res) 
{
  error("notimplemented yet.");  

}

void copy_mesh_to_vector(Mesh *mesh, double *y) {
  Element *e;
  Iterator *I = new Iterator(mesh);
  while ((e = I->next_active_element()) != NULL) {
    e->copy_coeffs_to_vector(y);
  }
  delete I;
}
void copy_vector_to_mesh(double *y, Mesh *mesh) {
  Element *e;
  Iterator *I = new Iterator(mesh);
  while ((e = I->next_active_element()) != NULL) {
    e->get_coeffs_from_vector(y);
  }
  delete I;
}

// Newton's iteration
int newton(int solver, DiscreteProblem *dp, Mesh *mesh, 
           double tol_newton, int &iter_num) 
{
  iter_num = 1;
  int n_dof = mesh->get_n_dof();
  double *y = new double[n_dof];
  if (y == NULL) error("vector y could not be allocated in newton().");
  double *res = new double[n_dof];
  if (res == NULL)
    error("res could not be allocated in newton().");

  // fill vector y using dof and coeffs arrays 
  // in elements
  copy_mesh_to_vector(mesh, y);

  // the iteration
  CooMatrix *mat = NULL;
  while (1) {
    // Reset the matrix:
    if (mat != NULL) delete mat;
    mat = new CooMatrix();

    // construct matrix and residual vector
    dp->assemble_matrix_and_vector(mesh, mat, res); 

    // calculate L2 norm of residual vector
    double res_norm = 0;
    for(int i=0; i<n_dof; i++) res_norm += res[i]*res[i];
    res_norm = sqrt(res_norm);

    // If residual norm less than 'tol_newton', quit
    // latest solution is in the vector y.
    // CAUTION: at least one full iteration forced
    //          here because sometimes the initial
    //          residual on fine mesh is too small
    printf("Residual norm: %.15f\n", res_norm);
    if(res_norm < tol_newton && iter_num > 1) break;

    // changing sign of vector res
    for(int i=0; i<n_dof; i++) res[i]*= -1;

    // solving the matrix system
    //solve_linear_system_umfpack((CooMatrix*)mat, res);
    if (solver == 0) solve_linear_system_umfpack((CooMatrix*)mat, res);
    else {
      if (solver == 1) solve_linear_system_cg((CooMatrix*)mat, res);
      else error("unknown matrix solver in newton().");
    }

    // updating vector y by new solution which is in res
    for(int i=0; i<n_dof; i++) y[i] += res[i];

    // copy coefficients from vector y to elements
    copy_vector_to_mesh(y, mesh);

    iter_num++;
    if (iter_num >= MAX_NEWTON_ITER_NUM) 
      return 0; // no success
  }

  if (mat != NULL) delete mat;
  if (y != NULL) delete [] y;
  if (res != NULL) delete [] res;

  // finished successfully
  return 1;
}




