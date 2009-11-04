// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "transforms.h"

// transform values from (-1, 0) to (-1, 1)
#define map_left(x) (2*x+1)
// transform values from (0, 1) to (-1, 1)
#define map_right(x) (2*x-1)

// returns values of Legendre polynomials transformed
// and normalized to be orthonormal on (-1,0)
double legendre_left(int i, double x) {
  return sqrt(2)*legendre_fn_tab_1d[i](map_left(x));
}

// returns values of Legendre polynomials transformed
// and normalized to be orthonormal on (0,1)
double legendre_right(int i, double x) {
  return sqrt(2)*legendre_fn_tab_1d[i](map_right(x));
}

void calc_elem_L2_error_squared_p(Element *e, Element *e_ref, 
                        double *y_prev, double *y_prev_ref, 
                        double bc_left_dir_values[MAX_EQN_NUM],
			double bc_right_dir_values[MAX_EQN_NUM]) 
{
  int n_eq = e->dof_size;
  double phys_x[MAX_PTS_NUM];          // quad points
  double phys_val_coarse[MAX_EQN_NUM][MAX_PTS_NUM]; // values of coarse mesh solution for all solution components
  double phys_val_fine[MAX_EQN_NUM][MAX_PTS_NUM];   // values of fine mesh solution for all solution components
  double phys_weights[MAX_PTS_NUM];    // quad weights

  // first process interval (-1, 0)
  int order = 2*e_ref->p;
  int pts_num = 0;

  // create Gauss quadrature on (-1, 1)
  // FIXME: this transforms the quadrature from (-1, 1) to (-1, 1)
  create_element_quadrature(-1, 1, order, phys_x, phys_weights,
                            &pts_num); 

  // evaluate coarse-mesh solution and its derivative 
  // at all quadrature points in (-1, 1), for every 
  // solution component
  double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
  e->get_coeffs(y_prev, coeffs, bc_left_dir_values,
                bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e->get_solution_point(phys_x[i], coeffs, val, der);
    for(int c=0; c<n_eq; c++) phys_val_coarse[c][i] = val[c];
  }

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (-1, 1), for every 
  // solution component
  e_ref->get_coeffs(y_prev, coeffs, bc_left_dir_values,
                         bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e_ref->get_solution_point(phys_x[i], coeffs, 
		       	           val, der);
    for(int c=0; c<n_eq; c++) phys_val_fine[c][i] = val[c];
  }

  // integrate over (-1, 1)
  double norm_squared[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) {
    norm_squared[c] = 0;
    for (int i=0; i<pts_num; i++) {
      double diff = phys_val_fine[c][i] - phys_val_coarse[c][i];
      norm_squared[c] += diff*diff*phys_weights[i];
    }
  }

  e->err_squared = 0;
  for (int c=0; c<n_eq; c++)  
    e->err_squared += norm_squared[c];
}

void calc_elem_L2_error_squared_hp(Element *e, 
				   Element *e_ref_left, Element *e_ref_right,
                                   double *y_prev, double *y_prev_ref, 
                                   double bc_left_dir_values[MAX_EQN_NUM],
			           double bc_right_dir_values[MAX_EQN_NUM]) 
{
  int n_eq = e->dof_size;
  double phys_x[MAX_PTS_NUM];          // quad points
  double phys_val_coarse[MAX_EQN_NUM][MAX_PTS_NUM]; // values of coarse mesh solution for all solution components
  double phys_val_fine[MAX_EQN_NUM][MAX_PTS_NUM];   // values of fine mesh solution for all solution components
  double phys_weights[MAX_PTS_NUM];    // quad weights

  // first process interval (-1, 0)
  int order = 2*e_ref_left->p;
  int pts_num = 0;

  // create Gauss quadrature on (-1, 0)
  create_element_quadrature(-1, 0, order, phys_x, phys_weights,
                            &pts_num); 

  // evaluate coarse-mesh solution and its derivative 
  // at all quadrature points in (-1, 0), for every 
  // solution component
  double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
  e->get_coeffs(y_prev, coeffs, bc_left_dir_values,
                bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e->get_solution_point(phys_x[i], coeffs, val, der);
    for(int c=0; c<n_eq; c++) phys_val_coarse[c][i] = val[c];
  }

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (-1, 0), for every 
  // solution component
  e_ref_left->get_coeffs(y_prev, coeffs, bc_left_dir_values,
                         bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e_ref_left->get_solution_point(phys_x[i], coeffs, 
		       	           val, der);
    for(int c=0; c<n_eq; c++) phys_val_fine[c][i] = val[c];
  }

  // integrate over (-1, 0)
  double norm_squared_left[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) {
    norm_squared_left[c] = 0;
    for (int i=0; i<pts_num; i++) {
      double diff = phys_val_fine[c][i] - phys_val_coarse[c][i];
      norm_squared_left[c] += diff*diff*phys_weights[i];
    }
  }

  // next process interval (0, 1)
  order = 2*e_ref_right->p;
  pts_num = 0;

  // create Gauss quadrature on (0, 1)
  create_element_quadrature(0, 1, order, phys_x, phys_weights,
                            &pts_num); 

  // evaluate coarse-mesh solution and its derivative 
  // at all quadrature points in (0, 1), for every 
  // solution component
  e->get_coeffs(y_prev, coeffs, bc_left_dir_values,
                bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e->get_solution_point(phys_x[i], coeffs, val, der);
    for(int c=0; c<n_eq; c++) phys_val_coarse[c][i] = val[c];
  }

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (0, 1), for every 
  // solution component
  e_ref_right->get_coeffs(y_prev, coeffs, bc_left_dir_values,
                         bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e_ref_right->get_solution_point(phys_x[i], coeffs, 
		       	           val, der);
    for(int c=0; c<n_eq; c++) phys_val_fine[c][i] = val[c];
  }

  // integrate over (0, 1)
  double norm_squared_right[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) {
    norm_squared_right[c] = 0;
    for (int i=0; i<pts_num; i++) {
      double diff = phys_val_fine[c][i] - phys_val_coarse[c][i];
      norm_squared_right[c] += diff*diff*phys_weights[i];
    }
  }

  e->err_squared = 0;
  for (int c=0; c<n_eq; c++)  
    e->err_squared += norm_squared_left[c] + norm_squared_right[c];
}

void calc_elem_L2_errors_squared(Mesh* mesh, Mesh* mesh_ref, 
				 double* y_prev, double* y_prev_ref) 
{
  Iterator *I = new Iterator(mesh);
  Iterator *I_ref = new Iterator(mesh_ref);

  // simultaneous traversal of 'mesh' and 'mesh_ref'
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    Element *e_ref = I_ref->next_active_element();
    if (e->level == e_ref->level) { // element 'e' was not refined in space
                                    // for reference solution
      calc_elem_L2_error_squared_p(e, e_ref, y_prev, y_prev_ref, 
                                   mesh->bc_left_dir_values,
			           mesh->bc_right_dir_values);
    }
    else { // element 'e' was refined in space for reference solution
      Element* e_ref_left = e_ref;
      Element* e_ref_right = I_ref->next_active_element();
      calc_elem_L2_error_squared_hp(e, e_ref_left, e_ref_right, 
                                    y_prev, y_prev_ref, 
                                    mesh->bc_left_dir_values,
			            mesh->bc_right_dir_values);
    }
  }
}

void print_element_errors(Mesh *mesh) 
{
  Iterator *I = new Iterator(mesh);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    printf("Elem %d, err_squared = %g\n", e->id, e->err_squared);
  }
}


// projects reference solution on element 'e' onto the space of 
// (discontinuous) polynomials of degree 'p_left' on (-1, 0)
// and degree 'p_right' on (0, 1). Returns the projection error,
// penalized by the number of degrees of freedom of the 
// candidate
double check_hp_candidate(Element *e, double *y_prev_ref, 
                          int p_left, int p_right, 
                          double bc_left_dir_values[MAX_EQN_NUM],
		          double bc_right_dir_values[MAX_EQN_NUM])
{
  int n_eq = e->dof_size;
  double phys_x[MAX_PTS_NUM];          // quad points
  double leg_pol_values[MAX_P+1][MAX_PTS_NUM];  // values of Legendre polynomials
  double phys_val_coarse[MAX_EQN_NUM][MAX_PTS_NUM];   // values of coarse mesh solution for all solution components
  double phys_val_fine[MAX_EQN_NUM][MAX_PTS_NUM];   // values of fine mesh solution for all solution components
  double phys_weights[MAX_PTS_NUM];    // quad weights
  double proj_coeffs_left[MAX_EQN_NUM][MAX_P+1]; // for every solution component
  double proj_coeffs_right[MAX_EQN_NUM][MAX_P+1];// for every solution component

  // first in interval (-1, 0): calculate L2 projection 
  // of the reference solution on polynomials of degree 
  // 'p_left' and calculate error in L2-norm
  int order = 2*e->sons[0]->p;
  int pts_num = 0;

  // create Gauss quadrature on (-1, 0)
  create_element_quadrature(-1, 0, order, phys_x, phys_weights,
                            &pts_num); 

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (-1, 0), for every 
  // solution component
  double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
  e->sons[0]->get_coeffs(y_prev_ref, coeffs, bc_left_dir_values,
                         bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e->sons[0]->get_solution_point(phys_x[i], coeffs, 
		       	           val, der);
    for(int c=0; c<n_eq; c++) phys_val_fine[c][i] = val[c];
  }

  // fill values of the transformed Legendre
  // polynomials in the integration points in (-1, 0)
  for(int m=0; m<p_left + 1; m++) { // loop over transf. Leg. polynomials
    for(int j=0; j<pts_num; j++) { // filling values at integration points
      leg_pol_values[m][j] = legendre_left(m, phys_x[j]);
    }
  }

  // calculate the projection coefficients for every 
  // transformed Legendre polynomial and every solution 
  // component. Since the basis is orthonormal, these 
  // are just integrals of the fine mesh solution with 
  // the transformed Legendre polynomials
  for(int m=0; m<p_left + 1; m++) { // loop over transf. Leg. polynomials
    for(int c=0; c<n_eq; c++) { // loop over solution components
      proj_coeffs_left[c][m] = 0;
      for(int j=0; j<pts_num; j++) { // loop over integration points
        proj_coeffs_left[c][m] += 
          phys_val_fine[c][j] * leg_pol_values[m][j] * phys_weights[j];
      }
    }
  }

  // evaluate the projection in (-1, 0) for every solution component
  // and every integration point
  for (int c=0; c<n_eq; c++) { // loop over solution components
    for (int j=0; j<pts_num; j++) { // loop over integration points
      phys_val_coarse[c][j] = 0;
      for (int m=0; m<p_left+1; m++) { // loop over transf. Leg. polynomials
        phys_val_coarse[c][j] += 
          leg_pol_values[m][j] + proj_coeffs_left[c][m];
      }
    }
  }

  // calculate the error squared in L2 norm for every solution  
  // component in (-1, 0)
  double err_squared_left[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) { // loop over solution components
    err_squared_left[c] = 0;
    for (int j=0; j<pts_num; j++) { // loop over integration points
      double diff = phys_val_fine[c][j] - phys_val_coarse[c][j];
      err_squared_left[c] += diff * diff * phys_weights[j]; 
    }
  }

  // next in interval (0, 1): calculate L2 projection 
  // of the reference solution on polynomials of degree 
  // 'p_right' and calculate error in L2-norm
  order = 2*e->sons[1]->p;
  pts_num = 0;

  // create Gauss quadrature on (0, 1)
  create_element_quadrature(0, 1, order, phys_x, phys_weights,
                            &pts_num); 

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (0, 1), for every 
  // solution component
  e->sons[1]->get_coeffs(y_prev_ref, coeffs, bc_left_dir_values,
                         bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e->sons[1]->get_solution_point(phys_x[i], coeffs, 
		       	           val, der);
    for(int c=0; c<n_eq; c++) phys_val_fine[c][i] = val[c];
  }

  // fill values of the transformed Legendre
  // polynomials in the integration points in (0, 1)
  // NOTE: these are the same as in (-1, 0), nothing to do!

  // calculate the projection coefficients for every 
  // transformed Legendre polynomial and every solution 
  // component. Since the basis is orthonormal, these 
  // are just integrals of the fine mesh solution with 
  // the transformed Legendre polynomials
  for(int m=0; m<p_right + 1; m++) { // loop over transf. Leg. polynomials
    for(int c=0; c<n_eq; c++) { // loop over solution components
      proj_coeffs_right[c][m] = 0;
      for(int j=0; j<pts_num; j++) { // loop over integration points
        proj_coeffs_right[c][m] += 
          phys_val_fine[c][j] * leg_pol_values[m][j] * phys_weights[j];
      }
    }
  }

  // evaluate the projection in (0, 1) for every solution component
  // and every integration point
  for (int c=0; c<n_eq; c++) { // loop over solution components
    for (int j=0; j<pts_num; j++) { // loop over integration points
      phys_val_coarse[c][j] = 0;
      for (int m=0; m<p_right+1; m++) { // loop over transf. Leg. polynomials
        phys_val_coarse[c][j] += 
          leg_pol_values[m][j] + proj_coeffs_right[c][m];
      }
    }
  }

  // calculate the error squared in L2 norm for every solution  
  // component in (0, 1)
  double err_squared_right[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) { // loop over solution components
    err_squared_right[c] = 0;
    for (int j=0; j<pts_num; j++) { // loop over integration points
      double diff = phys_val_fine[c][j] - phys_val_coarse[c][j];
      err_squared_right[c] += diff * diff * phys_weights[j]; 
    }
  }

  // calculating resulting error squared for 
  // every solution component
  double err_squared[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) { // loop over solution components
    err_squared[c] = err_squared_left[c] + err_squared_right[c];
  }

  // summing the errors squared over all components on that 
  // element, and taking a square root. NOTE: this is just 
  // one of many possible ways to merge the errors together
  double err_total = 0;
  for (int c=0; c<n_eq; c++) { // loop over solution components
    err_total += err_squared[c];
  }
  err_total = sqrt(err_total);

  // penalizing the error by the number of DOF that this 
  // candidate would bring to the system
  // NOTE: this is the most tricky part. Check how this is 
  // done in Hermes2D
  double error_scaled = -log(err_total); // this assumes that all
                                         // candidates bring the same 
                                         // number of new DOF 

  return error_scaled; 
}
