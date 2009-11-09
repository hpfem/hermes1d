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

// returns values of Legendre polynomials on (-1,1)
double legendre(int i, double x) {
  return legendre_fn_tab_1d[i](x);
}

// returns values of Legendre polynomials transformed
// and normalized to be orthonormal on (0,1)
double legendre_right(int i, double x) {
  return sqrt(2)*legendre_fn_tab_1d[i](map_right(x));
}

double calc_elem_L2_norm_squared(Element *e, double *y_prev, 
                        double bc_left_dir_values[MAX_EQN_NUM],
			double bc_right_dir_values[MAX_EQN_NUM]) 
{
  int n_eq = e->dof_size;
  double phys_x[MAX_PTS_NUM];          // quad points
  // values of coarse mesh solution for all solution components
  double phys_val_coarse[MAX_EQN_NUM][MAX_PTS_NUM]; 
  double phys_weights[MAX_PTS_NUM];    // quad weights

  // integration order
  int order = 2*e->p;
  int pts_num = 0;

  // create Gauss quadrature on (-1, 1)
  // FIXME: this only transforms the quadrature from (-1, 1) to (-1, 1),
  //        so it should be simplified
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

  // integrate square over (-1, 1)
  double norm_squared[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) {
    norm_squared[c] = 0;
    for (int i=0; i<pts_num; i++) {
      double val = phys_val_coarse[c][i];
      norm_squared[c] += val * val * phys_weights[i];
    }
  }

  double elem_norm_squared = 0;
  for (int c=0; c<n_eq; c++)  
    elem_norm_squared += norm_squared[c];
  return elem_norm_squared;
}

double calc_elem_L2_error_squared_p(Element *e, Element *e_ref, 
                        double *y_prev, double *y_prev_ref, 
                        double bc_left_dir_values[MAX_EQN_NUM],
			double bc_right_dir_values[MAX_EQN_NUM]) 
{
  int n_eq = e->dof_size;
  double phys_x[MAX_PTS_NUM];          // quad points
  // values of coarse mesh solution for all solution components
  double phys_val_coarse[MAX_EQN_NUM][MAX_PTS_NUM]; 
  // values of fine mesh solution for all solution components
  double phys_val_fine[MAX_EQN_NUM][MAX_PTS_NUM];   
  double phys_weights[MAX_PTS_NUM];    // quad weights

  int order = 2*e_ref->p;
  int pts_num = 0;

  // create Gauss quadrature on (-1, 1)
  // FIXME: this only transforms the quadrature from (-1, 1) to (-1, 1),
  // so it should be simplified
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

  double err_squared = 0;
  for (int c=0; c<n_eq; c++)  
    err_squared += norm_squared[c];
  return err_squared;
}

double calc_elem_L2_error_squared_hp(Element *e, 
				   Element *e_ref_left, Element *e_ref_right,
                                   double *y_prev, double *y_prev_ref, 
                                   double bc_left_dir_values[MAX_EQN_NUM],
			           double bc_right_dir_values[MAX_EQN_NUM]) 
{
  int n_eq = e->dof_size;
  double phys_x[MAX_PTS_NUM];          // quad points
  // values of coarse mesh solution for all solution components
  double phys_val_coarse[MAX_EQN_NUM][MAX_PTS_NUM]; 
  // values of fine mesh solution for all solution components
  double phys_val_fine[MAX_EQN_NUM][MAX_PTS_NUM];   
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
  e_ref_left->get_coeffs(y_prev_ref, coeffs, bc_left_dir_values,
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
  e_ref_right->get_coeffs(y_prev_ref, coeffs, bc_left_dir_values,
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

  double err_squared = 0;
  for (int c=0; c<n_eq; c++) {
    err_squared += norm_squared_left[c] + norm_squared_right[c];
  }
  return err_squared;
}

double calc_elem_L2_errors_squared(Mesh* mesh, Mesh* mesh_ref, 
				   double* y_prev, double* y_prev_ref, 
				   double *err_squared_array)
{
  double err_total_squared = 0;
  Iterator *I = new Iterator(mesh);
  Iterator *I_ref = new Iterator(mesh_ref);

  // simultaneous traversal of 'mesh' and 'mesh_ref'
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    Element *e_ref = I_ref->next_active_element();
    double err_squared;
    if (e->level == e_ref->level) { // element 'e' was not refined in space
                                    // for reference solution
      err_squared = calc_elem_L2_error_squared_p(e, e_ref, y_prev, y_prev_ref, 
                                         mesh->bc_left_dir_values,
			                 mesh->bc_right_dir_values);
    }
    else { // element 'e' was refined in space for reference solution
      Element* e_ref_left = e_ref;
      Element* e_ref_right = I_ref->next_active_element();
      err_squared = calc_elem_L2_error_squared_hp(e, e_ref_left, e_ref_right, 
                                          y_prev, y_prev_ref, 
                                          mesh->bc_left_dir_values,
			                  mesh->bc_right_dir_values);
    }
    err_squared_array[e->id] = err_squared;
    err_total_squared += err_squared;
  }
  return sqrt(err_total_squared);
}

double calc_solution_L2_norm(Mesh* mesh, double* y_prev)
{
  double norm_squared = 0;
  Iterator *I = new Iterator(mesh);

  // traverse 'mesh' and calculate squared solution L2 norm 
  // in every element 
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    norm_squared += calc_elem_L2_norm_squared(e, y_prev,
                                         mesh->bc_left_dir_values,
			                 mesh->bc_right_dir_values);
  }
  return sqrt(norm_squared);
}

/* qsort int comparison function */
int int_cmp(const void *a, const void *b)
{
    const double *ia = (const double *)a; // casting pointer types
    const double *ib = (const double *)b;
    return ia[0] - ib[0] < 0; 
	/* double comparison: returns negative if b[0] < a[0] 
	and positive if a[0] < b[0] */
}

// sorting err_array[] and returning array of sorted element indices
void sort_element_errors(int n, double *err_array, int *id_array) 
{
    double array[MAX_ELEM_NUM][2];
    for (int i=0; i<n; i++) {
      array[i][0] = err_array[i];
      array[i][1] = i;
    }

    qsort(array, n, 2*sizeof(double), int_cmp);

    for (int i=0; i<n; i++) id_array[i] = array[i][1];
}

// Assumes that reference solution is defined on two half-elements 'e_ref_left'
// and 'e_ref_right'. The reference solution is projected onto the space of 
// (discontinuous) polynomials of degree 'p_left' on (-1, 0)
// and degree 'p_right' on (0, 1). 
double check_refin_coarse_hp_fine_hp(Element *e, Element *e_ref_left, Element *e_ref_right, 
                                   double *y_prev_ref, int p_left, int p_right, 
                                   double bc_left_dir_values[MAX_EQN_NUM],
		                   double bc_right_dir_values[MAX_EQN_NUM])
{
  int n_eq = e->dof_size;
  double phys_x[MAX_PTS_NUM];                      // quad points
  double leg_pol_values[MAX_P+1][MAX_PTS_NUM];     // values of Legendre polynomials
  // values of coarse mesh solution for all solution components
  double phys_val_coarse[MAX_EQN_NUM][MAX_PTS_NUM];  
  // values of fine mesh solution for all solution components
  double phys_val_fine[MAX_EQN_NUM][MAX_PTS_NUM];   
  double phys_weights[MAX_PTS_NUM];                // quad weights
  double proj_coeffs_left[MAX_EQN_NUM][MAX_P+1];   // for every solution component
  double proj_coeffs_right[MAX_EQN_NUM][MAX_P+1];  // for every solution component

  // first in interval (-1, 0): calculate L2 projection
  // of the reference solution on transformed Legendre 
  // polynomials of degree 'p_left' and calculate error in L2-norm
  int order = 2*max(e_ref_left->p, p_left);
  int pts_num = 0;

  // create Gauss quadrature on (-1, 0)
  create_element_quadrature(-1, 0, order, phys_x, phys_weights, &pts_num); 

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (-1, 0), for every 
  // solution component
  double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
  e_ref_left->get_coeffs(y_prev_ref, coeffs, bc_left_dir_values,
                         bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e_ref_left->get_solution_point(phys_x[i], coeffs, val, der);
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
  // of the reference solution on transformed Legendre polynomials 
  // of degree 'p_right' and calculate error in L2-norm
  order = 2*max(e_ref_right->p, p_right);
  pts_num = 0;

  // create Gauss quadrature on (0, 1)
  create_element_quadrature(0, 1, order, phys_x, phys_weights, &pts_num); 

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (0, 1), for every 
  // solution component
  e_ref_right->get_coeffs(y_prev_ref, coeffs, bc_left_dir_values,
                          bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e_ref_right->get_solution_point(phys_x[i], coeffs, val, der);
    for(int c=0; c<n_eq; c++) phys_val_fine[c][i] = val[c];
  }

  // fill values of the transformed Legendre
  // polynomials in the integration points in (0, 1)
  // NOTE: these are the same as in (-1, 0), nothing to do here.

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
  for (int c=0; c<n_eq; c++) {           // loop over solution components
    for (int j=0; j<pts_num; j++) {      // loop over integration points
      phys_val_coarse[c][j] = 0;
      for (int m=0; m<p_right+1; m++) {  // loop over transf. Leg. polynomials
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
  double error_scaled = -log(err_total); // FIXME: this should involve 
                                         // penalization by the number of 
                                         // newly added degrees of freedom,
                                         // i.e., also poly order on element 'e'
  return error_scaled; 
}

// Assumes that reference solution is defined on two half-elements 'e_ref_left'
// and 'e_ref_right'. The reference solution is projected onto the space of 
// (discontinuous) polynomials of degree 'p_left' on (-1, 0)
// and degree 'p_right' on (0, 1). 
double check_refin_coarse_hp_fine_p(Element *e, Element *e_ref,
                                  double *y_prev_ref, int p_left, int p_right,
                                  double bc_left_dir_values[MAX_EQN_NUM],
		                  double bc_right_dir_values[MAX_EQN_NUM])
{
  int n_eq = e->dof_size;
  double phys_x[MAX_PTS_NUM];                     // quad points
  double leg_pol_values[MAX_P+1][MAX_PTS_NUM];     // values of Legendre polynomials
  // values of coarse mesh solution for all solution components
  double phys_val_coarse[MAX_EQN_NUM][MAX_PTS_NUM];  
  // values of fine mesh solution for all solution components
  double phys_val_fine[MAX_EQN_NUM][MAX_PTS_NUM];   
  double phys_weights[MAX_PTS_NUM];                // quad weights
  double proj_coeffs_left[MAX_EQN_NUM][MAX_P+1];   // for every solution component
  double proj_coeffs_right[MAX_EQN_NUM][MAX_P+1];  // for every solution component

  // first in interval (-1, 0): calculate L2 projection
  // of the reference solution on transformed Legendre 
  // polynomials of degree 'p_left' and calculate error in L2-norm
  int order = 2*max(e_ref->p, p_left);
  int pts_num = 0;

  // create Gauss quadrature on (-1, 0)
  create_element_quadrature(-1, 0, order, phys_x, phys_weights, &pts_num); 

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (-1, 0), for every 
  // solution component
  double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
  e_ref->get_coeffs(y_prev_ref, coeffs, bc_left_dir_values, bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e_ref->get_solution_point(phys_x[i], coeffs, val, der);
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
  // of the reference solution on transformed Legendre polynomials 
  // of degree 'p_right' and calculate error in L2-norm
  order = 2*max(e_ref->p, p_right);
  pts_num = 0;

  // create Gauss quadrature on (0, 1)
  create_element_quadrature(0, 1, order, phys_x, phys_weights, &pts_num); 

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (0, 1), for every 
  // solution component
  e_ref->get_coeffs(y_prev_ref, coeffs, bc_left_dir_values, bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e_ref->get_solution_point(phys_x[i], coeffs, val, der);
    for(int c=0; c<n_eq; c++) phys_val_fine[c][i] = val[c];
  }

  // fill values of the transformed Legendre
  // polynomials in the integration points in (0, 1)
  // NOTE: these are the same as in (-1, 0), nothing to do here.

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
  for (int c=0; c<n_eq; c++) {           // loop over solution components
    for (int j=0; j<pts_num; j++) {      // loop over integration points
      phys_val_coarse[c][j] = 0;
      for (int m=0; m<p_right+1; m++) {  // loop over transf. Leg. polynomials
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
  double error_scaled = -log(err_total); // FIXME: this should involve 
                                         // penalization by the number of 
                                         // newly added degrees of freedom,
                                         // i.e., also poly order on element 'e'
  return error_scaled; 
}

// Assumes that reference solution is defined on two half-elements 'e_ref_left'
// and 'e_ref_right'. The reference solution is projected onto the space of 
// polynomials of degree 'p' on (-1, 1). 
double check_refin_coarse_p_fine_hp(Element *e, Element *e_ref_left, Element *e_ref_right, 
                                    double *y_prev_ref, int p,
                                    double bc_left_dir_values[MAX_EQN_NUM],
		                    double bc_right_dir_values[MAX_EQN_NUM])
{
  int n_eq = e->dof_size;
  double phys_x_left[MAX_PTS_NUM];                          // quad points
  double phys_x_right[MAX_PTS_NUM];                         // quad points
  double leg_pol_values_left[MAX_P+1][MAX_PTS_NUM];         // values of Legendre polynomials
  double leg_pol_values_right[MAX_P+1][MAX_PTS_NUM];        // values of Legendre polynomials
  // values of coarse mesh solution for all solution components
  double phys_val_coarse_left[MAX_EQN_NUM][MAX_PTS_NUM];  
  double phys_val_coarse_right[MAX_EQN_NUM][MAX_PTS_NUM];  
  // values of fine mesh solution for all solution components
  double phys_val_fine_left[MAX_EQN_NUM][MAX_PTS_NUM];   
  double phys_val_fine_right[MAX_EQN_NUM][MAX_PTS_NUM];   
  double phys_weights_left[MAX_PTS_NUM];                    // quad weights
  double phys_weights_right[MAX_PTS_NUM];                   // quad weights
  double proj_coeffs_left[MAX_EQN_NUM][MAX_P+1];            // not the entire coefficient!
  double proj_coeffs_right[MAX_EQN_NUM][MAX_P+1];           // not the entire coefficient!
  double proj_coeffs[MAX_EQN_NUM][MAX_P+1];                 // sum of the left and right parts

  // first in interval (-1, 0): 
  int order_left = 2*max(e_ref_left->p, p);
  int pts_num_left = 0;
  // create Gauss quadrature on (-1, 0)
  create_element_quadrature(-1, 0, order_left, phys_x_left, phys_weights_left, &pts_num_left); 

  // first in interval (0, 1): 
  int order_right = 2*max(e_ref_right->p, p);
  int pts_num_right = 0;
  // create Gauss quadrature on (0, 1)
  create_element_quadrature(0, 1, order_right, phys_x_right, phys_weights_right, &pts_num_right); 

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (-1, 0), for every 
  // solution component
  double coeffs_left[MAX_EQN_NUM][MAX_COEFFS_NUM];
  e_ref_left->get_coeffs(y_prev_ref, coeffs_left, bc_left_dir_values,
                         bc_right_dir_values); 
  for (int i=0; i<pts_num_left; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e_ref_left->get_solution_point(phys_x_left[i], coeffs_left, val, der);
    for(int c=0; c<n_eq; c++) phys_val_fine_left[c][i] = val[c];
  }

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (0, 1), for every 
  // solution component
  double coeffs_right[MAX_EQN_NUM][MAX_COEFFS_NUM];
  e_ref_right->get_coeffs(y_prev_ref, coeffs_right, bc_left_dir_values,
                         bc_right_dir_values); 
  for (int i=0; i<pts_num_right; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e_ref_right->get_solution_point(phys_x_right[i], coeffs_right, val, der);
    for(int c=0; c<n_eq; c++) phys_val_fine_right[c][i] = val[c];
  }

  // fill values of the Legendre
  // polynomials in the integration points in (-1, 0)
  for(int m=0; m < p + 1; m++) {       // loop over Leg. polynomials
    for(int j=0; j<pts_num_left; j++) {     // filling values at integration points
      leg_pol_values_left[m][j] = legendre(m, phys_x_left[j]);
    }
  }

  // fill values of the Legendre
  // polynomials in the integration points in (0, 1)
  for(int m=0; m < p + 1; m++) {       // loop over Leg. polynomials
    for(int j=0; j<pts_num_right; j++) {     // filling values at integration points
      leg_pol_values_right[m][j] = legendre(m, phys_x_right[j]);
    }
  }

  // calculate the projection coefficients on the 
  // interval (-1, 1)... for every 
  // Legendre polynomial and every solution 
  // component. 
  for(int m=0; m<p + 1; m++) { // loop over Leg. polynomials
    for(int c=0; c<n_eq; c++) { // loop over solution components
      proj_coeffs[c][m] = 0;
      for(int j=0; j<pts_num_left; j++) { // loop over integration points left
        proj_coeffs[c][m] += 
          phys_val_fine_left[c][j] * leg_pol_values_left[m][j] * phys_weights_left[j];
      }
      for(int j=0; j<pts_num_right; j++) { // loop over integration points right
        proj_coeffs[c][m] += 
          phys_val_fine_right[c][j] * leg_pol_values_right[m][j] * phys_weights_right[j];
      }
    }
  }

  // evaluate the projection in (-1, 0) for every solution component
  // and every integration point
  for (int c=0; c<n_eq; c++) { // loop over solution components
    for (int j=0; j<pts_num_left; j++) { // loop over integration points left
      phys_val_coarse_left[c][j] = 0;
      for (int m=0; m < p+1; m++) { // loop over Leg. polynomials
        phys_val_coarse_left[c][j] += 
          leg_pol_values_left[m][j] + proj_coeffs[c][m];
      }
    }
  }

  // evaluate the projection in (0, 1) for every solution component
  // and every integration point
  for (int c=0; c<n_eq; c++) { // loop over solution components
    for (int j=0; j<pts_num_right; j++) { // loop over integration points right
      phys_val_coarse_right[c][j] = 0;
      for (int m=0; m < p+1; m++) { // loop over Leg. polynomials
        phys_val_coarse_right[c][j] += 
          leg_pol_values_right[m][j] + proj_coeffs[c][m];
      }
    }
  }

  // calculate the error squared in L2 norm for every solution  
  // component in (-1, 0)
  double err_squared_left[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) { // loop over solution components
    err_squared_left[c] = 0;
    for (int j=0; j<pts_num_left; j++) { // loop over integration points left
      double diff = phys_val_fine_left[c][j] - phys_val_coarse_left[c][j];
      err_squared_left[c] += diff * diff * phys_weights_left[j]; 
    }
  }

  // calculate the error squared in L2 norm for every solution  
  // component in (0, 1)
  double err_squared_right[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) { // loop over solution components
    err_squared_right[c] = 0;
    for (int j=0; j<pts_num_right; j++) { // loop over integration points right
      double diff = phys_val_fine_right[c][j] - phys_val_coarse_right[c][j];
      err_squared_right[c] += diff * diff * phys_weights_right[j]; 
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
  double error_scaled = -log(err_total); // FIXME: this should involve 
                                         // penalization by the number of 
                                         // newly added degrees of freedom,
                                         // i.e., also poly order on element 'e'
  return error_scaled; 
}

// Assumes that reference solution is defined on one single element 
// 'e_ref' (reference refinement did not split the element in space). 
// The reference solution is projected onto the space of 
// polynomials of degree 'p' on (-1, 1). 
double check_ref_coarse_p_fine_p(Element *e, Element *e_ref,
                                 double *y_prev_ref, int p,
                                 double bc_left_dir_values[MAX_EQN_NUM],
		                 double bc_right_dir_values[MAX_EQN_NUM])
{
  int n_eq = e_ref->dof_size;
  double phys_x[MAX_PTS_NUM];          // quad points
  double leg_pol_values[MAX_P+1][MAX_PTS_NUM];  // values of Legendre polynomials
  // values of coarse mesh solution for all solution components
  double phys_val_coarse[MAX_EQN_NUM][MAX_PTS_NUM];   
  // values of fine mesh solution for all solution components
  double phys_val_fine[MAX_EQN_NUM][MAX_PTS_NUM];   
  double phys_weights[MAX_PTS_NUM];    // quad weights
  double proj_coeffs[MAX_EQN_NUM][MAX_P+1]; // for every solution component

  // integration order
  int order = 2*e_ref->p;
  int pts_num = 0;

  // create Gauss quadrature on (-1, 1)
  // FIXME: this transformation is useless
  create_element_quadrature(-1, 1, order, phys_x, phys_weights, &pts_num); 

  // evaluate fine-mesh solution and its derivative 
  // at all quadrature points in (-1, 1), for every 
  // solution component
  double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
  e_ref->get_coeffs(y_prev_ref, coeffs, bc_left_dir_values,
                    bc_right_dir_values); 
  for (int i=0; i<pts_num; i++) {
    double val[MAX_EQN_NUM], der[MAX_EQN_NUM];
    e_ref->get_solution_point(phys_x[i], coeffs, val, der);
    for(int c=0; c<n_eq; c++) phys_val_fine[c][i] = val[c];
  }

  // fill values of the Legendre
  // polynomials in the integration points in (-1, 1)
  for(int m=0; m<p + 1; m++) { // loop over Leg. polynomials
    for(int j=0; j<pts_num; j++) {  // filling values at integration points
      leg_pol_values[m][j] = legendre(m, phys_x[j]);
    }
  }

  // calculate the projection coefficients for every 
  // Legendre polynomial and every solution 
  // component. Since the basis is orthonormal, these 
  // are just integrals of the fine mesh solution with 
  // the Legendre polynomials
  for(int m=0; m<p + 1; m++) { // loop over Leg. polynomials
    for(int c=0; c<n_eq; c++) { // loop over solution components
      proj_coeffs[c][m] = 0;
      for(int j=0; j<pts_num; j++) { // loop over integration points
        proj_coeffs[c][m] += 
          phys_val_fine[c][j] * leg_pol_values[m][j] * phys_weights[j];
      }
    }
  }

  // evaluate the projection in (-1, 1) for every solution component
  // and every integration point
  for (int c=0; c<n_eq; c++) { // loop over solution components
    for (int j=0; j<pts_num; j++) { // loop over integration points
      phys_val_coarse[c][j] = 0;
      for (int m=0; m<p+1; m++) { // loop over Leg. polynomials
        phys_val_coarse[c][j] += 
          leg_pol_values[m][j] + proj_coeffs[c][m];
      }
    }
  }

  // calculate the error squared in L2 norm for every solution  
  // component in (-1, 1)
  double err_squared[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) { // loop over solution components
    err_squared[c] = 0;
    for (int j=0; j<pts_num; j++) { // loop over integration points
      double diff = phys_val_fine[c][j] - phys_val_coarse[c][j];
      err_squared[c] += diff * diff * phys_weights[j]; 
    }
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
  double error_scaled = -log(err_total); // FIXME: this should involve 
                                         // penalization by the number of 
                                         // newly added degrees of freedom,
                                         // i.e., also poly order on element 'e'

  return error_scaled; 
}

// Refine all elements in the id_array list whose id_array >= 0
void refine_elements(Mesh *mesh, Mesh *mesh_ref, double *y_prev, double *y_prev_ref, 
                     int *id_array, double *err_squared_array) 
{
  int adapt_list[MAX_ELEM_NUM];
  int num_to_adapt = 0;

  // Create list of ids of elements to be refined, in increasing order
  for (int i=0; i<mesh->get_n_active_elem(); i++) {
    if (id_array[i] >= 0) {
      adapt_list[num_to_adapt] = i;
      num_to_adapt++;
    }
  }
 
  /*  
  // Debug: Printing list of elements to be refined
  printf("Elements to be adapted:\n");
  for (int i=0; i<num_to_adapt; i++) printf("Elem[%d]\n", adapt_list[i]);
  */

  Iterator *I = new Iterator(mesh);
  Iterator *I_ref = new Iterator(mesh_ref);

  // simultaneous traversal of 'mesh' and 'mesh_ref'
  Element *e;
  int counter = 0;
  while (counter != num_to_adapt) {
    e = I->next_active_element();
    Element *e_ref = I_ref->next_active_element();
    if (e->id == adapt_list[counter]) {
      counter++;
      if (e->level == e_ref->level) { // element 'e' was not refined in space
                                      // for reference solution
        //calc_best_hp_refinement_p(e, e_ref, y_prev, y_prev_ref, 
        //                          mesh->bc_left_dir_values,
	//  		          mesh->bc_right_dir_values);
      }
      else { // element 'e' was refined in space for reference solution
        Element* e_ref_left = e_ref;
        Element* e_ref_right = I_ref->next_active_element();
        //calc_best_hp_refinement_h(e, e_ref_left, e_ref_right, y_prev, y_prev_ref, 
        //                          mesh->bc_left_dir_values,
	//		          mesh->bc_right_dir_values);
      }
      // perform the refinement
      int p_only = 1;
      int p_new = 2;
      int p_new_left = 1, p_new_right = 1;
      if (p_only) e->p += 1;
      else e->refine(p_new_left, p_new_right);
    }
  }
}

