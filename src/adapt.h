// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _ADAPT_H_
#define _ADAPT_H_

#include "common.h"
#include "lobatto.h"
#include "legendre.h"
#include "mesh.h"
#include "matrix.h"
#include "iterator.h"

// Calculates the square in L2-norm of the difference between 
// the coarse and fine mesh solution, for all solution components.
// Assumes that the element was not refined in space for the 
// reference solution. 
// FIXME: to be moved to the Element class
double calc_elem_est_L2_error_squared_p(Element *e, Element *e_ref,
                        double *y_prev, double *y_prev_ref, 
                        double bc_left_dir_values[MAX_EQN_NUM],
			double bc_right_dir_values[MAX_EQN_NUM]);

// Calculates the square in L2-norm of the difference between 
// the coarse and fine mesh solution, for all solution components.
// Assumes that the element was refined in space for the 
// reference solution.
// FIXME: to be moved to the Element class
double calc_elem_est_L2_error_squared_hp(Element *e, 
                        Element *e_ref_left, Element *e_ref_right,
                        double *y_prev, double *y_prev_ref, 
                        double bc_left_dir_values[MAX_EQN_NUM],
			double bc_right_dir_values[MAX_EQN_NUM]);

// Calculates l2-norm (squared) of the difference between the coarse
// and reference solutions in all active elements of 'mesh'. Total
// error is returned.
double calc_elem_est_L2_errors_squared(Mesh* mesh, Mesh* mesh_ref, 
				   double* y_prev, double* y_prev_ref, 
				   double *err_squared_array);

// Can be used for both the coarse and reference solutions
double calc_approx_sol_L2_norm(Mesh* mesh, double* y_prev);

// Sort err_array[] and returning array of sorted element indices
void sort_element_errors(int n, double *err_array, int *id_array); 

// Assumes that reference solution is defined on two half-elements 'e_ref_left' 
// and 'e_ref_right'. The reference solution is projected onto the space of 
// (discontinuous) polynomials of degree 'p_left' on (-1, 0)
// and degree 'p_right' on (0, 1). 
double check_refin_coarse_hp_fine_hp(Element *e, Element *e_ref_left, Element *e_ref_right, 
                                   double *y_prev_ref, int p_left, int p_right, 
                                   double bc_left_dir_values[MAX_EQN_NUM],
	  	                   double bc_right_dir_values[MAX_EQN_NUM]);


// Assumes that reference solution is defined on two half-elements 'e_ref_left'
// and 'e_ref_right'. The reference solution is projected onto the space of 
// (discontinuous) polynomials of degree 'p_left' on (-1, 0)
// and degree 'p_right' on (0, 1). 
double check_refin_coarse_hp_fine_p(Element *e, Element *e_ref,
                                   double *y_prev_ref, int p_left, int p_right, 
                                   double bc_left_dir_values[MAX_EQN_NUM],
	  	                   double bc_right_dir_values[MAX_EQN_NUM]);

// Assumes that reference solution is defined on two half-elements 'e_ref_left'
// and 'e_ref_right'. The reference solution is projected onto the space of 
// polynomials of degree 'p' on (-1, 1). 
double check_refin_coarse_p_fine_hp(Element *e, Element *e_ref_left, Element *e_ref_right, 
                                    double *y_prev_ref, int p,
                                    double bc_left_dir_values[MAX_EQN_NUM],
		                    double bc_right_dir_values[MAX_EQN_NUM]);

// Assumes that reference solution is defined on one single element 
// 'e_ref' (reference refinement did not split the element in space). 
// The reference solution is projected onto the space of 
// polynomials of degree 'p' on (-1, 1). 
double check_refin_coarse_p_fine_p(Element *e, Element *e_ref, 
                                 double *y_prev_ref, int p, 
                                 double bc_left_dir_values[MAX_EQN_NUM],
	  	                 double bc_right_dir_values[MAX_EQN_NUM]);

// Error wrt. exact solution (if provided) on element 'e' 
double calc_elem_exact_L2_error_squared(exact_sol_type exact_sol,
                                        Element *e, double *y_prev, 
                                        double *bc_left_dir_values,
			                double *bc_right_dir_values,
                                        int order);

// Error wrt. exact solution (if provided) on the entire interval (A, B) 
double calc_exact_sol_L2_error(Mesh *mesh, double *y_prev, 
                               exact_sol_type exact_sol, int order); 

// Calculates L2 norm of function exact_sol in interval (A, B)
double calc_exact_sol_L2_norm(exact_sol_type exact_sol, int n_eq, 
                                   double A, double B, int subdivision, 
                                   int order);

#endif
