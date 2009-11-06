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
double calc_elem_L2_error_squared_p(Element *e, Element *e_ref,
                        double *y_prev, double *y_prev_ref, 
                        double bc_left_dir_values[MAX_EQN_NUM],
			double bc_right_dir_values[MAX_EQN_NUM]);

// Calculates the square in L2-norm of the difference between 
// the coarse and fine mesh solution, for all solution components.
// Assumes that the element was refined in space for the 
// reference solution.
// FIXME: to be moved to the Element class
double calc_elem_L2_error_squared_hp(Element *e, 
                        Element *e_ref_left, Element *e_ref_right,
                        double *y_prev, double *y_prev_ref, 
                        double bc_left_dir_values[MAX_EQN_NUM],
			double bc_right_dir_values[MAX_EQN_NUM]);

// Calculates l2-norm (squared) of the difference between the coarse
// and fine mesh solutions in all active elements of 'mesh'. The 
// calculated values are stored in the variable 'err_squared' in
// every element of 'mesh' 
void calc_elem_L2_errors_squared(Mesh* mesh, Mesh* mesh_ref, 
				 double* y_prev, double* y_prev_ref, 
                                 double *err_array);

// sorting err_array[] and returning array of sorted element indices
void sort_element_errors(int n, double *err_array, int *id_array); 

// Projects reference solution on element 'e' onto the space of 
// (discontinuous) polynomials of degree 'p_left' on (-1, 0)
// and degree 'p_right' on (0, 1)
double check_hp_candidate(Element *e, double *y_prev_ref, 
                        int p_left, int p_right, 
                        double bc_left_dir_values[MAX_EQN_NUM],
	  	        double bc_right_dir_values[MAX_EQN_NUM]);

#endif
