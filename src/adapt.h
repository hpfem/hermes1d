// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _ADAPT_H_
#define _ADAPT_H_

#include "common.h"
#include "lobatto.h"
#include "mesh.h"
#include "matrix.h"
#include "iterator.h"

// calculates the square in L2-norm of the difference between 
// the coarse and fine mesh solution, for all solution components
void calc_elem_L2_error_squared(Element *e, double *y_prev, double
			*y_prev_ref, 
                        double bc_left_dir_values[MAX_EQN_NUM],
			double bc_right_dir_values[MAX_EQN_NUM],
			double norm_squared[MAX_EQN_NUM]);

// projects reference solution on element 'e' onto a space of 
// (discontinuous) polynomials of degree p_left on the left
// son and degree p_right on the right son  
void try_hp_candidate(Element *e, double *y_prev_ref, 
                      int p_left, int p_right, double *error, 
                      int *new_ndof);

#endif
