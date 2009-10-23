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

// calculates the difference between the coarse and fine 
// mesh solution on element 'e' in L2 norm 
void calc_elem_L2_error(Element *e, double *y_prev, double
        *y_prev_ref);

// projects reference solution on element 'e' onto a space of 
// (discontinuous) polynomials of degree p_left on the left
// son and degree p_right on the right son  
void try_hp_candidate(Element *e, double *y_prev_ref, 
                      int p_left, int p_right, double *error, 
                      int *new_ndof);

#endif
