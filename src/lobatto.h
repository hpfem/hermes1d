// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef SHAPESET_LOBATTO_H_
#define SHAPESET_LOBATTO_H_

#include <math.h>

#include "common.h"

extern shape_fn_t lobatto_fn_tab_1d[];
extern shape_fn_t lobatto_der_tab_1d[];

extern int lobatto_order_1d[];

void fill_lobatto_array(double x, 
			double lobatto_array_val[MAX_P+1],
			double lobatto_array_der[MAX_P+1]);
double calc_lobatto_val(double x, int n);
double calc_lobatto_der(double x, int n);

#endif /* SHAPESET_LOBATTO_H_ */
