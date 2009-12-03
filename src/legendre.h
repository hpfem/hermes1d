// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef SHAPESET_LEGENDRE_H_
#define SHAPESET_LEGENDRE_H_

#include <math.h>

#include "common.h"
#include "quad_std.h"

extern double leg_norm_const_ref(int n);
extern void fill_legendre_array(double x, 
                                double val_array[MAX_P+1],
                                double der_array[MAX_P+1]);
extern double legendre_val_ref(double x, int n);
extern double legendre_der_ref(double x, int n);

// Precalculated values of Legendre polynomials and their derivatives 
// at all Gauss quadrature rules on the reference
// interval (-1, 1). The first index runs through Gauss quadrature 
// orders. The second index runs through the quadrature points of 
// the corresponding rule, and the third through the values of 
// Lobatto polynomials at that point. 
extern double legendre_val_ref_tab[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern double legendre_der_ref_tab[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern void precalculate_legendre_1d();
extern int legendre_order_1d[];

#endif /* SHAPESET_LEGENDRE_H_ */
