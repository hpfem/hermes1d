// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "quad_std.h"

// transformation of quadrature to physical element
void create_element_quadrature(double a, double b, 
                        int order, double *pts, double *weights, int *num) {
  double2 *ref_tab = g_quad_1d_std.get_points(order);
  *num = g_quad_1d_std.get_num_points(order);
  for (int i=0;i<*num;i++) {
    //change points and weights to interval (a, b)
    pts[i] = (b-a)/2.*ref_tab[i][0]+(b+a)/2.; 
    weights[i] = ref_tab[i][1]*(b-a)/2.;
  }
};

Quad1DStd::Quad1DStd()
{
  tables = std_tables_1d;
  np = std_np_1d;
  ref_vert[0] = -1.0;
  ref_vert[1] = 1.0;
  max_order = 99;
}

Quad1DStd g_quad_1d_std;
