// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "h1_polys.h"


static double h1_polys_fn_0(double x) { return pow(2,(1.0/2.0))/2; }

static double h1_polys_fn_1(double x) { return x*pow(6,(1.0/2.0))/4; }

static double h1_polys_fn_2(double x) { return -3*pow(10,(1.0/2.0))*(1.0/3.0 - pow(x,2))/16; }

static double h1_polys_fn_3(double x) { return 5*pow(2114,(1.0/2.0))*(-9*x/10 + pow(x,3))/302; }

static double h1_polys_fn_4(double x) { return 105*pow(2102,(1.0/2.0))*(27.0/140.0 - 33*pow(x,2)/28 + pow(x,4))/4204; }

static double h1_polys_fn_5(double x) { return 63*pow(50929582,(1.0/2.0))*(445*x/1057 - 1930*pow(x,3)/1359 + pow(x,5))/245296; }

static double h1_polys_fn_6(double x) { return -231*pow(4176478514,(1.0/2.0))*(11345.0/242781.0 - 745*pow(x,2)/1051 + 19215*pow(x,4)/11561 - pow(x,6))/4890848; }

static double h1_polys_fn_7(double x) { return pow(16006427306317305,(1.0/2.0))*(-988505*x/6576999 + 2316195*pow(x,3)/2192333 - 379911*pow(x,5)/199303 + pow(x,7))/24204608; }


shape_fn_t h1_polys_fn_tab_1d[] = {

h1_polys_fn_0,
h1_polys_fn_1,
h1_polys_fn_2,
h1_polys_fn_3,
h1_polys_fn_4,
h1_polys_fn_5,
h1_polys_fn_6,
h1_polys_fn_7,
};



static double h1_polys_der_0(double x) { return 0; }

static double h1_polys_der_1(double x) { return pow(6,(1.0/2.0))/4; }

static double h1_polys_der_2(double x) { return 3*x*pow(10,(1.0/2.0))/8; }

static double h1_polys_der_3(double x) { return -5*pow(2114,(1.0/2.0))*(9.0/10.0 - 3*pow(x,2))/302; }

static double h1_polys_der_4(double x) { return 105*pow(2102,(1.0/2.0))*(-33*x/14 + 4*pow(x,3))/4204; }

static double h1_polys_der_5(double x) { return 63*pow(50929582,(1.0/2.0))*(445.0/1057.0 - 1930*pow(x,2)/453 + 5*pow(x,4))/245296; }

static double h1_polys_der_6(double x) { return -231*pow(4176478514,(1.0/2.0))*(-1490*x/1051 + 76860*pow(x,3)/11561 - 6*pow(x,5))/4890848; }

static double h1_polys_der_7(double x) { return -pow(16006427306317305,(1.0/2.0))*(988505.0/6576999.0 - 6948585*pow(x,2)/2192333 + 1899555*pow(x,4)/199303 - 7*pow(x,6))/24204608; }


shape_fn_t h1_polys_der_tab_1d[] = {

h1_polys_der_0,
h1_polys_der_1,
h1_polys_der_2,
h1_polys_der_3,
h1_polys_der_4,
h1_polys_der_5,
h1_polys_der_6,
h1_polys_der_7,
};