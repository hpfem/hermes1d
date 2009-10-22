// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "lobatto.h"

int lobatto_order_1d[] = {
1, {% for f in functions[1:] %}
{{ f.id }},{% endfor %}
};
int legendre_order_1d[] = {
{% for f in functions %}
{{ f.id }},{% endfor %}
};

{% for f in functions %}
static double lobatto_fn_{{ f.id }}(double x) { return {{ f.expr }}; }
{% endfor %}

shape_fn_t lobatto_fn_tab_1d[] = {
{% for f in functions %}
lobatto_fn_{{ f.id }},{% endfor %}
};


{% for f in functions %}
static double lobatto_der_{{ f.id }}(double x) { return {{ f.expr_diff }}; }
{% endfor %}

shape_fn_t lobatto_der_tab_1d[] = {
{% for f in functions %}
lobatto_der_{{ f.id }},{% endfor %}
};


static double legendre_fn_0(double x) { return  legendre0(x); }
static double legendre_fn_1(double x) { return  legendre1(x); }
static double legendre_fn_2(double x) { return  legendre2(x); }
static double legendre_fn_3(double x) { return  legendre3(x); }
static double legendre_fn_4(double x) { return  legendre4(x); }
static double legendre_fn_5(double x) { return  legendre5(x); }
static double legendre_fn_6(double x) { return  legendre6(x); }
static double legendre_fn_7(double x) { return  legendre7(x); }
static double legendre_fn_8(double x) { return  legendre8(x); }
static double legendre_fn_9(double x) { return  legendre9(x); }
static double legendre_fn_10(double x) { return  legendre10(x); }
static double legendre_fn_11(double x) { return  legendre11(x); }

shape_fn_t legendre_fn_tab_1d[] = {
	legendre_fn_0, legendre_fn_1, legendre_fn_2, legendre_fn_3, legendre_fn_4, legendre_fn_5,
	legendre_fn_6, legendre_fn_7, legendre_fn_8, legendre_fn_9, legendre_fn_10, legendre_fn_11
};


static double legendre_der_0(double x) { return  legendre0x(x); }
static double legendre_der_1(double x) { return  legendre1x(x); }
static double legendre_der_2(double x) { return  legendre2x(x); }
static double legendre_der_3(double x) { return  legendre3x(x); }
static double legendre_der_4(double x) { return  legendre4x(x); }
static double legendre_der_5(double x) { return  legendre5x(x); }
static double legendre_der_6(double x) { return  legendre6x(x); }
static double legendre_der_7(double x) { return  legendre7x(x); }
static double legendre_der_8(double x) { return  legendre8x(x); }
static double legendre_der_9(double x) { return  legendre9x(x); }
static double legendre_der_10(double x) { return  legendre10x(x); }
static double legendre_der_11(double x) { return  legendre11x(x); }

shape_fn_t legendre_der_tab_1d[] = {
	legendre_der_0, legendre_der_1, legendre_der_2, legendre_der_3, legendre_der_4, legendre_der_5,
	legendre_der_6, legendre_der_7, legendre_der_8, legendre_der_9, legendre_der_10, legendre_der_11
};
