// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "lobatto.h"
#include "legendre.h"

int lobatto_order_1d[] = {
1,
 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
11,12,13,14,15,16,17,18,19,20,
21,22,23,24,25,26,27,28,29,30,
31,32,33,34,35,36,37,38,39,40,
41,42,43,44,45,46,47,48,49,50,
51,52,53,54,55,56,57,58,59,60,
61,62,63,64,65,66,67,68,69,70,
71,72,73,74,75,76,77,78,79,80,
81,82,83,84,85,86,87,88,89,90,
91,92,93,94,95,96,97,98,99,100,
};

static double lobatto_fn_0(double x) {
    return 1.0/2.0 - x/2;
}

static double lobatto_fn_1(double x) {
    return 1.0/2.0 + x/2;
}

static double lobatto_der_0(double x) {
    return -1.0/2.0;
}

static double lobatto_der_1(double x) {
    return 1.0/2.0;
}

// Fills an array of length MAX_P + 1 with Lobatto shape 
// functions (integrated normalized Legendre polynomials) 
// at point 'x'. 
extern void fill_lobatto_array(double x, 
                               double lobatto_array_val[MAX_P+1],
                               double lobatto_array_der[MAX_P+1]) {
    double legendre_array[MAX_P + 1];
    int max_fns_num = MAX_P + 1;
    // calculating (non-normalized) Legendre polynomials
    legendre_array[0] = 1.;
    legendre_array[1] = x;
    for (int i=1; i < max_fns_num; i++) {
      legendre_array[i+1]  = (2*i+1)*x*legendre_array[i] 
                             - i*legendre_array[i-1]; 
      legendre_array[i+1] /= i+1; 
    }
    // first fill the two linear Lobatto shape functions 
    lobatto_array_val[0] = lobatto_fn_0(x);
    lobatto_array_val[1] = lobatto_fn_1(x);
    lobatto_array_der[0] = lobatto_der_0(x);
    lobatto_array_der[1] = lobatto_der_1(x);
    // then fill the quadratic and higher which actually are 
    // the integrated Legendre polynomials
    for (int i=1; i < max_fns_num; i++) {
      lobatto_array_val[i+1] = 
        (legendre_array[i+1] - legendre_array[i-1]) / (2.*i + 1.);
      lobatto_array_val[i+1] /= leg_norm_const(i);
      lobatto_array_der[i+1] = legendre_array[i]; 
      lobatto_array_der[i+1] /= leg_norm_const(i);
    }
}

// FIXME - this function is inefficient, it fills the
// whole array
extern double calc_lobatto_val(double x, int n) 
{
    double val_array[MAX_P + 1];
    double der_array[MAX_P + 1];
    fill_lobatto_array(x, val_array, der_array);
    return val_array[n];
}

// FIXME - this function is inefficient, it fills the
// whole array
extern double calc_lobatto_der(double x, int n) 
{
    double val_array[MAX_P + 1];
    double der_array[MAX_P + 1];
    fill_lobatto_array(x, val_array, der_array);
    return der_array[n];
}

