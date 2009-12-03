// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES1D_COMMON_H
#define __HERMES1D_COMMON_H

// print general debug information
#define DEBUG 0

// printf debug information about the stiffness/Jacobi matrix
#define DEBUG_MATRIX 0

#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>

#define BOUNDARY_LEFT 0
#define BOUNDARY_RIGHT 1

const int MAX_P = 100;                 // this is the max poly degree allowed in elements
                                       // WARNING: projections taking place in 
                                       // transfer_solution()
const int MAX_NEWTON_ITER_NUM = 100;   // maximum allowed number of Newton's iterations

const int MAX_CAND_NUM = 100;          // maximum allowed number of hp-refinement
                                       // candidates of an element

const int MAX_ELEM_NUM = 100000;       // maximum number of elements
const int MAX_EQN_NUM = 10;            // maximum number of equations in the system
const int MAX_PTS_NUM = 501;           // refers both to plotting subdivision 
                                       // and Gauss quadrature points
const int MAX_COEFFS_NUM = MAX_P + 1;  // this is the maximum number of 
                                       // polynomial coefficients
const int MAX_STRING_LENGTH = 100;     // maximum string length 

typedef double (*exact_sol_type)(double x, 
        double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]);

void error(const char *msg);
void error(const char *msg, const char *msg1);
void info(const char *msg, const char *msg1);
void warning(const char *msg);

typedef double scalar;

typedef double double2[2];

typedef int int2[2];
typedef int int3[3];

typedef double (*shape_fn_t)(double);

void intro();

void throw_exception(char *text);

#define MEM_CHECK(var) if (var == NULL) { printf("Out of memory."); exit(1); }

#define verbose(msg)
#define warn(msg)

double max(double a, double b);

#endif
