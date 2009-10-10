#ifndef __HERMES1D_COMMON_H
#define __HERMES1D_COMMON_H

#define DEBUG 0

#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>

#define BOUNDARY_LEFT 0
#define BOUNDARY_RIGHT 1

#define BC_NATURAL 0
#define BC_DIRICHLET 1

const int MAX_EQN_NUM = 10;            // maximum number of equations in the system
const int MAX_PTS_NUM = 101;           // refers both to plotting subdivision and Gauss quadrature points
const int MAX_P = 50;                  // this is the maximum polynomial degree allowed in elements
const int MAX_COEFFS_NUM = MAX_P + 1;  // this is the maximum polynomial degree allowed in elements
const int MAX_STRING_LENGTH = 100;     // maximum string length 

void error(const char *msg);

typedef double scalar;

typedef double double2[2];

typedef int int2[2];

typedef double (*shape_fn_t)(double);

void intro();

void throw_exception(char *text);

#define MEM_CHECK(var) if (var == NULL) { printf("Out of memory."); exit(1); }

#define verbose(msg)
#define warn(msg)

double **dmalloc(int n_eq);
   
int **imalloc(int n_eq);


#endif
