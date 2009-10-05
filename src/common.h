#ifndef __HERMES1D_COMMON_H
#define __HERMES1D_COMMON_H

#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>

#define BOUNDARY_LEFT 0
#define BOUNDARY_RIGHT 1

void error(const char *msg);

typedef double scalar;

typedef double double2[2];

typedef int int2[2];

typedef double (*shape_fn_t)(double);

struct Vertex {
  double x;
};

class Element {
public:
  Vertex *v1, *v2; // endpoints
  int p;           // poly degree
  int *dof;        // connectivity array of length p+1
};

void intro();

void throw_exception(char *text);

#define MEM_CHECK(var) if (var == NULL) { printf("Out of memory."); exit(1); }

#define verbose(msg)
#define warn(msg)

#endif
