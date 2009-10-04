#ifndef __HERMES1D_COMMON_H
#define __HERMES1D_COMMON_H

#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>

void error(const char *msg);

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

#endif
