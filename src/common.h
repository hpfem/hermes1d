#include <stdlib.h>

void error(const char *msg) {
  printf("Error: %s\n", msg);
  exit(0);
};

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

