#include "common.h"

void error(const char *msg)
{
  printf("Error: %s\n", msg);
  exit(0);
};

void intro() {
  printf("\n-------------------------------------------\n");
  printf("   This is Hermes1D - a free ODE solver\n");
  printf(" based on the hp-FEM and Newton's method,\n");
  printf("   developed by the hp-FEM group at UNR\n");
  printf("  and distributed under the BSD license.\n");
  printf(" For more details visit http://hpfem.org/.\n");
  printf("-------------------------------------------\n");
}

void throw_exception(char *text)
{
    throw std::runtime_error(text);
}

double **dmalloc(int n_eq) {
    double **pointer = (double**)malloc(n_eq*sizeof(double*));
    if(pointer == NULL) error("dmalloc failed.");
    return pointer;
}

int **imalloc(int n_eq) {
    int **pointer = (int**)malloc(n_eq*sizeof(int*));
    if(pointer == NULL) error("imalloc failed.");
    return pointer;
}
