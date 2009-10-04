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

