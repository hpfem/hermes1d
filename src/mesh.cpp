#include "mesh.h"

void Mesh::create(double A, double B, int n_elem)
{
  this->n_elem = n_elem;
  this->vertices = new Vertex[n_elem+1];    // allocate array of vertices
  double h = (B - A)/n_elem;
  for(int i = 0; i < n_elem+1; i++) {
    this->vertices[i].x = A + i*h;             // so far equidistant division only
  }
  this->elems = new Element[n_elem];                // allocate array of elements
  for(int i=0; i<n_elem; i++) {
    this->elems[i].p = -1;
    this->elems[i].v1 = this->vertices + i;
    this->elems[i].v2 = this->vertices + i + 1;
    this->elems[i].dof = NULL;
  }
}

void Mesh::set_poly_orders(int poly_order)
{
  for(int i=0; i < this->n_elem; i++) {
    this->elems[i].p = poly_order;
    this->elems[i].dof = new int[poly_order+1];
  }
}

void Mesh::assign_dofs()
{
  int DEBUG=1;
  // define element connectivities
  // (so far only for zero Dirichlet conditions)
  // (a) enumerate vertex dofs
  elems[0].dof[0] = -1;        // Dirichlet left end of domain
  elems[0].dof[1] = 0;         // first vertex dof
  for(int i=1; i<n_elem-1; i++) {
    elems[i].dof[0] = i-1;     // Dirichlet left end of domain
    elems[i].dof[1] = i;       // first dof
  }
  elems[n_elem-1].dof[0] = n_elem-2;     // last vertex dof
  elems[n_elem-1].dof[1] = -1;      // Dirichlet left end of domain
  // (b) enumerate bubble dofs
  n_dof = n_elem-1; 
  for(int i=0; i<n_elem; i++) {
    for(int j=2; j<=elems[i].p; j++) {
      elems[i].dof[j] = n_dof++;     // enumerating higher-order dofs
    }
  }

  // test (print DOF in elements)
  int n_test = 0;
  for(int i=0; i<n_elem; i++) n_test += elems[i].p;
  n_test -= 1;
  if(n_dof != n_test) {
    printf("n_dof = %d, n_test = %d\n", n_dof, n_test);
    error("Internal: Test of total DOF number failed."); 
  }

  // test (print element connectivities)
  if(DEBUG) {
    printf("Printing element DOF arrays:\n");
    printf("Elements = %d\n", n_elem);
    printf("DOF = %d", n_dof);
    for (int i = 0; i < n_elem; i++) {
      printf("\nElement[%d]: ", i); 
      for(int j = 0; j<elems[i].p+1; j++) {
        printf("%d, ", elems[i].dof[j]);
      }
    }
    printf("\n"); 
  }
}
