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
