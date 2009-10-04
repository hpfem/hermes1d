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

// evaluate approximate solution at element 'm' at reference 
// point 'x_ref'. Here 'y' is the global vector of coefficients
void Linearizer::eval_approx(Element *e, double x_ref, double *y, double &x_phys, double &val) {
  val = 0;
  for(int i=0; i <= e->p; i++) {
    if(e->dof[i] >= 0) val += y[e->dof[i]]*lobatto_fn_tab_1d[i](x_ref);
  }
  double a = e->v1->x;
  double b = e->v2->x;
  x_phys = (a+b)/2 + x_ref*(b-a)/2;
  return;
}

void Linearizer::plot_solution(const char *out_filename, double *y_prev, int plotting_elem_subdivision)
{
  // Plot solution in Gnuplot format
  Element *Elems = this->mesh->get_elems();
  FILE *f = fopen(out_filename, "wb");
  for(int m=0; m<this->mesh->get_n_elems(); m++) {
    double h = 2./plotting_elem_subdivision;
    double x_phys, val;
    for(int i=0; i<=plotting_elem_subdivision; i++) {
      double x_ref = -1. + i*h;
      eval_approx(Elems + m, x_ref, y_prev, x_phys, val);
      fprintf(f, "%g %g\n", x_phys, val);
    }
  }
  fclose(f);
}
