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
  int count = 0;
  if (this->dir_bc_left_active[0])
      elems[0].dof[0] = -1;        // Dirichlet BC on the left
  else {
      elems[0].dof[0] = count;        // No Dirichlet BC on the left
      count++;
  }
  elems[0].dof[1] = count;         // first vertex dof
  for(int i=1; i<n_elem-1; i++) {
    elems[i].dof[0] = count;
    count++;
    elems[i].dof[1] = count;
  }
  elems[n_elem-1].dof[0] = count;
  count++;
  if (this->dir_bc_right_active[0])
      elems[n_elem-1].dof[1] = -1;      // Dirichlet BC on the right
  else {
      elems[n_elem-1].dof[1] = count;        // No Dirichlet BC on the right
      count++;
  }
  // (b) enumerate bubble dofs
  for(int i=0; i<n_elem; i++) {
    for(int j=2; j<=elems[i].p; j++) {
      elems[i].dof[j] = count;     // enumerating higher-order dofs
      count++;
    }
  }
  n_dof = count;

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
  Element *elems = this->mesh->get_elems();
  FILE *f = fopen(out_filename, "wb");
  // FIXME: this is a memory leak!!!
  double *phys_u_prev =    new double[plotting_elem_subdivision + 1];
  double *phys_du_prevdx = new double[plotting_elem_subdivision + 1];
  for(int m=0; m<this->mesh->get_n_elems(); m++) {
    double coeffs[100];
    if (m == 0 && elems[m].dof[0] == -1)
        coeffs[0] = this->mesh->dir_bc_left_values[0];
    else
        coeffs[0] = y_prev[elems[m].dof[0]];
    if (m == this->mesh->get_n_elems()-1 && elems[m].dof[1] == -1)
        coeffs[1] = this->mesh->dir_bc_right_values[0];
    else
        coeffs[1] = y_prev[elems[m].dof[1]];
    for (int j=2; j<=elems[m].p; j++)
        coeffs[j] = y_prev[elems[m].dof[j]];
    double pts_array[plotting_elem_subdivision+1];
    double h = 2./plotting_elem_subdivision;
    //double h = (elems[m].v2->x - elems[m].v1->x)/plotting_elem_subdivision;
    for (int j=0; j<plotting_elem_subdivision+1; j++)
        pts_array[j] = -1 + j*h;
    element_solution(elems + m, coeffs, plotting_elem_subdivision+1, pts_array, phys_u_prev, phys_du_prevdx); 
    double a = elems[m].v1->x;
    double b = elems[m].v2->x;
    for (int j=0; j<plotting_elem_subdivision+1; j++)
      fprintf(f, "%g %g\n", (a + b)/2 + pts_array[j] * (b-a)/2, phys_u_prev[j]);
  }
  fclose(f);
}

void Mesh::set_dirichlet_bc_left(int eq_n, double val)
{
    this->dir_bc_left_active[eq_n] = 1;
    this->dir_bc_left_values[eq_n] = val;
}

void Mesh::set_dirichlet_bc_right(int eq_n, double val)
{
    this->dir_bc_right_active[eq_n] = 1;
    this->dir_bc_right_values[eq_n] = val;
}
