#include "mesh.h"

void Mesh::create(double a, double a, int n_elem)
{
  this->n_elem = n_elem;
  this->vertices = new Vertex[this->n_elem+1];    // allocate array of vertices
  double h = (b - a)/this->n_elem;
  for(int i = 0; i < this->n_elem+1; i++) {
    this->vertices[i].x = a + i*h;             // so far equidistant division only
  }
  this->elems = new Element[this->n_elem];                // allocate array of elements
  for(int i=0; i<this->n_elem; i++) {
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
  // define element connectivities
  // (so far only for zero Dirichlet conditions)
  // (a) enumerate vertex dofs
  int count = 0;
  for(int c=0; c<this->n_eq; c++) { // loop over solution components
    if (this->bc_left_dir[c])
        elems[0].dof[c][0] = -1;        // Dirichlet BC on the left
    else {
        elems[0].dof[c][0] = count;        // No Dirichlet BC on the left
        count++;
    }
    elems[0].dof[c][1] = count;         // first vertex dof
    for(int i=1; i<this->n_elem-1; i++) {
      elems[i].dof[c][0] = count;
      count++;
      elems[i].dof[c][1] = count;
    }
    elems[this->n_elem-1].dof[c][0] = count;
    count++;
    if (this->bc_right_dir[c])
        elems[this->n_elem-1].dof[c][1] = -1;      // Dirichlet BC on the right
    else {
        elems[this->n_elem-1].dof[c][1] = count;        // No Dirichlet BC on the right
        count++;
    }
    // (b) enumerate bubble dofs
    for(int i=0; i<this->n_elem; i++) {
      for(int j=2; j<=elems[i].p; j++) {
        elems[i].dof[c][j] = count;     // enumerating higher-order dofs
        count++;
      }
    }
  }
  this->n_dof = count;

  // test (print element connectivities)
  if(DEBUG) {
    printf("Printing element DOF arrays:\n");
    printf("Elements = %d\n", this->n_elem);
    printf("DOF = %d", this->n_dof);
    for (int i = 0; i < this->n_elem; i++) {
      printf("\nElement[%d]:\n ", i); 
      for(int c = 0; c<this->n_eq; c++) {
        for(int j = 0; j<elems[i].p+1; j++) {
          printf("dof[%d][%d] = %d\n ", c, j, elems[i].dof[c][j]);
        }
      }
    }
    printf("\n"); 
  }
}

// evaluate approximate solution at element 'm' at reference 
// point 'x_ref'. Here 'y' is the global vector of coefficients
void Linearizer::eval_approx(Element *e, double x_ref, double *y, double *x_phys, double *val) {
  int n_eq = this->mesh->n_eq;
  for(int c=0; c<n_eq; c++) { // loop over solution components
    val[c] = 0;
    for(int i=0; i <= e->p; i++) { // loop over shape functions
      if(e->dof[c][i] >= 0) val[c] += y[e->dof[c][i]]*lobatto_fn_tab_1d[i](x_ref);
    }
  }
  double a = e->v1->x;
  double b = e->v2->x;
  *x_phys = (a+b)/2 + x_ref*(b-a)/2;
  return;
}

// c is solution component
void calculate_elem_coeffs(Mesh *mesh, int m, double *y_prev, double **coeffs, int c);

void Linearizer::plot_solution(const char *out_filename, double *y_prev, int plotting_elem_subdivision)
{
  int n_eq = this->mesh->n_eq;
  // Plot solution in Gnuplot format
  Element *elems = this->mesh->get_elems();
  FILE *f = fopen(out_filename, "wb");
  // FIXME: this is a memory leak!!!
  double *phys_u_prev =    new double[plotting_elem_subdivision + 1];
  double *phys_du_prevdx = new double[plotting_elem_subdivision + 1];
  for(int m=0; m<this->mesh->get_n_elems(); m++) {
    double coeffs[10][100];   // FIXME
    calculate_elem_coeffs(this->mesh, m, y_prev, coeffs, c); 

    double pts_array[plotting_elem_subdivision+1];
    double h = 2./plotting_elem_subdivision;

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

void Mesh::set_bc_left_dirichlet(int eq_n, double val)
{
    this->bc_left_dir[eq_n] = 1;
    this->bc_left_dir_values[eq_n] = val;
}

void Mesh::set_bc_right_dirichlet(int eq_n, double val)
{
    this->bc_right_dir[eq_n] = 1;
    this->bc_right_dir_values[eq_n] = val;
}

void Mesh::set_bc_left_natural(int eq_n)
{
    this->bc_left_dir[eq_n] = 0;
}

void Mesh::set_bc_right_natural(int eq_n)
{
    this->bc_right_dir[eq_n] = 0;
}

