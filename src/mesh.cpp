// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "mesh.h"

void Mesh::create(double a, double b, int n_elem)
{
  this->n_elem = n_elem;
  this->vertices = new Vertex[this->n_elem+1]; // allocate array of vertices
  double h = (b - a)/this->n_elem;
  for(int i = 0; i < this->n_elem+1; i++) {
    this->vertices[i].x = a + i*h;             // so far equidistant division only
  }
  this->elems = new Element[this->n_elem];     // allocate array of elements
  for(int i=0; i<this->n_elem; i++) {
    this->elems[i].p = -1;
    this->elems[i].v1 = this->vertices + i;
    this->elems[i].v2 = this->vertices + i + 1;
    this->elems[i].dof = imalloc(n_eq);
  }
}

// sets uniform polynomial degrees in the mesh
// and allocates elememnt dof arrays
void Mesh::set_uniform_poly_order(int poly_order)
{
  for(int i=0; i < this->n_elem; i++) {
    this->elems[i].p = poly_order;
    // c is solution component
    for(int c=0; c<this->n_eq; c++) {
      this->elems[i].dof[c] = new int[poly_order+1];
    }
  }
}

int Mesh::assign_dofs()
{
  // define element connectivities
  // (so far only for zero Dirichlet conditions)
  // (a) enumerate vertex dofs
  int count = 0;
  for(int c=0; c<this->n_eq; c++) { // loop over solution components
    if (this->bc_left_dir[c])
        elems[0].dof[c][0] = -1;        // Dirichlet BC on the left
    else {
        elems[0].dof[c][0] = count;     // No Dirichlet BC on the left
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
        elems[this->n_elem-1].dof[c][1] = -1;    // Dirichlet BC on the right
    else {
        elems[this->n_elem-1].dof[c][1] = count; // No Dirichlet BC on the right
        count++;
    }
    // (b) enumerate bubble dofs
    for(int i=0; i<this->n_elem; i++) {
      for(int j=2; j<=elems[i].p; j++) {
        elems[i].dof[c][j] = count;              // enumerating higher-order dofs
        count++;
      }
    }
  }
  this->n_dof = count;

  // test (print element connectivities)
  if(0) {
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
    exit(0);
  }

  return this->n_dof;
}

// return coefficients for all shape functions on the element m,
// for all solution components
void Mesh::calculate_elem_coeffs(int m, double *y_prev, 
                                 double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM])
{
  for(int c=0; c<n_eq; c++) {
    if (m == 0 && elems[m].dof[c][0] == -1) {
        coeffs[c][0] = bc_left_dir_values[c];
    }
    else {
        coeffs[c][0] = y_prev[elems[m].dof[c][0]];
    }
    if (m == n_elem-1 && elems[m].dof[c][1] == -1) {
        coeffs[c][1] = bc_right_dir_values[c];
    }
    else {
        coeffs[c][1] = y_prev[elems[m].dof[c][1]];
    }
    for (int j=2; j<=elems[m].p; j++) {
        coeffs[c][j] = y_prev[elems[m].dof[c][j]];
    }
  }
}

// evaluate previous solution and its derivative 
// in the "pts_array" points
void Mesh::element_solution(Element *e, double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], int pts_num, 
		            double pts_array[MAX_PTS_NUM], double val[MAX_EQN_NUM][MAX_PTS_NUM], 
                            double der[MAX_EQN_NUM][MAX_PTS_NUM])
{
  double a = e->v1->x;
  double b = e->v2->x;
  double jac = (b-a)/2.; 
  int p = e->p;
  for(int c=0; c<n_eq; c++) { 
    for (int i=0 ; i<pts_num; i++) {
      der[c][i] = val[c][i] = 0;
      for(int j=0; j<=p; j++) {
        val[c][i] += coeff[c][j]*lobatto_fn_tab_1d[j](pts_array[i]);
        der[c][i] += coeff[c][j]*lobatto_der_tab_1d[j](pts_array[i]);
      }
      der[c][i] /= jac;
    }
  }
} 

// evaluate previous solution and its derivative 
// at the reference point x_ref to element 'e'.
void Mesh::element_solution_point(double x_ref, Element *e, 
			    double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], double *val, double *der)
{
  double a = e->v1->x;
  double b = e->v2->x;
  double jac = (b-a)/2.; 
  int p = e->p;
  for(int c=0; c<n_eq; c++) {
    der[c] = val[c] = 0;
    for(int j=0; j<=p; j++) {
      val[c] += coeff[c][j]*lobatto_fn_tab_1d[j](x_ref);
      der[c] += coeff[c][j]*lobatto_der_tab_1d[j](x_ref);
    }
    der[c] /= jac;
  }
} 

// transformation of k-th shape function defined on Gauss points 
// corresponding to 'order' to physical interval (a,b)
void Mesh::element_shapefn(double a, double b, 
		     int k, int order, double *val, double *der) {
  double2 *ref_tab = g_quad_1d_std.get_points(order);
  int pts_num = g_quad_1d_std.get_num_points(order);
  for (int i=0 ; i<pts_num; i++) {
    // change function values and derivatives to interval (a, b)
    val[i] = lobatto_fn_tab_1d[k](ref_tab[i][0]);
    double jac = (b-a)/2.; 
    der[i] = lobatto_der_tab_1d[k](ref_tab[i][0]) / jac; 
  }
};

// transformation of k-th shape function at the reference 
// point x_ref to physical interval (a,b).
void Mesh::element_shapefn_point(double x_ref, double a, double b, 
		                 int k, double *val, double *der) {
    // change function values and derivatives to interval (a, b)
    *val = lobatto_fn_tab_1d[k](x_ref);
    double jac = (b-a)/2.; 
    *der = lobatto_der_tab_1d[k](x_ref) / jac; 
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

// Tvaluate (vector-valued) approximate solution at reference 
// point 'x_ref' in element 'm'. Here 'y' is the global vector 
// of coefficients. The result is a vector of length mesh->n_eq
void Linearizer::eval_approx(Element *e, double x_ref, double *y, 
                             double *x_phys, double *val) {
  int n_eq = this->mesh->get_n_eq();
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

// Plot solution in Gnuplot format
void Linearizer::plot_solution(const char *out_filename, 
                               double *y_prev, int plotting_elem_subdivision)
{
  int n_eq = this->mesh->get_n_eq();
  int n_elem = this->mesh->get_n_elems();  
  Element *elems = this->mesh->get_elems();
  // FIXME:
  if(n_eq > MAX_EQN_NUM)
      error("number of equations too high in plot_solution().");
  FILE *f[MAX_EQN_NUM];
  char final_filename[MAX_EQN_NUM][MAX_STRING_LENGTH];
  for(int c=0; c<n_eq; c++) {
    if(n_eq == 1) sprintf(final_filename[c], "%s", out_filename);
    else sprintf(final_filename[c], "%s_%d", out_filename, c);
    f[c] = fopen(final_filename[c], "wb");
    if(f[c] == NULL) error("problem opening file in plot_solution().");
  }
  // FIXME
  if(plotting_elem_subdivision > 100) 
  error ("plotting_elem_subdivision too high in plot_solution()."); 
  double phys_u_prev[MAX_EQN_NUM][MAX_PTS_NUM];
  double phys_du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM];
  for(int m=0; m<n_elem; m++) {
    // FIXME:
    if(elems[m].p > 90) error ("element degree too hign in plot(solution)."); 
    double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
    this->mesh->calculate_elem_coeffs(m, y_prev, coeffs); 

    double pts_array[MAX_PTS_NUM];
    double h = 2./plotting_elem_subdivision;

    for (int j=0; j<plotting_elem_subdivision+1; j++) pts_array[j] = -1 + j*h;
    this->mesh->element_solution(elems + m, coeffs, plotting_elem_subdivision+1, 
                       pts_array, phys_u_prev, phys_du_prevdx); 
    double a = elems[m].v1->x;
    double b = elems[m].v2->x;
    // loop over solution components
    for(int c=0; c<n_eq; c++) {
      for (int j=0; j<plotting_elem_subdivision+1; j++) {
        fprintf(f[c], "%g %g\n", (a + b)/2  
                + pts_array[j] * (b-a)/2, phys_u_prev[c][j]);
      }
    }
  }
  for(int c=0; c<n_eq; c++) {
    printf("Output written to %s.\n", final_filename[c]);
    fclose(f[c]);
  }
}

// Returns pointers to x and y coordinates in **x and **y
// you should free it yourself when you don't need it anymore
// y_prev --- the solution coefficients (all equations)
// comp --- which component you want to process
// plotting_elem_subdivision --- the number of subdivision of the element
// x, y --- the doubles list of x,y
// n --- the number of points

void Linearizer::get_xy(double *y_prev, int comp,
        int plotting_elem_subdivision,
        double **x, double **y, int *n)
{
    int n_eq = this->mesh->get_n_eq();
    int n_elem = this->mesh->get_n_elems();
    Element *elems = this->mesh->get_elems();

    *n = n_elem * (plotting_elem_subdivision+1);
    double *x_out = new double[*n];
    double *y_out = new double[*n];

    // FIXME:
    if(n_eq > MAX_EQN_NUM)
        error("number of equations too high in plot_solution().");
    // FIXME
    if(plotting_elem_subdivision > MAX_PTS_NUM)
        error("plotting_elem_subdivision too high in plot_solution().");
    double phys_u_prev[MAX_EQN_NUM][MAX_PTS_NUM];
    double phys_du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM];
    for(int m=0; m<n_elem; m++) {
        // FIXME:
        if(elems[m].p > 90)
            error("element degree too hign in plot(solution).");
        double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
        this->mesh->calculate_elem_coeffs(m, y_prev, coeffs);

        double pts_array[MAX_PTS_NUM];
        double h = 2./plotting_elem_subdivision;

        for (int j=0; j<plotting_elem_subdivision+1; j++)
            pts_array[j] = -1 + j*h;
        this->mesh->element_solution(elems + m, coeffs,
                plotting_elem_subdivision+1, pts_array,
                phys_u_prev, phys_du_prevdx);
        double a = elems[m].v1->x;
        double b = elems[m].v2->x;
        for (int j=0; j<plotting_elem_subdivision+1; j++) {
            x_out[m*(plotting_elem_subdivision+1) + j] =
                (a + b)/2 + pts_array[j] * (b-a)/2;
            y_out[m*(plotting_elem_subdivision+1) + j] =
                phys_u_prev[comp][j];
        }
    }
    *x = x_out;
    *y = y_out;
}
