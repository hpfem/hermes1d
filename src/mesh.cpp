// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "mesh.h"
#include "iterator.h"

Element::Element() 
{
  x1 = x2 = 0;
  p = 0; 
  dof = NULL;
  //sons = NULL; 
  active = 1;
  level = 0;
  id = -1;
}

Element::Element(double x_left, double x_right, int deg, int n_eq) 
{
  x1 = x_left;
  x2 = x_right;
  p = deg; 
  this->dof_alloc(n_eq);
  if (dof == NULL) error("Not enough memory in Element().");
  //sons = NULL; 
  active = 1;
  level = 0;
  id = -1;
}

unsigned Element::is_active() 
{
  return this->active;
}

void Element::refine(int p_left, int p_right, int n_eq) 
{
  double x1 = this->x1;
  double x2 = this->x2;
  double midpoint = (x1 + x2)/2.; 
  this->sons[0] = new Element(x1, midpoint, p_left, n_eq);
  this->sons[1] = new Element(midpoint, x2, p_right, n_eq);
  this->sons[0]->level = this->level + 1; 
  this->sons[1]->level = this->level + 1; 
  // copying negative dof to sons if any
  for(int c=0; c<n_eq; c++) {
    if (this->dof[c][0] < 0) this->sons[0]->dof[c][0] = this->dof[c][0];
    if (this->dof[c][1] < 0) this->sons[1]->dof[c][1] = this->dof[c][1];
  }
  this->active = 0;
}

void Element::dof_alloc(int n_eq) 
{
  this->dof = (int**)malloc(n_eq*sizeof(int*));
  if(this->dof == NULL) error("Element dof_alloc() failed.");
  // c is solution component
  for(int c=0; c<n_eq; c++) {
    this->dof[c] = new int[MAX_POLYORDER + 1];
    // important for th etreatment of boundary conditions
    for(int i=0; i<MAX_POLYORDER + 1; i++) this->dof[c][i] = 0;
  }
}

Mesh::Mesh(double a, double b, int n_base_elem, int p_init, int n_eq)
{
  // domain end points
  left_endpoint = a;
  right_endpoint = b;
  // number of equations
  this->n_eq = n_eq;
  // print the banner (only once)
  static int n_calls = 0;
  n_calls++;
  if (n_calls == 1) intro();
  // check maximum number of equations
  if(n_eq > MAX_EQN_NUM) 
  error("Maximum number of equations exceeded (set in common.h)");
  // arrays for boundary conditions
  this->bc_left_dir_values = new double[n_eq];
  this->bc_right_dir_values = new double[n_eq];
  for (int i=0; i<n_eq; i++) {
    this->bc_left_dir_values[i] = 0;
    this->bc_right_dir_values[i] = 0;
  }
  // number of elements
  this->n_base_elem = n_base_elem;
  this->n_active_elem = n_base_elem;
  // allocate element array
  this->base_elems = new Element[this->n_base_elem];     
  if (base_elems == NULL) error("Not enough memory in Mesh::create().");
  if (p_init > MAX_POLYORDER) 
    error("Max element order exceeded (set in common.h).");
  // element length
  double h = (b - a)/this->n_base_elem;          
  // fill initial element array
  for(int i=0; i<this->n_base_elem; i++) {         
    // polynomial degree
    this->base_elems[i].p = p_init;
    // allocate element dof arrays
    this->base_elems[i].dof_alloc(n_eq);
    // define element end points
    this->base_elems[i].x1 = a + i*h;
    this->base_elems[i].x2 = this->base_elems[i].x1 + h;
  }
  this->assign_elem_ids();
}

int Mesh::get_n_active_elems()
{
    return this->n_active_elem;
}

void Mesh::refine_single_elem(int id, int p_left, int p_right)
{
    Iterator I(this);
    Element *e;
    while ((e = I.next_active_element()) != NULL) {
        printf("%d %d\n", e->id, id);
        if (e->id == id) {
            e->refine(p_left, p_right, this->n_eq);
            this->n_active_elem++;
            return;
        }
    }
    error("refine_single_elem: Element not found.");
}

void Mesh::refine_multi_elems(int n, int *id_array, int2 *p_pair_array)
{
    Iterator *I = new Iterator(this);
    Element *e;
    int count = 0;
    while ((e = I->next_active_element()) != NULL) {
        if (e->id == id_array[count]) {
            if (count >= n)
                error("refine_multi_elems: not enough elems specified");
            e->refine(p_pair_array[count][0], p_pair_array[count][1], this->n_eq);
            this->n_active_elem++;
            count++;
        }
    }
}

void Mesh::set_bc_left_dirichlet(int eq_n, double val)
{
  this->bc_left_dir_values[eq_n] = val;
  // deactivate the corresponding dof for the left-most
  // element and all his descendants adjacent to the 
  // left boundary
  Element *e = &(this->base_elems[0]);
  do {
    e->dof[eq_n][0] = -1;
    e = e->sons[0];
  } while (e != NULL);
}

void Mesh::set_bc_right_dirichlet(int eq_n, double val)
{
  this->bc_right_dir_values[eq_n] = val;
  // deactivate the corresponding dof for the right-most
  // element and all his descendants adjacent to the 
  // right boundary
  Element *e = &(this->base_elems[this->n_base_elem-1]);
  do {
    e->dof[eq_n][1] = -1;
    e = e->sons[1];
  } while (e != NULL);

}

// define element connectivities (global dof)
int Mesh::assign_dofs()
{
  Iterator *I = new Iterator(this);
  // (1) enumerate vertex dofs
  int count_dof = 0;
  // loop over solution components
  for(int c=0; c<this->n_eq; c++) {    
    Element *e;
    I->reset();
    while ((e = I->next_active_element()) != NULL) {
      if (e->dof[c][0] != -1) e->dof[c][0] = count_dof++; 
      if (e->dof[c][1] != -1) e->dof[c][1] = count_dof; 
      else count_dof--;
    }
    count_dof++;
    // (2) enumerate bubble dofs
    I->reset();
    while ((e = I->next_active_element()) != NULL) {
      for(int j=2; j<=e->p; j++) {
        e->dof[c][j] = count_dof;
        count_dof++;
      }
    }
  }
  this->n_dof = count_dof;

  // enumerate elements
  this->assign_elem_ids();

  // print element connectivities
  if(DEBUG_ELEM_DOF) {
    printf("Printing element DOF arrays:\n");
    printf("Elements = %d\n", this->n_base_elem);
    printf("DOF = %d", this->n_dof);
    for(int c = 0; c<this->n_eq; c++) {
      I->reset();
      Element *e;
      while ((e = I->next_active_element()) != NULL) {
        printf("\nElement (%g, %g), id = %d\n ", e->x1, e->x2, e->id); 
        for(int j = 0; j<e->p + 1; j++) {
          printf("dof[%d][%d] = %d\n ", c, j, e->dof[c][j]);
        }
      }
    }
    printf("\n"); 
  }

  return this->n_dof;
}

int Mesh::assign_elem_ids()
{
    Iterator *I = new Iterator(this);
    int count_id = 0;
    Element *e;
    I->reset();
    while ((e = I->next_active_element()) != NULL) {
        e->id = count_id++;
    }
}

// return coefficients for all shape functions on the element m,
// for all solution components
void Mesh::calculate_elem_coeffs(Element *e, double *y_prev, 
                                 double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM])
{
  if (!e->is_active()) error("Internal in calculate_elem_coeffs().");
  for(int c=0; c<n_eq; c++) {
    // coeff of the left vertex function
    if (e->dof[c][0] == -1) coeffs[c][0] = bc_left_dir_values[c];
    else coeffs[c][0] = y_prev[e->dof[c][0]];
    // coeff of the right vertex function
    if (e->dof[c][1] == -1) coeffs[c][1] = bc_right_dir_values[c];
    else coeffs[c][1] = y_prev[e->dof[c][1]];
    //completing coeffs of bubble functions
    for (int j=2; j<=e->p; j++) {
        coeffs[c][j] = y_prev[e->dof[c][j]];
    }
  }
}

// evaluate previous solution and its derivative 
// in the "pts_array" points
void Mesh::element_solution(Element *e, 
                            double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], 
                            int pts_num, double pts_array[MAX_PTS_NUM], 
                            double val[MAX_EQN_NUM][MAX_PTS_NUM], 
                            double der[MAX_EQN_NUM][MAX_PTS_NUM])
{
  double a = e->x1;
  double b = e->x2;
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
			    double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], 
                            double *val, double *der)
{
  double a = e->x1;
  double b = e->x2;
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

// Evaluate (vector-valued) approximate solution at reference 
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
  double a = e->x1;
  double b = e->x2;
  *x_phys = (a+b)/2 + x_ref*(b-a)/2;
  return;
}

// Plot solution in Gnuplot format
void Linearizer::plot_solution(const char *out_filename, 
                               double *y_prev, int plotting_elem_subdivision)
{
    int n_eq = this->mesh->get_n_eq();
    FILE *f[MAX_EQN_NUM];
    char final_filename[MAX_EQN_NUM][MAX_STRING_LENGTH];
    for(int c=0; c<n_eq; c++) {
        if(n_eq == 1)
            sprintf(final_filename[c], "%s", out_filename);
        else
            sprintf(final_filename[c], "%s_%d", out_filename, c);
        f[c] = fopen(final_filename[c], "wb");
        if(f[c] == NULL) error("problem opening file in plot_solution().");
        int n;
        double *x, *y;
        this->get_xy(y_prev, c, plotting_elem_subdivision, &x, &y, &n);
        for (int i=0; i < n; i++)
            fprintf(f[c], "%g %g\n", x[i], y[i]);
        delete[] x;
        delete[] y;
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
    int n_active_elem = this->mesh->get_n_active_elems();
    Iterator *I = new Iterator(this->mesh);

    *n = n_active_elem * (plotting_elem_subdivision+1);
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
        
    Element *e;
    int counter = 0;
    while ((e = I->next_active_element()) != NULL) {
        if (counter >= n_active_elem)
            error("Internal error: wrong n_active_elem");
        // FIXME:
        if(e->p > MAX_POLYORDER)
            error("element degree too hign in plot(solution).");
        double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
        this->mesh->calculate_elem_coeffs(e, y_prev, coeffs);

        double pts_array[MAX_PTS_NUM];
        double h = 2./plotting_elem_subdivision;

        for (int j=0; j<plotting_elem_subdivision+1; j++)
            pts_array[j] = -1 + j*h;
        this->mesh->element_solution(e, coeffs,
                plotting_elem_subdivision+1, pts_array,
                phys_u_prev, phys_du_prevdx);
        double a = e->x1;
        double b = e->x2;
        for (int j=0; j<plotting_elem_subdivision+1; j++) {
            x_out[counter*(plotting_elem_subdivision+1) + j] =
                (a + b)/2 + pts_array[j] * (b-a)/2;
            y_out[counter*(plotting_elem_subdivision+1) + j] =
                phys_u_prev[comp][j];
        }
        counter++;
    }
    *x = x_out;
    *y = y_out;
}
