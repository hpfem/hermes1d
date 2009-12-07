// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "mesh.h"
#include "iterator.h"
#include "adapt.h"
#include "transforms.h"
#include "linearizer.h"

// debug - prints element errors as they come to adapt()
int PRINT_ELEM_ERRORS = 1;

// debug - prints element dof arrays in assign_dofs()
int DEBUG_ELEM_DOF = 0;

Element::Element() 
{
  x1 = x2 = 0;
  p = 0; 
  dof = NULL;
  sons[0] = sons[1] = NULL; 
  active = 1;
  level = 0;
  id = -1;
  dof_size = 0;
}

Element::Element(double x_left, double x_right, int lev, int deg, int n_eq) 
{
  x1 = x_left;
  x2 = x_right;
  p = deg; 
  dof_size = n_eq;
  this->dof_alloc();
  if (dof == NULL) error("Not enough memory in Element().");
  sons[0] = sons[1] = NULL; 
  active = 1;
  level = lev;
  id = -1;
}

unsigned Element::is_active() 
{
  return this->active;
}

void Element::refine(int type, int p_left, int p_right) 
{
  if(type == 0) {         // p-refinement
    this->p = p_left;
  }
  else {
    double x1 = this->x1;
    double x2 = this->x2;
    double midpoint = (x1 + x2)/2.; 
    this->sons[0] = new Element(x1, midpoint, this->level + 1, p_left, dof_size);
    this->sons[1] = new Element(midpoint, x2, this->level + 1, p_right, dof_size);
    // copying negative dof to sons if any
    for(int c=0; c<dof_size; c++) {
      if (this->dof[c][0] < 0) this->sons[0]->dof[c][0] = this->dof[c][0];
      if (this->dof[c][1] < 0) this->sons[1]->dof[c][1] = this->dof[c][1];
    }
    this->active = 0;
  }
}

void Element::refine(int3 cand) 
{
  this->refine(cand[0], cand[1], cand[2]);
}

// initialize element and allocate dof arrays for all
// solution components
void Element::init(double x1, double x2, int p_init, 
                   int id, int active, int level, int n_eq)
{
  this->x1 = x1;
  this->x2 = x2;
  this->p = p_init;
  this->id = id;
  this->active = active;
  this->level = level;
  this->dof_size = n_eq;
  this->dof_alloc();
}

void Element::dof_alloc() 
{
  this->dof = (int**)malloc(dof_size*sizeof(int*));
  if(this->dof == NULL) error("Element dof_alloc() failed.");
  // c is solution component
  for(int c=0; c<dof_size; c++) {
    this->dof[c] = new int[MAX_P + 1];
    // important for th etreatment of boundary conditions
    for(int i=0; i<MAX_P + 1; i++) this->dof[c][i] = 0;
  }
}

// return coefficients for all shape functions on the element m,
// for all solution components
void Element::get_coeffs(double *y_prev, 
                         double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM],
                         double bc_left_dir_values[MAX_EQN_NUM],
                         double bc_right_dir_values[MAX_EQN_NUM])
{
  if (!this->is_active()) error("Internal in calculate_elem_coeffs().");
  int dof_size = this->dof_size;
  for(int c=0; c<dof_size; c++) {
    // coeff of the left vertex function
    if (this->dof[c][0] == -1) coeffs[c][0] = bc_left_dir_values[c];
    else coeffs[c][0] = y_prev[this->dof[c][0]];
    // coeff of the right vertex function
    if (this->dof[c][1] == -1) coeffs[c][1] = bc_right_dir_values[c];
    else coeffs[c][1] = y_prev[this->dof[c][1]];
    //completing coeffs of bubble functions
    for (int j=2; j<=this->p; j++) {
        coeffs[c][j] = y_prev[this->dof[c][j]];
    }
  }
}

// Evaluate solution and its derivatives in quadrature points 'x_phys' 
// in the element (coeffs[][] array provided).
void Element::get_solution_quad(int flag, double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], 
                                int quad_order, 
                                double val_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
				double der_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM])
{
  double phys_x[MAX_QUAD_PTS_NUM];          // quad points
  double phys_w[MAX_QUAD_PTS_NUM];          // quad weights
  int pts_num;
  create_phys_element_quadrature(this->x1, this->x2, 
                                 quad_order, phys_x, phys_w, 
                                 &pts_num); 

  double jac = (this->x2 - this->x1)/2.; // Jacobian of reference map
  int dof_size = this->dof_size;
  int p = this->p;
  double x_ref[MAX_QUAD_PTS_NUM];
  // filling the values and derivatives
  if (flag == 0) { // integration points in the whole element
    for(int c=0; c<dof_size; c++) { 
      for (int i=0 ; i < pts_num; i++) {
        der_phys[c][i] = val_phys[c][i] = 0;
        for(int j=0; j<=p; j++) {
          val_phys[c][i] += coeff[c][j]*lobatto_val_ref_tab[quad_order][i][j];
          der_phys[c][i] += coeff[c][j]*lobatto_der_ref_tab[quad_order][i][j];
        }
        der_phys[c][i] /= jac;
      }
    }
  }
  if (flag == -1) { // integration points in the left half of element
    for(int c=0; c<dof_size; c++) { 
      for (int i=0 ; i < pts_num; i++) {
        der_phys[c][i] = val_phys[c][i] = 0;
        for(int j=0; j<=p; j++) {
          val_phys[c][i] += coeff[c][j]*lobatto_val_ref_tab_left[quad_order][i][j];
          der_phys[c][i] += coeff[c][j]*lobatto_der_ref_tab_left[quad_order][i][j];
        }
        der_phys[c][i] /= jac;
      }
    }
  }
  if (flag == 1) { // integration points in the right half of element
    for(int c=0; c<dof_size; c++) { 
      for (int i=0 ; i < pts_num; i++) {
        der_phys[c][i] = val_phys[c][i] = 0;
        for(int j=0; j<=p; j++) {
          val_phys[c][i] += coeff[c][j]*lobatto_val_ref_tab_right[quad_order][i][j];
          der_phys[c][i] += coeff[c][j]*lobatto_der_ref_tab_right[quad_order][i][j];
        }
        der_phys[c][i] /= jac;
      }
    }
  }
} 

// Evaluate solution and its derivatives in plotting points 'x_phys' 
// in the element (coeffs[][] array provided).
void Element::get_solution_plot(double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], 
                           int pts_num, double x_phys[MAX_PLOT_PTS_NUM], 
                           double val_phys[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
                           double der_phys[MAX_EQN_NUM][MAX_PLOT_PTS_NUM])
{
  double x1 = this->x1;
  double x2 = this->x2;
  double jac = (x2-x1)/2.; // Jacobian of reference map
  int dof_size = this->dof_size;
  int p = this->p;
  double x_ref[MAX_PLOT_PTS_NUM];
  // transforming points to (-1, 1)
  for (int i=0 ; i < pts_num; i++) x_ref[i] = inverse_map(x1, x2, x_phys[i]);
  // filling the values and derivatives
  for(int c=0; c<dof_size; c++) { 
    for (int i=0 ; i < pts_num; i++) {
      der_phys[c][i] = val_phys[c][i] = 0;
      for(int j=0; j<=p; j++) {
        val_phys[c][i] += coeff[c][j]*lobatto_val_ref(x_ref[i], j);
        der_phys[c][i] += coeff[c][j]*lobatto_der_ref(x_ref[i], j);
      }
      der_phys[c][i] /= jac;
    }
  }
} 

// Evaluate solution and its derivatives in quadrature points 'x_phys' 
// in the element (coeffs[][] array not provided).
void Element::get_solution_quad(int flag, int quad_order, double *y_prev, 
                                double val_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                                double der_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                                double bc_left_dir_values[MAX_EQN_NUM],
                                double bc_right_dir_values[MAX_EQN_NUM])
{

  double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM];
  this->get_coeffs(y_prev, coeff, bc_left_dir_values, bc_right_dir_values);
  this->get_solution_quad(flag, coeff, quad_order, val_phys, der_phys);
} 

// Calculates square of the L2 or H1 norm of the solution in element
double Element::calc_elem_norm_squared(int norm, double *y_prev, 
                                double bc_left_dir_values[MAX_EQN_NUM],
                                double bc_right_dir_values[MAX_EQN_NUM])
{
  double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM];
  this->get_coeffs(y_prev, coeff, bc_left_dir_values, bc_right_dir_values);
  double val_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM];
  double der_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM];
  int pts_num;
  int order = 2*this->p;
  this->get_solution_quad(0, coeff, order, val_phys, der_phys);

  double phys_x[MAX_QUAD_PTS_NUM];
  double phys_weights[MAX_QUAD_PTS_NUM];
  create_phys_element_quadrature(this->x1, this->x2, order, phys_x, phys_weights,
                                 &pts_num); 

  // integrate square over (-1, 1)
  int n_eq = this->dof_size;
  double norm_squared[MAX_EQN_NUM];
  for (int c=0; c<n_eq; c++) {
    norm_squared[c] = 0;
    for (int i=0; i<pts_num; i++) {
      double val = val_phys[c][i];
      if (norm == 1) {
        double der;
        der = der_phys[c][i];
        norm_squared[c] += (val*val + der*der) * phys_weights[i];
      }
      else norm_squared[c] += val*val * phys_weights[i];
    }
  }

  double elem_norm_squared = 0;
  for (int c=0; c<n_eq; c++)  
    elem_norm_squared += norm_squared[c];

  return elem_norm_squared;
} 

// Evaluate solution and its derivatives in plotting points 'x_phys' 
// in the element (coeffs[][] array not provided).
void Element::get_solution_plot(double x_phys[MAX_PLOT_PTS_NUM], int pts_num,  
                           double val_phys[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
                           double der_phys[MAX_EQN_NUM][MAX_PLOT_PTS_NUM],
                           double *y_prev, 
                           double bc_left_dir_values[MAX_EQN_NUM],
                           double bc_right_dir_values[MAX_EQN_NUM])
{
  double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM];
  this->get_coeffs(y_prev, coeff, bc_left_dir_values, bc_right_dir_values);
  this->get_solution_plot(coeff, pts_num, x_phys, val_phys, der_phys);
} 

// evaluate solution and its derivative 
// at point x_phys (coeffs[][] array provided)
void Element::get_solution_point(double x_phys,
			         double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], 
                                 double val[MAX_EQN_NUM], double der[MAX_EQN_NUM])
{
  double x1 = this->x1;
  double x2 = this->x2;
  double jac = (x2-x1)/2.;
  int dof_size = this->dof_size; 
  int p = this->p;
  // transforming point x_phys to (-1, 1)
  double x_ref = inverse_map(x1, x2, x_phys);
  for(int c=0; c<dof_size; c++) {
    der[c] = val[c] = 0;
    for(int j=0; j<=p; j++) {
      val[c] += coeff[c][j]*lobatto_val_ref(x_ref, j);
      der[c] += coeff[c][j]*lobatto_der_ref(x_ref, j);
    }
    der[c] /= jac;
  }
} 

// evaluate solution and its derivative 
// at the point x_phys (coeffs[][] array not provided)
void Element::get_solution_point(double x_phys,
				 double val_phys[MAX_EQN_NUM], 
                                 double der_phys[MAX_EQN_NUM],
                                 double *y_prev, double *bc_left_dir_values,
                                 double *bc_right_dir_values)
{
  double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM];
  this->get_coeffs(y_prev, coeff, bc_left_dir_values, bc_right_dir_values); 
  this->get_solution_point(x_phys, coeff, val_phys, der_phys);
} 

// Copying elements including all their variables, dof arrays, and 
// refinement trees, 
// FIXME - the recursive version is slow and should be improved.
void Element::copy_element_recursively(Element *e_trg) 
{
  // copy all variables of Element class
  e_trg->init(this->x1, this->x2, this->p, this->id, 
              this->active, this->level, this->dof_size);

  // copy dof arrays for all solution components
  for(int c=0; c < dof_size; c++) {
    for(int i=0; i < MAX_P + 1; i++) e_trg->dof[c][i] = this->dof[c][i];
  }

  // replicate sons if relevant
  if(this->sons[0] != NULL) {          // element was split in space (sons will be replicated)
    e_trg->sons[0] = new Element();
    e_trg->sons[1] = new Element();
    // left son
    this->sons[0]->copy_element_recursively(e_trg->sons[0]);
    // right son
    this->sons[1]->copy_element_recursively(e_trg->sons[1]);
  }
  else {                               // element was split in space (sons are NULL)
    e_trg->sons[0] = NULL;
    e_trg->sons[1] = NULL;
  }
}

// gets physical coordinate of a reference point x_ref \in [-1,1]
double Element::get_x_phys(double x_ref) 
{
  double a = this->x1;
  double b = this->x2;
  return (a+b)/2. + x_ref*(b-a)/2.;
}

Mesh::Mesh() {
  n_eq = 0;
  n_base_elem = 0;
  n_active_elem = 0;
  n_dof = 0;
  base_elems = NULL;
}

// creates equidistant mesh with uniform polynomial degree of elements
Mesh::Mesh(double a, double b, int n_base_elem, int p_init, int n_eq)
{
  // print the banner (only once)
  static int n_calls = 0;
  n_calls++;
  if (n_calls == 1) intro();

  // check maximum number of equations
  if(n_eq > MAX_EQN_NUM) 
  error("Maximum number of equations exceeded (set in common.h)");

  // all Mesh class variables
  this->left_endpoint = a;
  this->right_endpoint = b;
  this->n_base_elem = n_base_elem;
  this->n_eq = n_eq;
  this->n_active_elem = n_base_elem;

  // erase arrays of Dirichlet boundary conditions
  for (int i=0; i<MAX_EQN_NUM; i++) {
    this->bc_left_dir_values[i] = 0;
    this->bc_right_dir_values[i] = 0;
  }

  // allocate element array
  this->base_elems = new Element[this->n_base_elem];     
  if (base_elems == NULL) error("Not enough memory in Mesh::create().");
  if (p_init > MAX_P) 
    error("Max element order exceeded (set in common.h).");
  // element length
  double h = (b - a)/this->n_base_elem;          
  // fill initial element array
  for(int i=0; i<this->n_base_elem; i++) { 
    int id = i;         
    int active = 1;
    int level = 0; 
    this->base_elems[i].init(a + i*h, a + i*h + h, p_init, 
                             id, active, level, n_eq);
  }
}

// creates mesh using a given array of n_base_elem+1 points (pts_array)
// and an array of n_base_elem polynomial degrees (p_array)
Mesh::Mesh(int n_base_elem, double *pts_array, int *p_array, int n_eq)
{
  // print the banner (only once)
  static int n_calls = 0;
  n_calls++;
  if (n_calls == 1) intro();

  // check maximum number of equations
  if(n_eq > MAX_EQN_NUM) 
  error("Maximum number of equations exceeded (set in common.h)");

  // all Mesh class variables
  this->left_endpoint = pts_array[0];
  this->right_endpoint = pts_array[n_base_elem];
  this->n_base_elem = n_base_elem;
  this->n_eq = n_eq;
  this->n_active_elem = n_base_elem;

  // erase arrays of Dirichlet boundary conditions
  for (int i=0; i<MAX_EQN_NUM; i++) {
    this->bc_left_dir_values[i] = 0;
    this->bc_right_dir_values[i] = 0;
  }

  // allocate element array
  this->base_elems = new Element[this->n_base_elem];     
  if (base_elems == NULL) error("Not enough memory in Mesh::create().");

  // initialize element array
  for(int i=0; i<this->n_base_elem; i++) {
    if (p_array[i] > MAX_P) {
      error("Max element order exceeded (set in common.h).");
    }
    int id = i;
    int active = 1;
    int level = 0; 
    this->base_elems[i].init(pts_array[i], pts_array[i+1], 
                             p_array[i], id, active, level, n_eq);
  }
}

// CAUTION - this is expensive (traverses the entire tree 
// from the beginning until the element is found)
void Mesh::refine_single_elem(int id, int3 cand)
{
    Iterator I(this);
    Element *e;
    while ((e = I.next_active_element()) != NULL) {
        printf("%d %d\n", e->id, id);
        if (e->id == id) {
            e->refine(cand);
            if (cand[0] == 1) this->n_active_elem++; // hp-refinement
            return;
        }
    }
    error("refine_single_elem: Element not found.");
}

// performs mesh refinement using a list of elements to be 
// refined and a list of corresponding polynomial degree 
// pairs for the sons
void Mesh::refine_elems(int elem_num, int *id_array, int3 *cand_array)
{
    Iterator *I = new Iterator(this);
    Element *e;
    int count = 0;
    while ((e = I->next_active_element()) != NULL) {
        if (e->id == id_array[count]) {
            if (count >= elem_num)
                error("refine_multi_elems: not enough elems specified");
            e->refine(cand_array[count]);
            if (cand_array[count][0] == 1) this->n_active_elem++;
            count++;
        }
    }
}

// splits the indicated elements and 
// increases poly degree in sons by one
void Mesh::reference_refinement(int start_elem_id, int elem_num)
{
    Iterator *I = new Iterator(this);
    Element *e;
    int count = 0;
    while ((e = I->next_active_element()) != NULL) {
        if (e->id >= start_elem_id && e->id < start_elem_id + elem_num) {
	    if (count >= elem_num) return;
            int3 cand = {1, e->p + 1, e->p + 1};
            e->refine(cand);
            
            if (cand[0] == 1) this->n_active_elem++; // if hp-refinement
            count++;
        }
    }
    this->assign_dofs();
}

void Mesh::set_bc_left_dirichlet(int eq_n, double val)
{
  this->bc_left_dir_values[eq_n] = val;
  // deactivate the corresponding dof for the left-most
  // element and all his descendants adjacent to the 
  // left boundary
  Element *e = this->base_elems + 0;
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
  Element *e = this->base_elems + this->n_base_elem - 1;
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
      for(int j=2; j <= e->p; j++) {
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
        printf("\nElement (%g, %g), id = %d, p = %d\n ", 
               e->x1, e->x2, e->id, e->p); 
        for(int j = 0; j<e->p + 1; j++) {
          printf("dof[%d][%d] = %d\n ", c, j, e->dof[c][j]);
        }
      }
    }
    printf("\n"); 
  }

  delete I;
  return this->n_dof;
}

int Mesh::assign_elem_ids()
{
    Iterator *I = new Iterator(this);
    int count_id = 0;
    Element *e;
    I->reset();
    while ((e = I->next_active_element()) != NULL) {
        e->id = count_id;
        count_id++;
    }
    delete I;
}

Element* Mesh::first_active_element()
{
  Element *e = base_elems;
  while(!e->is_active()) {
    e = e->sons[0];
  }
  return e;
}

Element* Mesh::last_active_element()
{
  Element *e = base_elems + n_base_elem - 1;
  while(!e->is_active()) {
    e = e->sons[1];
  }
  return e;
}

// defining macro to check whether hp candidates are admissible
#define add_hp_candidate_if_ok(new_p_left, new_p_right) \
    candidate_ok = 1; \
    if (new_p_left >= MAX_P) candidate_ok = 0; \
    if (new_p_right >= MAX_P) candidate_ok = 0; \
    if (ref_solution_p_refined) { \
        if (new_p_left >= p_ref && new_p_right >= p_ref) \
            candidate_ok = 0; \
    } else { \
        if (new_p_left == p_ref_left && new_p_right == p_ref_right) \
            candidate_ok = 0; \
        if (new_p_left > p_ref_left || new_p_right > p_ref_right) \
            candidate_ok = 0; \
    } \
    if (candidate_ok) { \
        cand_list[counter][0] = 1; \
        cand_list[counter][1] = new_p_left; \
        cand_list[counter][2] = new_p_right; \
        counter++; \
    }

// defining macro to check whether p candidates are admissible
#define add_p_candidate_if_ok(new_p) \
    candidate_ok = 1; \
    if (new_p >= MAX_P) candidate_ok = 0; \
    if (ref_solution_p_refined) \
        if (new_p >= p_ref) \
            candidate_ok = 0; \
    if (candidate_ok) { \
          cand_list[counter][0] = 0; \
          cand_list[counter][1] = new_p; \
          cand_list[counter][2] = -1; \
          counter++; \
        } \

/*
   p_ref_left and p_ref_right are the polynomial orders of the reference
   solution (left and right refinement). If p_ref_right is -1, then the
   reference solution is only "p" refined.
   adapt_type = 0 ... hp-adaptivity (full list of candidates)
   adapt_type = 1 ... h-adaptivity (h-candidates only)
   */
int Element::create_cand_list(int adapt_type, int p_ref_left, 
                              int p_ref_right, int3 *cand_list) 
{
    int ref_solution_p_refined = (p_ref_right == -1);
    int p_ref = p_ref_left; // use only if ref_solution_p_refined == 1
    int counter = 0;
    int candidate_ok;

    if (adapt_type == 0) {
      // p-adaptivity candidates: 
      add_p_candidate_if_ok(this->p + 1);
      add_p_candidate_if_ok(this->p + 2);

      // hp-adaptivity candidates: 
      // we first divide the poly degree by two
      int base_p = this->p / 2; 
      if (base_p < 1) base_p = 1;

      add_hp_candidate_if_ok(base_p, base_p); // p -> (p/2, p/2) 
      add_hp_candidate_if_ok(base_p + 1, base_p); // p -> (p/2+1, p/2) 
      add_hp_candidate_if_ok(base_p, base_p + 1); // p -> (p/2, p/2+1) 
      add_hp_candidate_if_ok(base_p + 1, base_p + 1); // p -> (p/2+1, p/2+1) 
      add_hp_candidate_if_ok(base_p + 2, base_p); // p -> (p/2+2, p/2) 
      add_hp_candidate_if_ok(base_p, base_p + 2); // p -> (p/2, p/2+2) 
      add_hp_candidate_if_ok(base_p + 1, base_p + 2); // p -> (p/2+1, p/2+2) 
      add_hp_candidate_if_ok(base_p + 2, base_p + 1); // p -> (p/2+2, p/2+1) 
      add_hp_candidate_if_ok(base_p + 2, base_p + 2); // p -> (p/2+2, p/2+2) 
    }
    if (adapt_type == 1) {
      add_hp_candidate_if_ok(this->p, this->p); // p -> (p, p) 
    }
    if (adapt_type == 2) {
      add_p_candidate_if_ok(this->p + 1);
      add_p_candidate_if_ok(this->p + 2);
    }

    return counter;
}

void Element::print_cand_list(int num_cand, int3 *cand_list) 
{
  printf("Element (%g, %g): refinement candidates:\n", this->x1, this->x2);
  for (int i=0; i < num_cand; i++) { 
    printf("%d %d %d\n", cand_list[i][0], cand_list[i][1], cand_list[i][2]);
  }
}

// transformation of k-th shape function defined on Gauss points 
// corresponding to 'order' to physical interval (a,b)
void element_shapefn(double a, double b, 
		     int k, int order, double *val, double *der) {
  double2 *ref_tab = g_quad_1d_std.get_points(order);
  int pts_num = g_quad_1d_std.get_num_points(order);
  for (int i=0 ; i<pts_num; i++) {
    // change function values and derivatives to interval (a, b)
    val[i] = lobatto_val_ref(ref_tab[i][0], k);
    double jac = (b-a)/2.; 
    der[i] = lobatto_der_ref(ref_tab[i][0], k) / jac; 
  }
};

// transformation of k-th shape function at the reference 
// point x_ref to physical interval (a,b).
void element_shapefn_point(double x_ref, double a, double b, 
		           int k, double *val, double *der) {
    // change function values and derivatives to interval (a, b)
  *val = lobatto_val_ref(x_ref, k);
    double jac = (b-a)/2.; 
    *der = lobatto_der_ref(x_ref, k) / jac; 
}

// Replicate mesh including dof arrays in all elements
Mesh *Mesh::replicate()
{
  // copy base mesh element array, use dummy poly degree first 
  // (poly degrees in base mesh may have changed since initialization)
  int p_dummy = -1; 
  Mesh *mesh_new = new Mesh(this->left_endpoint, this->right_endpoint, 
			    this->n_base_elem, p_dummy, this->n_eq);

  // copy all Mesh class variables
  mesh_new->set_n_eq(this->n_eq);
  mesh_new->set_n_base_elem(this->n_base_elem);
  mesh_new->set_n_active_elem(this->n_active_elem);
  mesh_new->set_left_endpoint(this->left_endpoint);
  mesh_new->set_right_endpoint(this->right_endpoint);
  mesh_new->set_n_dof(this->n_dof);

  // copy links to Dirichlet boundary condition arrays
  for(int c = 0; c<this->n_eq; c++) {
    mesh_new->bc_left_dir_values[c] = this->bc_left_dir_values[c]; 
    mesh_new->bc_right_dir_values[c] = this->bc_right_dir_values[c]; 
  }

  // replicate all base mesh elements including all their 
  // variables, dof arrays, and tree-structure
  Element *base_elems_new = mesh_new->get_base_elems();
  for(int i=0; i<this->n_base_elem; i++) {
    Element *e_src = this->base_elems + i;
    Element *e_trg = base_elems_new + i;
    e_src->copy_element_recursively(e_trg);
  }

  return mesh_new;
}

void Mesh::plot(const char* filename) 
{
    FILE *f = fopen(filename, "wb");
    if(f == NULL) error("problem opening file in Mesh::plot().");
    Iterator I(this);
    Element *e;
    while ((e = I.next_active_element()) != NULL) {
      fprintf(f, "%g %d\n", e->x1, 0);
      fprintf(f, "%g %d\n", e->x1, e->p);
      fprintf(f, "%g %d\n", e->x2, e->p);
      fprintf(f, "%g %d\n\n", e->x2, 0);
    }
    fclose(f);
    printf("Mesh written to %s.\n", filename);
}

// Plots the error between the reference and coarse mesh solutions
// if the reference refinement was p-refinement
void Mesh::plot_element_error_p(int norm, FILE *f, Element *e, Element *e_ref, 
                                double *y_prev, double *y_prev_ref,
                                int subdivision)
{
  int n_eq = this->get_n_eq();
  int pts_num = subdivision + 1;
  if (pts_num > MAX_PLOT_PTS_NUM) {
    printf("Try to increase MAX_PLOT_PTS_NUM in common.h\n");
    error("MAX_PLOT_PTS_NUM exceeded in plot_element_error_p().");
  }
  double x1 = e->x1;
  double x2 = e->x2;

  // calculate point array
  double x_phys[MAX_PLOT_PTS_NUM];  
  double h = (x2 - x1)/subdivision; 
  for (int i=0; i < pts_num; i++) {
    x_phys[i] = x1 + i*h;
  }

  // get coarse mesh solution values and derivatives
  double phys_u[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
         phys_dudx[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];
  e->get_solution_plot(x_phys, pts_num, phys_u, phys_dudx, y_prev, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  // get fine mesh solution values and derivatives
  double phys_u_ref[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
         phys_dudx_ref[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];
  e->get_solution_plot(x_phys, pts_num, phys_u_ref, phys_dudx_ref, y_prev_ref, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  for (int i=0; i < pts_num; i++) {
    double diff_squared_pt = 0;
    for (int c=0; c < n_eq; c++) {
      double val = phys_u_ref[c][i] - phys_u[c][i]; 
      double der = phys_dudx_ref[c][i] - phys_dudx[c][i]; 
      diff_squared_pt += val*val;
      if (norm == 1) diff_squared_pt += der*der;
    }
    fprintf(f, "%g %g\n", x_phys[i], sqrt(diff_squared_pt));
  }
  fprintf(f, "\n");
}

// Plots the error between the reference and coarse mesh solutions
// if the reference refinement was hp-refinement
void Mesh::plot_element_error_hp(int norm, FILE *f, Element *e, 
                                 Element *e_ref_left, 
                                 Element *e_ref_right, 
                                 double *y_prev, double *y_prev_ref,
                                 int subdivision)
{
  // We will be plotting the error separately in the 
  // elements e_ref_left (left half of 'e') and 
  // e_ref_right (right half of 'e').
  subdivision /= 2;
  int pts_num = subdivision + 1;
  if (pts_num > MAX_PLOT_PTS_NUM) {
    printf("Try to increase MAX_PLOT_PTS_NUM in common.h\n");
    error("MAX_PLOT_PTS_NUM exceeded in plot_element_error_hp().");
  }

  // First: element e_ref_left
  double x1 = e_ref_left->x1;
  double x2 = e_ref_left->x2;

  // calculate array of x-coordinates
  double x_phys_left[MAX_PLOT_PTS_NUM];  
  double h = (x2 - x1)/subdivision; 
  for (int i=0; i < pts_num; i++) {
    x_phys_left[i] = x1 + i*h;
  }

  // get coarse mesh solution values and derivatives
  double phys_u_left[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
         phys_dudx_left[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];
  e->get_solution_plot(x_phys_left, pts_num, phys_u_left, phys_dudx_left, y_prev, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  // get fine mesh solution values and derivatives
  double phys_u_ref_left[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
         phys_dudx_ref_left[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];
  e_ref_left->get_solution_plot(x_phys_left, pts_num, phys_u_ref_left, 
                  phys_dudx_ref_left, y_prev_ref, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  for (int i=0; i < pts_num; i++) {
    double diff_squared_pt = 0;
    for (int c=0; c < n_eq; c++) {
      double val = phys_u_ref_left[c][i] - phys_u_left[c][i]; 
      double der = phys_dudx_ref_left[c][i] - phys_dudx_left[c][i]; 
      diff_squared_pt += val*val;
      if (norm == 1) diff_squared_pt += der*der;
    }
    fprintf(f, "%g %g\n", x_phys_left[i], sqrt(diff_squared_pt));
    fprintf(f, "\n");
  }

  // Second: element e_ref_right
  x1 = e_ref_right->x1;
  x2 = e_ref_right->x2;

  // calculate array of x-coordinates
  double x_phys_right[MAX_PLOT_PTS_NUM];  
  h = (x2 - x1)/subdivision; 
  for (int i=0; i < pts_num; i++) {
    x_phys_right[i] = x1 + i*h;
  }

  // get coarse mesh solution values and derivatives
  double phys_u_right[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
         phys_dudx_right[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];
  e->get_solution_plot(x_phys_right, pts_num, phys_u_right, phys_dudx_right, y_prev, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  // get fine mesh solution values and derivatives
  double phys_u_ref_right[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
         phys_dudx_ref_right[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];
  e_ref_right->get_solution_plot(x_phys_right, pts_num, phys_u_ref_right, 
    phys_dudx_ref_right, y_prev_ref, this->bc_left_dir_values, 
    this->bc_right_dir_values); 

  for (int i=0; i < pts_num; i++) {
    double diff_squared_pt = 0;
    for (int c=0; c < n_eq; c++) {
      double val = phys_u_ref_right[c][i] - phys_u_right[c][i]; 
      double der = phys_dudx_ref_right[c][i] - phys_dudx_right[c][i]; 
      diff_squared_pt += val*val;
      if (norm == 1) diff_squared_pt += der*der;
    }
    fprintf(f, "%g %g\n", x_phys_right[i], sqrt(diff_squared_pt));
    fprintf(f, "\n");
  }
}

// Plots the error wrt. the exact solution (if available)
void Mesh::plot_element_error_exact(int norm, FILE *f, Element *e, 
                                    double *y_prev, exact_sol_type exact_sol, int subdivision)
{
  int pts_num = subdivision + 1;
  double x1 = e->x1;
  double x2 = e->x2;

  // calculate array of x-coordinates
  double x_phys[MAX_PLOT_PTS_NUM];  
  double h = (x2 - x1)/subdivision; 
  for (int i=0; i < pts_num; i++) {
    x_phys[i] = x1 + i*h;
  }

  // get coarse mesh solution values and derivatives
  double phys_u[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
         phys_dudx[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];
  e->get_solution_plot(x_phys, pts_num, phys_u, phys_dudx, y_prev, 
                  this->bc_left_dir_values, this->bc_right_dir_values); 

  for (int i=0; i < pts_num; i++) {
    double exact_sol_point[MAX_EQN_NUM];
    double exact_der_point[MAX_EQN_NUM];
    exact_sol(x_phys[i], exact_sol_point, exact_der_point);
    double diff_squared_pt = 0; 
    for (int c=0; c < n_eq; c++) {
      double val = exact_sol_point[c] - phys_u[c][i];
      double der = exact_der_point[c] - phys_dudx[c][i]; 
      diff_squared_pt += val*val;
      if (norm == 1) diff_squared_pt += der*der;
    }
    fprintf(f, "%g %g\n", x_phys[i], sqrt(diff_squared_pt));
  }
  fprintf(f, "\n");
}

// Plots the error between the reference and coarse mesh solutions
void Mesh::plot_error_est(int norm, const char *filename, Mesh* mesh_ref, 
		          double* y_prev, double* y_prev_ref, int subdivision)
{
  char final_filename[MAX_STRING_LENGTH];
  sprintf(final_filename, "%s", filename);
  FILE *f = fopen(final_filename, "wb");
  if(f == NULL) error("problem opening file in plot_error_est().");

  // simultaneous traversal of 'this' and 'mesh_ref'
  Element *e;
  Iterator *I = new Iterator(this);
  Iterator *I_ref = new Iterator(mesh_ref);
  while ((e = I->next_active_element()) != NULL) {
    Element *e_ref = I_ref->next_active_element();
    if (e->level == e_ref->level) { // element 'e' was not refined in space
                                    // for reference solution
      if (e_ref->p >= MAX_P) {
        printf("Try to increase MAX_P in common.h.\n");
        error("Max poly degree exceeded in plot_error_est().");
      }
      plot_element_error_p(norm, f, e, e_ref, y_prev, y_prev_ref, subdivision);
    }
    else { // element 'e' was refined in space for reference solution
      Element* e_ref_left = e_ref;
      Element* e_ref_right = I_ref->next_active_element();
      if (e_ref_left->p >= MAX_P || e_ref_right->p >= MAX_P) {
        printf("Try to increase MAX_P in common.h.\n");
        error("Max poly degree exceeded in plot_error_est().");
      }
      plot_element_error_hp(norm, f, e, e_ref_left, e_ref_right, 
                            y_prev, y_prev_ref, subdivision);
    }
  }

  fclose(f);
  printf("Error function written to %s.\n", final_filename);
}

// Plots the error wrt. the exact solution (if available)
void Mesh::plot_error_exact(int norm, const char *filename,  
		            double* y_prev, exact_sol_type exact_sol, int subdivision)
{
  char final_filename[MAX_STRING_LENGTH];
  sprintf(final_filename, "%s", filename);
  FILE *f = fopen(final_filename, "wb");
  if(f == NULL) error("problem opening file in plot_error_exact().");

  // traversal of 'this'
  Element *e;
  Iterator *I = new Iterator(this);
  while ((e = I->next_active_element()) != NULL) {
    if (e->p >= MAX_P) {
      printf("Try to increase MAX_P in common.h.\n");
      error("Max poly degree exceeded in plot_error_est().");
    }
    plot_element_error_exact(norm, f, e, y_prev, exact_sol, subdivision);
  }

  fclose(f);
  printf("Exact solution error written to %s.\n", final_filename);
}

// Refine coarse mesh elements whose id_array >= 0, and 
// adjust the reference mesh accordingly.  
// Returns updated coarse and reference meshes, with the last 
// coarse and reference mesh solutions on them, respectively. 
// The coefficient vectors and numbers of degrees of freedom 
// on both meshes are also updated. 
void adapt(int norm, int adapt_type, double threshold, 
           double *err_array,
           Mesh* &mesh, Mesh* &mesh_ref, 
           double * &y_prev, double* &y_prev_ref, 
           int &n_dof, int &n_dof_ref) 
{
  int n_elem = mesh->get_n_active_elem();

  // debug
  if (PRINT_ELEM_ERRORS) {
    printf("  Element errors (est) in adapt():\n");
    for(int i=0; i < n_elem; i++) {
      printf("  Elem [%d]: %g\n", i, err_array[i]);
    }
  }

  // Find element with largest error
  double max_elem_error = 0;
  for(int i=0; i < n_elem; i++) {
    double elem_error = err_array[i];
    if (elem_error > max_elem_error) {
      max_elem_error = elem_error;
    }
  }

  // Create auxiliary array of element indices
  int id_array[MAX_ELEM_NUM];
  for(int i=0; i < n_elem; i++) {
   if(err_array[i] < threshold*max_elem_error) id_array[i] = -1; 
   else id_array[i] = i;
  }

  // Print elements to be refined
  /*
  printf("Elements to be refined:\n");
  for (int i=0; i<this->get_n_active_elem(); i++) {
    if (id_array[i] >= 0) printf("Elem[%d], error = %g\n", id_array[i], 
                                 err_array[i]);
  }
  */

  int adapt_list[MAX_ELEM_NUM];
  int num_to_adapt = 0;

  // Create list of elements to be refined, in increasing order
  for (int i=0; i < n_elem; i++) {
    if (id_array[i] >= 0) {
      adapt_list[num_to_adapt] = id_array[i];
      num_to_adapt++;
    }
  }
 
  // Debug: Printing list of elements to be refined
  //printf("Elements to be refined: ");
  //for (int i=0; i<num_to_adapt; i++) printf("%d ", adapt_list[i]);
  //printf("\n");

  // Replicate the coarse and fine meshes. The original meshes become
  // backup and refinements will only be done in the new ones.
  Mesh *mesh_new = mesh->replicate();
  Mesh *mesh_ref_new = mesh_ref->replicate();

  // Simultaneous traversal of all meshes.
  // For each element in 'mesh_new', create a list of refinement 
  // candidates and select the one that best resembles the reference 
  // solution on 'mesh_ref_new'. While elements in 'mesh_new' are refined, 
  // corresponding refinements are also done in 'mesh_ref_new'.
  Iterator *I = new Iterator(mesh);
  Iterator *I_new = new Iterator(mesh_new);
  Iterator *I_ref = new Iterator(mesh_ref);
  Iterator *I_ref_new = new Iterator(mesh_ref_new);
  Element *e = I->next_active_element();
  Element *e_new = I_new->next_active_element();
  Element *e_ref = I_ref->next_active_element();
  Element *e_ref_new = I_ref_new->next_active_element();
  int counter = 0;
  while (counter != num_to_adapt) {
    if (e->id == adapt_list[counter]) {
      counter++;
      int choice = 0;
      int3 cand_list[MAX_CAND_NUM];   // Every refinement candidates consists of three
                                      // numbers: 1/0 whether it is a p-candidate or not,
                                      // and then either one or two polynomial degrees
      Element* e_ref_left;
      Element* e_ref_right;
      Element* e_ref_new_left;
      Element* e_ref_new_right;
      if (e->level == e_ref->level) { // element 'e' was not refined in space
                                      // for reference solution
        e_ref_left = e_ref;
        e_ref_new_left = e_ref_new;
        e_ref_right = NULL;
        e_ref_new_right = NULL;
        int num_cand = e->create_cand_list(adapt_type, e_ref->p, -1, cand_list);
        // debug:
        //e->print_cand_list(num_cand, cand_list);
        // reference element was p-refined
        choice = select_hp_refinement(e, e_ref, NULL, num_cand, cand_list, 
                                      0, norm, y_prev_ref, 
                                      mesh->bc_left_dir_values,
	  		              mesh->bc_right_dir_values);
      }
      else { // element 'e' was refined in space for reference solution
        e_ref_left = e_ref;
        e_ref_new_left = e_ref_new;
        e_ref_right = I_ref->next_active_element();
        e_ref_new_right = I_ref_new->next_active_element();
        int num_cand = e->create_cand_list(adapt_type, e_ref_left->p, 
                                           e_ref_right->p, cand_list);
        // reference element was hp-refined
        choice = select_hp_refinement(e, e_ref_left, e_ref_right, 
                                      num_cand, cand_list, 1, norm, y_prev_ref, 
                                      mesh->bc_left_dir_values,
			              mesh->bc_right_dir_values);
      }

      // Next we perform the refinement defined by cand_list[choice]
      // e_new_last... element in mesh_new that will be refined,
      // e_ref_left... corresponding element in fine mesh (if reference 
      //               refinement was p-refinement). In this case 
      //               e_ref_right == NULL
      // e_ref_left, e_ref_right... corresponding pair of elements 
      //               in the fine mesh if reference refinement was hp-refinement
      Element *e_new_last = e_new;
      // moving pointers 'e' and 'e_new' to the next position in coarse meshes
      // and 'e_ref' and 'e_ref_new' to the next position in fine meshes 
      e = I->next_active_element();
      e_new = I_new->next_active_element();
      e_ref = I_ref->next_active_element();
      e_ref_new = I_ref_new->next_active_element();
      // perform the refinement of element e_new_last
      e_new_last->refine(cand_list[choice]);
      printf("  Refined element (%g, %g), cand = (%d %d %d)\n", 
             e_new_last->x1, e_new_last->x2, cand_list[choice][0], 
             cand_list[choice][1], cand_list[choice][2]);
      if(cand_list[choice][0] == 1) mesh_new->n_active_elem++; 
      // perform corresponding refinement(s) in the new fine mesh
      if (e_new_last->level == e_ref_left->level) { // ref. refinement of 'e_last' was 
                                                    // p-refinement so also future ref. 
                                                    // refinements will be p-refinements
        if (cand_list[choice][0] == 0) { // e_last is being p-refined, thus also
                                         // e_ref_new_left needs to be p-refined
          int new_p = cand_list[choice][1];
	  e_ref_new_left->refine(0, new_p + 1, -1);
        }
        else { // e_new_last is being split, thus e_ref_new_left needs to be 
               // split as well
          int new_p_left = cand_list[choice][1];
          int new_p_right = cand_list[choice][2];
	  e_ref_new_left->refine(1, new_p_left + 1, new_p_right + 1);
          mesh_ref_new->n_active_elem++;
        }
      }
      else { // ref. refinement was hp-refinement, so also future
             // ref. refinements will be hp-refinements
        if (cand_list[choice][0] == 0) { // e_new_last is being p-refined, thus also 
                                         // e_ref_new_left and e_ref_new_right  
                                         // will just be p-refined
          int new_p = cand_list[choice][1];
	  e_ref_new_left->refine(0, new_p + 1, -1);
	  e_ref_new_right->refine(0, new_p + 1, -1);
        }
        else { // e_new_last is being hp-refined, so we need to 
               // split both e_ref_new_left and e_ref_new_right
          int new_p_left = cand_list[choice][1];
          int new_p_right = cand_list[choice][2];
	  e_ref_new_left->refine(1, new_p_left + 1, new_p_left + 1);
          mesh_ref_new->n_active_elem++;
	  e_ref_new_right->refine(1, new_p_right + 1, new_p_right + 1);
          mesh_ref_new->n_active_elem++;
       }
      }
    }
    else {
      e = I->next_active_element();
      e_new = I_new->next_active_element();
      e_ref = I_ref->next_active_element();
      e_ref_new = I_ref_new->next_active_element();
      if (e->level != e_ref->level) { 
        e_ref = I_ref->next_active_element();
        e_ref_new = I_ref_new->next_active_element();
      }    
    }
  }
  // enumerate dofs in both new meshes
  int n_dof_new = mesh_new->assign_dofs();
  int n_dof_ref_new = mesh_ref_new->assign_dofs();
  printf("Coarse mesh refined (%d elem, %d DOF)\n", 
         mesh_new->get_n_active_elem(), n_dof_new);
  printf("Fine mesh updated (%d elem, %d DOF)\n", 
  	 mesh_ref_new->get_n_active_elem(), n_dof_ref_new);

  // Transfer last coarse and fine mesh solutions to the new coarse and 
  // fine meshes, respectively. 
  double *y_prev_new = new double[n_dof_new];
  double *y_prev_ref_new = new double[n_dof_ref_new];
  transfer_solution(mesh, mesh_new, y_prev, y_prev_new);
  printf("Last coarse mesh solution copied to new coarse mesh.\n");
  transfer_solution(mesh_ref, mesh_ref_new, y_prev_ref, y_prev_ref_new);
  printf("Last fine mesh solution copied to new fine mesh.\n");

  // Delete old meshes and copy the new ones in place of them
  delete mesh;
  delete mesh_ref;
  mesh = mesh_new;
  mesh_ref = mesh_ref_new;

  // Do the same with the coefficients vectors.
  delete [] y_prev;
  delete [] y_prev_ref;
  y_prev = y_prev_new;
  y_prev_ref = y_prev_ref_new;

  // Last adjust the number of dofs in each mesh
  n_dof = n_dof_new;
  n_dof_ref = n_dof_ref_new;
}

void adapt_plotting(Mesh *mesh, Mesh *mesh_ref, double *y_prev, 
              double *y_prev_ref, int norm, int exact_sol_provided, 
              exact_sol_type exact_sol) 
{
  // Plot the coarse mesh solution
  Linearizer l(mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, y_prev);

  // Plot the fine mesh solution
  Linearizer lxx(mesh_ref);
  const char *out_filename2 = "solution_ref.gp";
  lxx.plot_solution(out_filename2, y_prev_ref);

  // Plot the coarse and fine meshes
  const char *mesh_filename = "mesh.gp";
  mesh->plot(mesh_filename);
  const char *mesh_ref_filename = "mesh_ref.gp";
  mesh_ref->plot(mesh_ref_filename);

  // Plot the error estimate (difference between 
  // coarse and fine mesh solutions)
  const char *err_est_filename = "error_est.gp";
  mesh->plot_error_est(norm, err_est_filename, mesh_ref, y_prev, y_prev_ref);

  // Plot error wrt. exact solution (if available)
  if (exact_sol_provided) {   
    const char *err_exact_filename = "error_exact.gp";
    mesh->plot_error_exact(norm, err_exact_filename, y_prev, exact_sol);
  }
}
