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
  sons[0] = sons[1] = NULL; 
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
  sons[0] = sons[1] = NULL; 
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
  this->dof_size = n_eq;
  // c is solution component
  for(int c=0; c<n_eq; c++) {
    this->dof[c] = new int[MAX_POLYORDER + 1];
    // important for th etreatment of boundary conditions
    for(int i=0; i<MAX_POLYORDER + 1; i++) this->dof[c][i] = 0;
  }
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
  // number of base elements
  this->n_base_elem = n_base_elem;
  // number of active elements
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
    // allocate element dof arrays for all solution components 
    // and length MAX_POLYORDER
    this->base_elems[i].dof_alloc(n_eq);
    // define element end points
    this->base_elems[i].x1 = a + i*h;
    this->base_elems[i].x2 = this->base_elems[i].x1 + h;
  }
  this->assign_elem_ids();
}

// creates mesh using a given array of n_base_elem+1 points (pts_array)
// and an array of n_base_elem polynomial degrees (p_array)
Mesh::Mesh(int n_base_elem, double *pts_array, int *p_array, int n_eq)
{
  // domain end points
  left_endpoint = pts_array[0];
  right_endpoint = pts_array[n_base_elem];
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
  // number of base elements
  this->n_base_elem = n_base_elem;
  // number of active elements
  this->n_active_elem = n_base_elem;
  // allocate element array
  this->base_elems = new Element[this->n_base_elem];     
  if (base_elems == NULL) error("Not enough memory in Mesh::create().");
  // fill initial element array
  for(int i=0; i<this->n_base_elem; i++) {
    if (p_array[i] > MAX_POLYORDER) 
      error("Max element order exceeded (set in common.h).");
    // polynomial degree
    this->base_elems[i].p = p_array[i];
    // allocate element dof arrays for all solution components 
    // and length MAX_POLYORDER
    this->base_elems[i].dof_alloc(n_eq);
    // define element end points
    this->base_elems[i].x1 = pts_array[i];
    this->base_elems[i].x2 = pts_array[i+1];
  }
  this->assign_elem_ids();
}

// caution - this is expensive (traverses the entire tree 
// from the beginning until the element is found)
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

// performs mesh refinement using a list of elements to be 
// refined and a list of corresponding polynomial degree 
// pairs for the sons
void Mesh::refine_elems(int elem_num, int *id_array, int2 *p_pair_array)
{
    Iterator *I = new Iterator(this);
    Element *e;
    int count = 0;
    while ((e = I->next_active_element()) != NULL) {
        if (e->id == id_array[count]) {
            if (count >= elem_num)
                error("refine_multi_elems: not enough elems specified");
            e->refine(p_pair_array[count][0], p_pair_array[count][1], this->n_eq);
            this->n_active_elem++;
            count++;
        }
    }
}

// splits the indicated elements and 
// increases poly degree in sons by one
void Mesh::refine_elems(int start_elem_id, int elem_num)
{
    Iterator *I = new Iterator(this);
    Element *e;
    int count = 0;
    while ((e = I->next_active_element()) != NULL) {
        if (e->id >= start_elem_id && e->id < start_elem_id + elem_num) {
	    if (count >= elem_num) return;
            e->refine(e->p + 1, e->p + 1, this->n_eq);
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
        e->id = count_id++;
    }
    delete I;
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

// copying elements including their refinement trees 
// inactive Dirichlet DOF are replicated
// FIXME - the recursive version is slow, improve it!
void copy_element_sons_recursively(Element *e_src, Element *e_trg, int n_eq) {
  // if element has been refined
  if(e_src->sons[0] != NULL) {
    int p_left = e_src->sons[0]->p;
    int p_right = e_src->sons[1]->p;
    e_trg->refine(p_left, p_right, n_eq);
    // left son
    copy_element_sons_recursively(e_src->sons[0], e_trg->sons[0], n_eq);
    // right son
    copy_element_sons_recursively(e_src->sons[1], e_trg->sons[1], n_eq);
  }
}

// Replicate mesh
Mesh *Mesh::replicate()
{
  // copy base mesh element array, use dummy poly degree (p_init
  // cannot be used since p-refinements on base mesh may have taken 
  // place)
  int p_dummy = -1; 
  Mesh *mesh_ref = new Mesh(this->left_endpoint, this->right_endpoint, 
			    this->n_base_elem, p_dummy, this->n_eq);

  // copy element degrees on base mesh
  for(int i=0; i<this->n_base_elem; i++) {
    mesh_ref->get_base_elems()[i].p = this->get_base_elems()[i].p;
  }

  // copy number of base elements
  mesh_ref->set_n_base_elem(this->get_n_base_elem());
    
  // copy arrays of Dirichlet boundary conditions
  for(int c = 0; c<this->n_eq; c++) {
    mesh_ref->bc_left_dir_values[c] = this->bc_left_dir_values[c]; 
    mesh_ref->bc_right_dir_values[c] = this->bc_right_dir_values[c]; 
  }

  // copy inactive Dirichlet DOF on first and last element of base mesh
  for(int c=0; c<this->n_eq; c++) {
    if(this->base_elems[0].dof[c][0] == -1) {
      Element *e_left = mesh_ref->get_base_elems();
      e_left->dof[c][0] = -1;
    }
    if(this->base_elems[this->n_base_elem-1].dof[c][1] == -1) {
      Element *e_right = mesh_ref->get_base_elems()+this->n_base_elem-1;
      e_right->dof[c][1] = -1;
    }
  }

  // replicate tree structure on all base mesh elements,
  // this includes replication of inactive Dirichlet DOF
  for(int i=0; i<this->n_base_elem; i++) {
    Element *e_src = this->base_elems + i;
    Element *e_trg = mesh_ref->get_base_elems() + i;
    copy_element_sons_recursively(e_src, e_trg, this->n_eq);
  }

  // enumerate elements in reference mesh
  mesh_ref->assign_elem_ids();

  return mesh_ref;
}

