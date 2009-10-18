// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _MESH_H_
#define _MESH_H_

#include "common.h"
#include "lobatto.h"
#include "quad_std.h"

class Element {
public:
  Element();
  Element(double x_left, double x_right, int deg, int n_eq);
  void dof_alloc(int n_eq);
  void refine(int p_left, int p_right, int n_eq);
  unsigned is_active();
  unsigned active;   // flag used by assembling algorithm
  double x1, x2;     // endpoints
  int p;             // poly degree
  int **dof;         // connectivity array of length p+1 for every solution component
  int id;
  unsigned level;    // refinement level (zero for initial mesh elements) 
  Element *sons[2];  // for refinement
};

class Mesh {
    public:
        Mesh(double a, double b, int n_elem, int p_init, int n_eq);
        int assign_dofs();
        Element *get_base_elems() {
	  return this->base_elems;
        }
        int get_n_base_elems() {
	  return this->n_base_elem;
        }
        int get_n_active_elems();
        int get_n_dof() {
	  return this->n_dof;
        }
        int get_n_eq() {
	  return this->n_eq;
        }
        void calculate_elem_coeffs(Element *e, double *y_prev, double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM]);
        void element_solution(Element *e, double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], int pts_num, 
		      double pts_array[MAX_PTS_NUM], double val[MAX_EQN_NUM][MAX_PTS_NUM], 
                      double der[MAX_EQN_NUM][MAX_PTS_NUM]);
        void element_solution_point(double x_ref, Element *e, 
			    double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], double *val, double *der);
        void element_shapefn(double a, double b, 
			     int k, int order, double *val, double *der);
        void element_shapefn_point(double x_ref, double a, double b, 
				   int k, double *val, double *der);
        void set_bc_left_dirichlet(int eq_n, double val);
        void set_bc_right_dirichlet(int eq_n, double val);

        void refine_single_elem(int id, int p_left, int p_right);
        void refine_multi_elems(int n, int *id_array, int2 *p_id_array);

        double *bc_left_dir_values;  // values for the Dirichlet condition left
        double *bc_right_dir_values; // values for the Dirichlet condition right
        double left_endpoint, right_endpoint;

    private:
        int n_eq;
        int n_base_elem;
        int n_active_elem;
        int n_dof;
        Element *base_elems;

        int assign_elem_ids();
};

class Linearizer {
    public:
        Linearizer(Mesh *mesh) {
            this->mesh = mesh;
        }

        // evaluate approximate solution at element 'm' at reference
        // point 'x_ref'. Here 'y' is the global vector of coefficients
        void eval_approx(Element *e, double x_ref, double *y, double *x_phys,
			 double *val);

        void plot_solution(const char *out_filename, double *y_prev, int
                // FIXME: code needs to be fixed to allow
                // plotting_elem_subdivision to be 100 and more
                plotting_elem_subdivision=50);

        void get_xy(double *y_prev, int comp, int plotting_elem_subdivision,
                double **x, double **y, int *n);

    private:
        Mesh *mesh;
};

#endif
