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
    Element(double x_left, double x_right, int lev, int deg, int n_eq);
    void free_element() {
        if (this->sons[0] != NULL) delete this->sons[0];
        if (this->sons[1] != NULL) delete this->sons[1];
        if (this->dof != NULL) {
            for(int c=0; c < this->dof_size; c++)
                delete[] this->dof[c];
            if (this->dof != NULL)
                free(this->dof);
        }
        /*
        this->dof = NULL;
        this->sons[0] = NULL;
        this->sons[1] = NULL;
        */
    }
    ~Element() {
        this->free_element();
    }
    virtual void dof_alloc();
    void init(double x1, double x2, int p_init, int n_eq);
    void copy_sons_recursively(Element *e_trg);
    double get_x_phys(double x_ref); // gets physical coordinate of a reference poin
    void get_coeffs(double *y_prev, 
		    double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM],
                    double bc_left_dir_values[MAX_EQN_NUM],
                    double bc_right_dir_values[MAX_EQN_NUM]);
    void get_solution(double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], 
         int pts_num, double x_phys[MAX_PTS_NUM], 
         double val_phys[MAX_EQN_NUM][MAX_PTS_NUM], 
         double der_phys[MAX_EQN_NUM][MAX_PTS_NUM]);
    void get_solution(double x_phys[MAX_PTS_NUM], int pts_num,  
                      double val_phys[MAX_EQN_NUM][MAX_PTS_NUM], 
                      double der_phys[MAX_EQN_NUM][MAX_PTS_NUM],
                      double *y_prev, 
                      double bc_left_dir_values[MAX_EQN_NUM],
                      double bc_right_dir_values[MAX_EQN_NUM]);
    void get_solution_point(double x_ref,
         double coeff[MAX_EQN_NUM][MAX_COEFFS_NUM], 
	 double val[MAX_EQN_NUM], double der[MAX_EQN_NUM]);
    void get_solution_point(double x_phys,
			    double val[MAX_EQN_NUM], double der[MAX_EQN_NUM], 
                            double *y_prev, double *bc_left_dir_values,
                            double *bc_right_dir_values);
    int create_cand_list(int p_ref_left, int p_ref_right, int3 *cand_list);
    void print_cand_list(int num_cand, int3 *cand_list);
    void refine(int3 cand);
    unsigned is_active();
    unsigned active;   // flag used by assembling algorithm
    double x1, x2;     // endpoints
    int p;             // poly degrees
    int dof_size;      // size of the dof[] array
    int **dof;         // connectivity array of length p+1 
                       // for every solution component
    int id;
    unsigned level;    // refinement level (zero for initial mesh elements) 
    Element *sons[2];  // for refinement
};

class Mesh {
    public:
        Mesh();
        Mesh(double a, double b, int n_elem, int p_init, int n_eq);
        Mesh(int n_base_elem, double *pts_array, int *p_array, int n_eq);
        ~Mesh() {
            if (this->base_elems != NULL) {
                delete[] this->base_elems;
            }
        }
        int assign_dofs();
        Element *get_base_elems() {
            return this->base_elems;
        }

        int get_n_base_elem() {
            return this->n_base_elem;
        }
        void set_n_base_elem(int n_base_elem) {
            this->n_base_elem = n_base_elem;
        }
        int get_n_active_elem() {
            return this->n_active_elem;
        }
        void set_n_active_elem(int n_active_elem) {
            this->n_active_elem = n_active_elem;
        }
        int get_n_dof() {
            return this->n_dof;
        }
        int get_n_eq() {
            return this->n_eq;
        }
        void set_n_eq(int n_eq) {
            this->n_eq = n_eq;
        }
        double get_left_endpoint() {
            return left_endpoint; 
        }
        double get_right_endpoint() {
            return right_endpoint; 
        }
        Element* first_active_element();
        Element* last_active_element();
        void set_bc_left_dirichlet(int eq_n, double val);
        void set_bc_right_dirichlet(int eq_n, double val);
        void refine_single_elem(int id, int3 cand);
        void refine_elems(int elem_num, int *id_array, int3 *cand_array);
        void reference_refinement(int start_elem_id, int elem_num);
        void adapt(int norm, double threshold, Mesh *mesh_ref, 
                   double *y_prev, double *y_prev_ref, 
                   double *err_squared_array);
        Mesh *replicate(); 
        void plot(const char* filename); // plots the mesh and polynomial degrees of elements
        void plot_element_error_p(FILE *f[MAX_EQN_NUM], Element *p, Element *e_ref, 
				  double* y_prev, double* y_prev_ref, 
                                  int subdivision = 20); // plots error wrt. reference solution
        void plot_element_error_hp(FILE *f[MAX_EQN_NUM], Element *p, Element *e_ref_left, Element *e_ref_right, 
				   double* y_prev, double* y_prev_ref, 
                                   int subdivision = 20); // plots error wrt. reference solution
                                                          // if ref. refinement was hp-refinement
        void plot_element_error_exact(FILE *f[MAX_EQN_NUM], Element *p, 
				   double* y_prev, exact_sol_type exact_sol,
                                   int subdivision = 20); // plots error wrt. reference solution
                                                          // if ref. refinement was hp-refinement
        void plot_error_est(const char *filename, Mesh* mesh_ref, 
			double* y_prev, double* y_prev_ref, 
                        int subdivision = 100);  // plots error wrt. reference solution
        void plot_error_exact(const char *filename, double* y_prev, exact_sol_type exact_sol, 
                        int subdivision = 100); // plots error wrt. exact solution
        int assign_elem_ids();
        double bc_left_dir_values[MAX_EQN_NUM];  // values for the Dirichlet condition left
        double bc_right_dir_values[MAX_EQN_NUM]; // values for the Dirichlet condition right

    private:
        double left_endpoint, right_endpoint;
        int n_eq;
        int n_base_elem;
        int n_active_elem;
        int n_dof;
        Element *base_elems;

};

// transforms point 'x_phys' from element (x1, x2) to (-1, 1)
double inverse_map(double x1, double x2, double x_phys);


#endif
