#ifndef _MESH_H_
#define _MESH_H_

#include "common.h"
#include "lobatto.h"

class Mesh {
    public:
        Mesh(int n_eqn) {
            this->n_eqn = n_eqn;
            this->dir_bc_left_active = new int[n_eqn];
            this->dir_bc_left_values = new double[n_eqn];
            this->dir_bc_right_active = new int[n_eqn];
            this->dir_bc_right_values = new double[n_eqn];
            for (int i=0; i<n_eqn; i++) {
                this->dir_bc_left_active[i] = 0;
                this->dir_bc_left_values[i] = 0;
                this->dir_bc_right_active[i] = 0;
                this->dir_bc_right_values[i] = 0;
            }
        }
        void create(double A, double B, int n);
        void set_poly_orders(int poly_order);
        void assign_dofs();
        Vertex *get_vertices() {
            return this->vertices;
        }
        Element *get_elems() {
            return this->elems;
        }
        int get_n_elems() {
            return this->n_elem;
        }
        int get_n_dof() {
            return this->n_dof;
        }
        void set_dirichlet_bc_left(int eq_n, double val);
        void set_dirichlet_bc_right(int eq_n, double val);

    private:
        int n_eqn;
        int n_elem;
        int n_dof;
        Vertex *vertices;
        Element *elems;

        // Dirichlet boundary conditions at both endpoints
        // (first integer in pair indicates whether there is 
        // a Dirichlet condition for that equation, the other
        // one tells the value)
        int *dir_bc_left_active;
        double *dir_bc_left_values;
        int *dir_bc_right_active;
        double *dir_bc_right_values;

};

class Linearizer {
    public:
        Linearizer(Mesh *mesh) {
            this->mesh = mesh;
        }
        // evaluate approximate solution at element 'm' at reference 
        // point 'x_ref'. Here 'y' is the global vector of coefficients
        void eval_approx(Element *e, double x_ref, double *y, double &x_phys,
                double &val);
        void plot_solution(const char *out_filename, double *y_prev, int
                plotting_elem_subdivision=100);

    private:
        Mesh *mesh;
};

#endif
