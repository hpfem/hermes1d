#ifndef _MESH_H_
#define _MESH_H_

#include "common.h"
#include "lobatto.h"


class Mesh {
    public:
        Mesh(int n_eq) {
            // Print the banner (only once)
            static int n_calls = 0;
            n_calls++;
            if (n_calls == 1) intro();
            this->n_eq = n_eq;
            this->bc_left_dir = new int[n_eq];
            this->bc_left_dir_values = new double[n_eq];
            this->bc_right_dir = new int[n_eq];
            this->bc_right_dir_values = new double[n_eq];
            for (int i=0; i<n_eq; i++) {
                this->bc_left_dir[i] = BC_INVALID;
                this->bc_left_dir_values[i] = 0;
                this->bc_right_dir[i] = BC_INVALID;
                this->bc_right_dir_values[i] = 0;
            }
        }
        void create(double a, double b, int n_elem);
        void set_poly_orders(int poly_order);
        void assign_dofs(int n_eq);
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
        void set_bc_left_dirichlet(int eq_n, double val);
        void set_bc_left_natural(int eqn);
        void set_bc_right_dirichlet(int eq_n, double val);
        void set_bc_right_natural(int eqn);

        // Dirichlet boundary conditions at both endpoints
        // (first integer in pair indicates whether there is 
        // a Dirichlet condition for that equation, the other
        // one tells the value)
        int *bc_left_dir; //1...Dirichlet (essential) BC at left end point
                          //0...natural BC (Neumann, Newton, none)
        double *bc_left_dir_values; // values for the Dirichlet condition
        int *bc_right_dir; //1...Dirichlet (essential) BC at right end point
                           //0...natural BC (Neumann, Newton, none)
        double *bc_right_dir_values; // values for the Dirichlet condition

    private:
        int n_eq;
        int n_elem;
        int n_dof;
        Vertex *vertices;
        Element *elems;

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
                plotting_elem_subdivision=100);

    private:
        Mesh *mesh;
};

void element_solution(Element *e, double *coeff, int pts_num, 
        double *pts_array, double *val, double *der);

void element_solution_point(double x_ref, Element *e, 
			    double *coeff, double *val, double *der);

#endif
