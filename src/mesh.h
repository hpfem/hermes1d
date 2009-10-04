#ifndef _MESH_H_
#define _MESH_H_

#include "common.h"
#include "lobatto.h"

class Mesh {
    public:
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

    private:
        int n_elem;
        int n_dof;
        Vertex *vertices;
        Element *elems;

    // Dirichlet boundary conditions at both endpoints
    // (first integer in pair indicates whether there is 
    // a Dirichlet condition for that equation, the other
    // one tells the value)
    int2 DIR_BC_LEFT[];
    int2 DIR_BC_RIGHT[];

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
