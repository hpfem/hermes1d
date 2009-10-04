#ifndef _MESH_H_
#define _MESH_H_

#include "common.h"

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
};

#endif
