#ifndef _MESH_H_
#define _MESH_H_

#include "common.h"

class Mesh {
    public:
        void create(double A, double B, int n);
        void set_poly_orders(int poly_order);

    private:
        int n_elem;
        Vertex *vertices;
        Element *elems;
};

#endif
