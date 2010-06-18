# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

cdef extern from "math.h":

    double c_sqr "sqr"(double x)
    double c_sqrt "sqrt"(double x)
    double c_atan "atan"(double x)
    double c_pi "M_PI"

cdef extern from "stdlib.h":

    ctypedef unsigned long size_t
    void *malloc (size_t size)
    void free(void *mem)
    void *memcpy(void *dst, void *src, long n)

    void exit(int exit_code)

cdef extern from "hermes1d.h":

    ctypedef double double4[4]
    ctypedef double double3[3]
    ctypedef int int3[3]
    ctypedef int int2[2]

    cdef struct c_Element "Element":
        double x1, x2
        int p
        int *dof

    cdef struct c_Mesh "Mesh":
        void create(double A, double B, int n)
        int get_n_base_elems()
        int get_n_dofs()
        void set_poly_orders(int poly_order)
        void assign_dofs()
        c_Element *get_base_elems()
        void set_dirichlet_bc_left(int eq_n, double val)
        void set_dirichlet_bc_right(int eq_n, double val)
    c_Mesh *new_Mesh "new Mesh" (double a, double b, int n_elem, int p_init,
            int eq_num)

    cdef struct c_Linearizer "Linearizer":
        void plot_solution(char *out_filename,
                int plotting_elem_subdivision)
        void get_xy_mesh(int comp, int plotting_elem_subdivision,
                double **x, double **y, int *n)
    c_Linearizer *new_Linearizer "new Linearizer" (c_Mesh *mesh)
