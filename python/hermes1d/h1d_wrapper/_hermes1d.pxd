# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

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
        void copy_vector_to_mesh(double *y, int sln)
        void copy_mesh_to_vector(double *y, int sln)
    c_Mesh *new_Mesh "new Mesh" (double a, double b, int n_elem, int p_init,
            int eq_num)
    c_Mesh *new_Mesh2 "new Mesh" (int n_macro_elem, double *pts_array,
            int *p_array, int *m_array, int *div_array, int n_eq)

    cdef struct c_Linearizer "Linearizer":
        void plot_solution(char *out_filename,
                int plotting_elem_subdivision)
        void get_xy_mesh(int comp, int plotting_elem_subdivision,
                double **x, double **y, int *n)
    c_Linearizer *new_Linearizer "new Linearizer" (c_Mesh *mesh)
