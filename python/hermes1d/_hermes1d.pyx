# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

from _hermes_common cimport c2numpy_double, delete

cdef class Element:
    cdef c_Element *thisptr

    @property
    def p(self):
        return self.thisptr.p

    @property
    def x1(self):
        return self.thisptr.x1

    @property
    def x2(self):
        return self.thisptr.x2

cdef class Mesh:
    cdef c_Mesh *thisptr

    def __init__(self, double a, double b, int n_elem, int p_init, int eq_num):
        self.thisptr = new_Mesh(a, b, n_elem, p_init, eq_num)

    def __dealloc__(self):
        delete(self.thisptr)

cdef class Linearizer:
    cdef c_Linearizer *thisptr

    def __cinit__(self, Mesh mesh):
        self.thisptr = new_Linearizer(mesh.thisptr)

    def plot_solution(self, out_filename, plotting_elem_subdivision):
        #cdef double *A
        cdef int n
        #numpy2c_double_inplace(y_prev, &A, &n)
        self.thisptr.plot_solution(out_filename, plotting_elem_subdivision)

    def get_xy(self, int comp, int plotting_elem_subdivision):
        """
        Returns (x, y), where x, y are arrays of points.

        y_prev is the input array of points.
        """
        #cdef double *A
        #cdef int nA
        #numpy2c_double_inplace(y_prev, &A, &nA)
        cdef double *x
        cdef double *y
        cdef int n
        self.thisptr.get_xy_mesh(comp, plotting_elem_subdivision,
                &x, &y, &n)
        x_numpy = c2numpy_double(x, n)
        y_numpy = c2numpy_double(y, n)
        return x_numpy, y_numpy
