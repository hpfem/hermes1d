# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

from hermes_common._hermes_common cimport c2numpy_double, delete, PY_NEW, \
    numpy2c_double_inplace, numpy2c_int_inplace

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

    def __init__(self, *args):
        from numpy import array
        cdef double a, b
        cdef int n_elem, p_init, eq_num, n
        cdef double *pts_array
        cdef int *p_array, *m_array, *div_array
        if len(args) == 5:
            a, b, n_elem, p_init, eq_num = args
            self.thisptr = new_Mesh(a, b, n_elem, p_init, eq_num)
        elif len(args) == 4:
            a, b, n_elem, p_init = args
            eq_num = 1
            self.thisptr = new_Mesh(a, b, n_elem, p_init, eq_num)
        elif len(args) == 2:
            pts, p = args
            pts = array(pts, dtype="double")
            p = array(p, dtype="int32")
            if not (len(pts) == len(p) + 1):
                raise ValueError("len(pts) must be equal to len(p) + 1")
            n_elem = len(p)
            m = array([1]*len(p), dtype="int32")
            div = array([1]*len(p), dtype="int32")
            eq_num = 1
            numpy2c_double_inplace(pts, &pts_array, &n)
            numpy2c_int_inplace(p, &p_array, &n)
            numpy2c_int_inplace(m, &m_array, &n)
            numpy2c_int_inplace(div, &div_array, &n)
            self.thisptr = new_Mesh2(n_elem, pts_array, p_array, m_array, div_array,
                    eq_num)
        else:
            raise ValueError("Don't understand the arguments")

    def __dealloc__(self):
        delete(self.thisptr)

    def copy_vector_to_mesh(self, sol, int comp):
        cdef double *Y
        cdef int n
        numpy2c_double_inplace(sol, &Y, &n)
        self.thisptr.copy_vector_to_mesh(Y, comp)

    def assign_dofs(self):
        return self.thisptr.assign_dofs()

    def plot_to_file(self, filename):
        self.thisptr.plot(filename)

cdef class Linearizer:
    cdef c_Linearizer *thisptr
    cdef Mesh mesh

    def __cinit__(self, Mesh mesh):
        self.thisptr = new_Linearizer(mesh.thisptr)
        self.mesh = mesh

    def plot_solution(self, out_filename, plotting_elem_subdivision):
        #cdef double *A
        cdef int n
        #numpy2c_double_inplace(y_prev, &A, &n)
        self.thisptr.plot_solution(out_filename, plotting_elem_subdivision)

    def get_xy(self, sol, int comp, int plotting_elem_subdivision):
        """
        Returns (x, y), where x, y are arrays of points.

        sol ... is the solution to linearize
        comp ... solution component
        """
        cdef double *x
        cdef double *y
        cdef int n
        self.mesh.copy_vector_to_mesh(sol, comp)
        self.thisptr.get_xy_mesh(comp, plotting_elem_subdivision,
                &x, &y, &n)
        x_numpy = c2numpy_double(x, n)
        y_numpy = c2numpy_double(y, n)
        return x_numpy, y_numpy

cdef api object c2py_Mesh(c_Mesh *h):
    cdef Mesh n
    n = <Mesh>PY_NEW(Mesh)
    n.thisptr = h
    return n
