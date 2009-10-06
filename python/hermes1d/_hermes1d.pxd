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

cdef extern from "arrayobject.h":

    cdef enum NPY_TYPES:
        NPY_BOOL
        NPY_BYTE
        NPY_UBYTE
        NPY_SHORT
        NPY_USHORT
        NPY_INT
        NPY_UINT
        NPY_LONG
        NPY_ULONG
        NPY_LONGLONG
        NPY_ULONGLONG
        NPY_FLOAT
        NPY_DOUBLE
        NPY_LONGDOUBLE
        NPY_CFLOAT
        NPY_CDOUBLE
        NPY_CLONGDOUBLE
        NPY_OBJECT
        NPY_STRING
        NPY_UNICODE
        NPY_VOID
        NPY_NTYPES
        NPY_NOTYPE

    ctypedef int npy_intp

    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef npy_intp *dimensions
        cdef npy_intp *strides
        cdef int flags

    object PyArray_SimpleNewFromData(int nd, npy_intp* dims, int typenum,
            void* data)
    void import_array()


cdef extern from "Python.h":
    ctypedef void PyObject
    void Py_INCREF(PyObject *x)
    void Py_DECREF(PyObject *x)

cdef extern from "stdcython.h":
    void init_global_empty_tuple()
    object PY_NEW(object t)

cdef extern from "hermes1d.h":

    # This is just the C++ "delete" statement
    void delete(...)

    void throw_exception(char *msg)

    ctypedef double double4[4]
    ctypedef double double3[3]
    ctypedef int int3[3]
    ctypedef int int2[2]

    cdef struct c_Vertex "Vertex":
        double x

    cdef struct c_Element "Element":
        c_Vertex* *v1, *v2
        int p
        int *dof

    cdef struct c_Mesh "Mesh":
        void create(double A, double B, int n)
        int get_n_elems()
        int get_n_dofs()
        void set_poly_orders(int poly_order)
        void assign_dofs()
        c_Vertex *get_vertices()
        c_Element *get_elems()
        void set_dirichlet_bc_left(int eq_n, double val)
        void set_dirichlet_bc_right(int eq_n, double val)
    c_Mesh *new_Mesh "new Mesh" (int eq_num)
