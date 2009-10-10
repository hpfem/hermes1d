cdef class Vertex:
    cdef c_Vertex *thisptr

    @property
    def x(self):
        return self.thisptr.x

    def __str__(self):
        return "Vertex(%r)" % (self.x)

cdef class Element:
    cdef c_Element *thisptr

    @property
    def p(self):
        return self.thisptr.p

cdef class Mesh:
    cdef c_Mesh *thisptr

    def __init__(self, int eq_num):
        self.thisptr = new_Mesh(eq_num)

    def __dealloc__(self):
        delete(self.thisptr)

cdef class Linearizer:
    cdef c_Linearizer *thisptr

    def __cinit__(self, Mesh mesh):
        self.thisptr = new_Linearizer(mesh.thisptr)

    def plot_solution(self, out_filename, y_prev, plotting_elem_subdivision):
        cdef double *A
        cdef int n
        numpy2c_double_inplace(y_prev, &A, &n)
        self.thisptr.plot_solution(out_filename, A, plotting_elem_subdivision)

    def get_xy(self, y_prev, int comp, int plotting_elem_subdivision):
        """
        Returns (x, y), where x, y are arrays of points.

        y_prev is the input array of points.
        """
        cdef double *A
        cdef int nA
        numpy2c_double_inplace(y_prev, &A, &nA)
        cdef double *x
        cdef double *y
        cdef int n
        self.thisptr.get_xy(A, comp, plotting_elem_subdivision,
                &x, &y, &n)
        x_numpy = c2numpy_double(x, n)
        y_numpy = c2numpy_double(y, n)
        return x_numpy, y_numpy

#-----------------------------------------------------------------------
# Common C++ <-> Python+NumPy conversion tools:

import sys
import traceback
# this is important to be called here, otherwise we can't use the NumPy C/API:
import_array()

global_namespace = {"verbose": False}

cdef api void cmd(char *text):
    n = run_cmd(text, global_namespace)
    global_namespace.update(n)

cdef api void set_verbose_cmd(int verbose):
    global_namespace["verbose"] = verbose

cdef api void insert_object(char *name, object o):
    """
    Inserts an object into the global namespace.

    Example 1:

    insert_object("a", c2py_int(3));
    cmd("print a");

    This prints "3".

    Example 2:

    int a[3] = {1, 5, 3};
    insert_object("A", c2numpy_int(a, 3));
    cmd("print A");

    This prints "[1  5  3]" (this is how the NumPy array is printed).

    Example 3:

    double a[3] = {1, 5, 3};
    insert_object("A", c2numpy_double(a, 3));
    cmd("print A");

    This prints "[ 1.  5.  3.]" (this is how the NumPy array is printed).
    """
    global_namespace.update({name: o})

cdef api object get_object(char *name):
    """
    Retrieves an object from the Python namespace.

    Example:

    // put into python:
    int a[3] = {1, 5, 3};
    insert_object("A", c2numpy_int(a, 3));

    // retrieve from python:
    double *A;
    int n;
    numpy2c_double_inplace(get_object("A"), &A, &n);
    """
    return global_namespace.get(name)

cdef api object c2py_int(int i):
    return i

cdef api int py2c_int(object i):
    return i

cdef api double py2c_double(object i):
    return i

cdef api object c2py_mesh(c_Mesh *m):
    cdef Mesh mesh
    print "ok1"
    mesh = <Mesh>PY_NEW(Mesh)
    print "ok2"
    mesh.thisptr = <c_Mesh *>m
    print "ok3"
    return mesh

cdef api object c2numpy_int(int *A, int len):
    """
    Construct the integer NumPy array by copying the data.
    """
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="int32")
    cdef int *pvec = <int *>vec.data
    memcpy(pvec, A, len*sizeof(int))
    return vec

cdef api object c2numpy_int_inplace(int *A, int len):
    """
    Construct the integer NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_INT, A)

cdef api object c2numpy_double(double *A, int len):
    """
    Construct the double NumPy array by copying the data.
    """
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="double")
    cdef double *pvec = <double *>vec.data
    memcpy(pvec, A, len*sizeof(double))
    return vec

cdef api object c2numpy_double_inplace(double *A, int len):
    """
    Construct the double NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, A)

cdef api void numpy2c_double_inplace(object A_n, double **A_c, int *n):
    cdef ndarray A = A_n
    if not (A.nd == 1 and A.strides[0] == sizeof(double)):
        from numpy import array
        A = array(A.flat, dtype="double")
    n[0] = len(A)
    A_c[0] = <double *>(A.data)

cdef api object run_cmd(char *text, object namespace):
    try:
        verbose = namespace.get("verbose")
        if verbose:
            print "got a text:", text
        if verbose:
            print "evaluting in the namespace:"
            print namespace
        code = compile(text, "", "exec")
        eval(code, {}, namespace)
        if verbose:
            print "new namespace:"
            print namespace
        return namespace
    except SystemExit, e:
        try:
            exit_code = int(e)
        except:
            exit_code = -1
        exit(exit_code)
    except:
        etype, value, tb = sys.exc_info()
        s = "".join(traceback.format_exception(etype, value, tb))
        s = "Exception raised in the Python code:\n" + s
        throw_exception(s)
