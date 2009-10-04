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

    def __cinit__(self, int eq_num):
        self.thisptr = new_Mesh(eq_num)

    def __dealloc__(self):
        delete(self.thisptr)


import sys
import traceback

global_namespace = {"verbose": False}

cdef api void cmd(char *text):
    n = run_cmd(text, global_namespace)
    global_namespace.update(n)

cdef api void set_verbose_cmd(int verbose):
    global_namespace["verbose"] = verbose

cdef api void insert_int(char *name, int i):
    """
    Inserts the int "i" into the global namespace.

    Example:

    insert_int("a", 34);
    cmd("print a");

    This prints "34".
    """
    global_namespace.update({name: i})

cdef api void insert_double_array(char *name, double *A, int len):
    """
    Inserts an array of doubles into the global namespace as a NumPy array.

    Example:

    double a[3] = {1, 5, 3};
    insert_double_array("A", a, 3);
    cmd("print A");

    This prints "[ 1.  5.  3.]" (this is how the NumPy array is printed).
    """
    global_namespace.update({name: array_double_c2numpy(A, len)})

cdef api void insert_int_array(char *name, int *A, int len):
    """
    Inserts an array of ints into the global namespace as a NumPy array.

    Example:

    int a[3] = {1, 5, 3};
    insert_double_array("A", a, 3);
    cmd("print A");

    This prints "[1  5  3]" (this is how the NumPy array is printed).
    """
    global_namespace.update({name: array_int_c2numpy(A, len)})

cdef api void insert_object(char *name, object o):
    global_namespace.update({name: o})

cdef api object get_symbol(char *name):
    return global_namespace.get(name)

cdef ndarray array_int_c2numpy(int *A, int len):
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="int32")
    cdef int *pvec = <int *>vec.data
    memcpy(pvec, A, len*sizeof(int))
    return vec

cdef ndarray array_double_c2numpy(double *A, int len):
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="double")
    cdef double *pvec = <double *>vec.data
    memcpy(pvec, A, len*sizeof(double))
    return vec

cdef api void array_double_numpy2c_inplace(object A_n, double **A_c, int *n):
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
