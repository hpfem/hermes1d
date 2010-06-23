"""
Module for handling Fekete points approximations.
"""

from gauss_lobatto_points import points

def get_x_phys(x_ref, a, b):
    return (a+b)/2. + x_ref*(b-a)/2.;

class Mesh1D(object):

    def __init__(self, points, orders):
        if not (len(points) == len(orders) + 1):
            raise Exception("points vs order mismatch")
        self._points = points
        self._orders = orders

    def iter_elems(self):
        for i in range(len(self._orders)):
            yield (self._points[i], self._points[i+1], self._orders[i])

class Function(object):
    """
    Represents a function on a mesh.

    The values are given in the Fekete points.
    """

    def __init__(self, obj, mesh=None):
        if isinstance(obj, Function):
            raise NotImplementedError()
        else:
            if not isinstance(mesh, Mesh1D):
                raise Exception("You need to specify a mesh.")
            self._mesh = mesh
            self._values = []
            for a, b, order in mesh.iter_elems():
                fekete_points = points[order]
                elem_values = []
                for p in fekete_points:
                    p = get_x_phys(p, a, b)
                    val = obj(p)
                    elem_values.append(val)
                self._values.append(elem_values)
            print self._values

    def __call__(self, x):
        print x
        return 5



def test1():
    m = Mesh1D((-5, -4, 3, 10), (1, 5, 1))

def test2():
    eps = 1e-12
    func = lambda x: x**2
    f = Function(func, Mesh1D((-5, -4, 3, 10), (1, 5, 1)))
    for x in [-5, -4, 3, 10, 0, 0.01, 0.0001]:
        assert abs(f(x) - func(x)) < eps

def main():
    test1()
    test2()

if __name__ == "__main__":
    main()
