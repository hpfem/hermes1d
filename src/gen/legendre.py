#! /usr/bin/env python

print "Importing..."
import os
from jinja2 import Environment, FileSystemLoader
from sympy import Symbol, integrate, legendre, factor, sqrt, ccode

n_functions = 100
precision = 25

def legendre_norm(i, x):
    """
    Returns the normalized integrated Legendre polynomial.
    """
    f = legendre(i, x)
    n = sqrt(integrate(f**2, (x, -1, 1)))
    return f/n

def legendre_shape_function(i, x):
    """
    Returns the Lobatto shape function.
    """
    #if i == 0:
    #    return 1-legendre_shape_function(1, x)
    f = legendre_norm(i, x)
    #if i == 1:
    #    f /= sqrt(2)
    return f.expand()

def horner_scheme(p, x):
    """
    Rewrites the polynomial using the Horner scheme.
    """
    a = p.subs(x, 0)
    if p == a:
        return p
    rest = ((p-a)/x).expand()
    return x*horner_scheme(rest, x)+a

x = Symbol("x")
env = Environment(loader=FileSystemLoader('.'))

functions = []
print "Calculating shape functions..."
for i in range(n_functions):
    print "  i=%d" % i
    lob = legendre_shape_function(i, x)
    lob_diff = lob.diff(x)
    lob = horner_scheme(lob.n(precision), x)
    lob_diff = horner_scheme(lob_diff.n(precision), x)
    functions.append({"id": i,
        "expr": ccode(lob),
        "expr_diff": ccode(lob_diff),
        })

print "Generating the C file..."
template = "legendre.cpp"
t = env.get_template(template)
open(os.path.join("..", template), "w").write(t.render({
    "functions": functions,
    }))
