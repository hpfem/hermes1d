#! /usr/bin/env python

print "Importing..."
import os
from jinja2 import Environment, FileSystemLoader
from sympy import Symbol, integrate, legendre, factor, sqrt, ccode

n_functions = 100
precision = 25
factor_const = False

def legendre_int(i, x):
    """
    Returns the normalized integrated Legendre polynomial.
    """
    y = Symbol("y", dummy=True)
    f = legendre(i, y)
    n = sqrt(integrate(f**2, (y, -1, 1)))
    return integrate(f, (y, -1, x))/n

def lobatto(i, x):
    """
    Returns the Lobatto shape function.
    """
    if i == 0:
        return 1-lobatto(1, x)
    f = legendre_int(i-1, x)
    if i == 1:
        f /= sqrt(2)
    return f.expand()

def horner_scheme(p, x, factor_const=True):
    """
    Rewrites the polynomial using the Horner scheme.
    """
    a = p.subs(x, 0)
    if p == a:
        return p
    rest = ((p-a)/x).expand()
    if factor_const:
        const = rest.subs(x, 0)
        if const != 0:
            rest = (rest/const).expand()
        else:
            const = 1
    else:
        const = 1
    return const*x*horner_scheme(rest, x, factor_const)+a

x = Symbol("x")
env = Environment(loader=FileSystemLoader('.'))

functions = []
print "Calculating shape functions..."
for i in range(n_functions):
    print "  i=%d" % i
    lob = lobatto(i, x)
    lob_diff = lob.diff(x)
    lob = horner_scheme(lob.n(precision), x, factor_const=factor_const)
    print lob
    lob_diff = horner_scheme(lob_diff.n(precision), x,
            factor_const=factor_const)
    functions.append({"id": i,
        "expr": ccode(lob),
        "expr_diff": ccode(lob_diff),
        })

print "Generating the C file..."
template = "lobatto.cpp"
t = env.get_template(template)
open(os.path.join("..", template), "w").write(t.render({
    "functions": functions,
    }))
