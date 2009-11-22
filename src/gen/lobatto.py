#! /usr/bin/env python

print "Importing..."
import os
from jinja2 import Environment, FileSystemLoader
from sympy import Symbol, ccode

def ccode2(s):
    s = ccode(s)
    s = s.replace("pow(x,2)", "((x)*(x))")
    return s


from common import lobatto, horner_scheme

n_functions = 30
precision = 25
factor_const = False

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
        "expr": ccode2(lob),
        "expr_diff": ccode2(lob_diff),
        })

print "Generating the C file..."
template = "lobatto.cpp"
t = env.get_template(template)
open(os.path.join("..", template), "w").write(t.render({
    "functions": functions,
    }))
