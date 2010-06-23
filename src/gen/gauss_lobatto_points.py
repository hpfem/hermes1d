# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

# AUTOGENERATED BY gauss_lobatto.py, DO NOT EDIT

points = {
    {% for s in data %}
    "{{ s.p }}": [
        {% for point in s.points %}{{ point }},
        {% endfor %} ],
    {% endfor %}
}