// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "transforms.h"

#define N_chebyshev (3+2*(P_MAX-1))

double chebyshev_points[N_chebyshev]

void calculate_chebyshev_points()
{
    for (int i; i < N_chebyshev, i++)
        chebyshev_points[i] = cos(i*M_PI/(N_chebyshev-1));
}

double phi(int i, double x)
{
}
