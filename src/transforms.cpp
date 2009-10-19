// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "transforms.h"

#define N_chebyshev (3+2*(MAX_P-1))

double chebyshev_points[N_chebyshev];

void calculate_chebyshev_points()
{
    for (int i; i < N_chebyshev; i++)
        chebyshev_points[i] = cos(i*M_PI/(N_chebyshev-1));
}

// transform values from (-1, 0) to (-1, 1)
#define map_left(x) (2*x+1)
// transform values from (0, 1) to (-1, 1)
#define map_right(x) (2*x-1)

double phi(int i, double x)
{
    if (x < 0) {
        if (i == 0)
            return lobatto_fn_tab_1d[0](map_left(x));
        else if (i % 2 == 0)
            return 0;
        else
            return lobatto_fn_tab_1d[(i+1)/2](map_left(x));
    } else {
        if (i == 0)
            return 0;
        else if (i == 1)
            return lobatto_fn_tab_1d[0](map_left(x));
        else if (i % 2 == 1)
            return 0;
        else
            return lobatto_fn_tab_1d[i/2](map_left(x));
    } 
}
