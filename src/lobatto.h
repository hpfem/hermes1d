// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef SHAPESET_LOBATTO_H_
#define SHAPESET_LOBATTO_H_

#include <math.h>

#include "common.h"

// Legendre polynomials
#define legendre0(x) (1.0)
#define legendre1(x) (x)
#define legendre2(x) (1.0 / 2.0 * (3 * (x) * (x) - 1))
#define legendre3(x) (1.0 / 2.0 * (5 * (x) * (x) - 3) * (x))
#define legendre4(x) (1.0 / 8.0 * ((35 * (x) * (x) - 30) * (x) * (x) + 3))
#define legendre5(x) (1.0 / 8.0 * ((63 * (x) * (x) - 70) * (x) * (x) + 15) * (x))
#define legendre6(x) (1.0 / 16.0 * (((231 * (x) * (x) - 315) * (x) * (x) + 105) * (x) * (x) - 5))
#define legendre7(x) (1.0 / 16.0 * (((429 * (x) * (x) - 693) * (x) * (x) + 315) * (x) * (x) - 35) * (x))
#define legendre8(x) (1.0 / 128.0 * ((((6435 * (x) * (x) - 12012) * (x) * (x) + 6930) * (x) * (x) - 1260) * (x) * (x) + 35))
#define legendre9(x) (1.0 / 128.0 * ((((12155 * (x) * (x) - 25740) * (x) * (x) + 18018) * (x) * (x) - 4620) * (x) * (x) + 315) * (x))
#define legendre10(x) (1.0 / 256.0 * (((((46189 * (x) * (x) - 109395) * (x) * (x) + 90090) * (x) * (x) - 30030) * (x) * (x) + 3465) * (x) * (x) - 63))
#define legendre11(x) (1.0 / 256.0 * (((((88179 * (x) * (x) - 230945) * (x) * (x) + 218790) * (x) * (x) - 90090) * (x) * (x) - 15015) * (x) * (x) - 693) * (x))

// Derivatives of Legendre polynomials

#define legendre0x(x) (0.0)
#define legendre1x(x) (1.0)
#define legendre2x(x) (3.0 * (x))
#define legendre3x(x) (15.0 / 2.0 * (x) * (x) - 3.0 / 2.0)
#define legendre4x(x) (5.0 / 2.0 * (x) * (7.0 * (x) * (x) - 3.0))
#define legendre5x(x) ((315.0 / 8.0 * (x) * (x) - 105.0 / 4.0) * (x) * (x) + 15.0 / 8.0)
#define legendre6x(x) (21.0 / 8.0 * (x) * ((33.0 * (x) * (x) - 30.0) * (x) * (x) + 5.0))
#define legendre7x(x) (((3003.0 / 16.0 * (x) * (x) - 3465.0 / 16.0) * (x) * (x) + 945.0 / 16.0) * (x) * (x) - 35.0 / 16.0)
#define legendre8x(x) (9.0 / 16.0 * (x) * (((715.0 * (x) * (x) - 1001.0) * (x) * (x) + 385.0) * (x) * (x) - 35.0))
#define legendre9x(x) ((((109395.0 / 128.0 * (x) * (x) - 45045.0 / 32.0) * (x) * (x) + 45045.0 / 64.0) * (x) * (x) - 3465.0 / 32.0) * (x) * (x) + 315.0 / 128.0)
#define legendre10x(x) (1.0 / 256.0 * ((((461890 * (x) * (x) - 875160) * (x) * (x) + 540540) * (x) * (x) - 120120) * (x) * (x) + 6930) * (x))
#define legendre11x(x) (1.0 / 256.0 * (((((969969 * (x) * (x) - 2078505) * (x) * (x) + 1531530) * (x) * (x) - 450450) * (x) * (x) + 45045) * (x) * (x) - 693))


//

extern shape_fn_t lobatto_fn_tab_1d[];
extern shape_fn_t lobatto_der_tab_1d[];
extern shape_fn_t legendre_fn_tab_1d[];
extern shape_fn_t legendre_der_tab_1d[];

extern int lobatto_order_1d[];
extern int legendre_order_1d[];

#endif /* SHAPESET_LOBATTO_H_ */
