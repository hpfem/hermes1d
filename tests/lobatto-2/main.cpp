#include "hermes1d.h"

#include "legendre.h"
#include "lobatto.h"
#include "quad_std.h"

// This test makes sure that the derivatives of 
// the Lobatto shape functions starting with the 
// quadratic one are the Legendre polynomials
// at all possible quadrature points

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int main(int argc, char* argv[])
{
  // maximum index of Lobatto function tested
  int max_n = 30; //FIXME: should be 99

  // maximum allowed error at an integration point
  double max_allowed_error = 1e-12;

  // loop over Lobatto shape functions starting with
  // the quadratic one
  for (int n = 2; n < max_n; n++) {
    // looking at the difference at integration points using 
    // Gauss quadratures of orders 1, 2, ... max_n
    for (int quad_order=0; quad_order < max_n; quad_order++) {
      int num_pts = g_quad_1d_std.get_num_points(quad_order);
      double2 *quad_tab = g_quad_1d_std.get_points(quad_order);
      for (int i=0; i<num_pts; i++) {
        double point_i = quad_tab[i][0];
        double val = fabs(legendre_fn_tab_1d[n-1](point_i) -
                          lobatto_der_tab_1d[n](point_i));
        //printf("%g %g\n", legendre_fn_tab_1d[n-1](point_i), 
        //       lobatto_der_tab_1d[n](point_i));
        if(val > max_allowed_error) {
          printf("n = %d, quad_order = %d, i = %d, difference = %g\n", 
                 n, quad_order, i, val);
          return ERROR_FAILURE;
        }
      }
    }
  }

  printf("Success!\n");
  return ERROR_SUCCESS;
}
