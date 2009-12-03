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
  // maximum poly degree of Lobatto function tested
  int max_test_poly_degree = MAX_P;
  int ok = 1;

  // maximum allowed error at an integration point
  double max_allowed_error = 1e-12;

  // loop over Lobatto shape functions starting with
  // the quadratic one
  for (int n = 2; n < max_test_poly_degree; n++) {
    // looking at the difference at integration points using 
    // Gauss quadratures of orders 1, 2, ... MAX_QUAD_ORDER
    for (int quad_order=0; quad_order < MAX_QUAD_ORDER; quad_order++) {
      int num_pts = g_quad_1d_std.get_num_points(quad_order);
      double2 *quad_tab = g_quad_1d_std.get_points(quad_order);
      for (int i=0; i<num_pts; i++) {
        double point_i = quad_tab[i][0];
        double val = fabs(calc_legendre_val(point_i, n-1) -
                          calc_lobatto_der(point_i, n));
        printf("poly_deg = %d, quad_order = %d, x = %g, difference = %g\n", 
               n, quad_order, point_i, val);
        if(val > max_allowed_error) {
          printf("Failure!\n");
          return ERROR_FAILURE;
        }
      }
    }
  }

  printf("Success!\n");
  return ERROR_SUCCESS;
}
