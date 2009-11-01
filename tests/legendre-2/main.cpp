#include "hermes1d.h"

#include "legendre.h"
#include "quad_std.h"

// This test makes sure that Legendre polynomials
// are orthonormal.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int main(int argc, char* argv[])
{
  // maximum poly degree of Legendre polynomials tested
  int max_poly_degree = 20; 

  // maximum allowed error
  double max_allowed_error = 1e-10;

  // loop over polynomial degrees, starting with 1
  double max_actual_error = 0; 
  for (int poly_deg_1=0; poly_deg_1<max_poly_degree; poly_deg_1++) {
    for (int poly_deg_2=0; poly_deg_2<max_poly_degree; poly_deg_2++) {
      // integrating the Legendre polynomial of degree 'poly_deg'
      // from -1 to 1 using Gauss quadratures of orders 1, 2, ...
      // MAX_POLYORDER
      for (int quad_order=poly_deg_1 + poly_deg_2; 
           quad_order<2*max_poly_degree; quad_order++) {
        int num_pts = g_quad_1d_std.get_num_points(quad_order);
        double2 *quad_tab = g_quad_1d_std.get_points(quad_order);
        double val = 0;
        for (int i=0; i<num_pts; i++) {
          double point_i = quad_tab[i][0];
          double weight_i = quad_tab[i][1];
          val += legendre_fn_tab_1d[poly_deg_1](point_i) * 
                 legendre_fn_tab_1d[poly_deg_2](point_i) * weight_i;
        }
        if (poly_deg_1 == poly_deg_2) { // val must be one
          if (fabs(val - 1.) > max_actual_error) 
            max_actual_error = fabs(val - 1.);
          if (max_actual_error > max_allowed_error) {
            printf("poly_deg_1 = %d, poly_deg_2 = %d, quad_order = %d, val - 1. = %g\n", 
                   poly_deg_1, poly_deg_2, quad_order, val - 1.);      
            return ERROR_FAILURE;
          }
        }
        else { // val must be zero
          if (fabs(val) > max_actual_error) 
            max_actual_error = fabs(val);
          if (max_actual_error > max_allowed_error) {
            printf("poly_deg_1 = %d, poly_deg_2 = %d, quad_order = %d, val = %g\n", 
                   poly_deg_1, poly_deg_2, quad_order, val);      
            return ERROR_FAILURE;
          }
        }
      }
    }
  }

  printf("Success!\n");
  return ERROR_SUCCESS;
}
