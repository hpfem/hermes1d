// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "mesh.h"
#include "linearizer.h"
#include "iterator.h"

// Evaluate (vector-valued) approximate solution at reference 
// point 'x_ref' in element 'm'. Here 'y' is the global vector 
// of coefficients. The result is a vector of length mesh->n_eq
void Linearizer::eval_approx(Element *e, double x_ref, double *y, 
                             double *x_phys, double *val) {
  int n_eq = this->mesh->get_n_eq();
  for(int c=0; c<n_eq; c++) { // loop over solution components
    val[c] = 0;
    for(int i=0; i <= e->p; i++) { // loop over shape functions
      if(e->dof[c][i] >= 0) val[c] += y[e->dof[c][i]]*lobatto_fn_tab_1d[i](x_ref);
    }
  }
  double a = e->x1;
  double b = e->x2;
  *x_phys = (a+b)/2 + x_ref*(b-a)/2;
  return;
}

// Plot solution in Gnuplot format
void Linearizer::plot_solution(const char *out_filename, 
                               double *y_prev, int plotting_elem_subdivision)
{
    int n_eq = this->mesh->get_n_eq();
    FILE *f[MAX_EQN_NUM];
    char final_filename[MAX_EQN_NUM][MAX_STRING_LENGTH];
    for(int c=0; c<n_eq; c++) {
        if(n_eq == 1)
            sprintf(final_filename[c], "%s", out_filename);
        else
            sprintf(final_filename[c], "%s_%d", out_filename, c);
        f[c] = fopen(final_filename[c], "wb");
        if(f[c] == NULL) error("problem opening file in plot_solution().");
        int n;
        double *x, *y;
        this->get_xy(y_prev, c, plotting_elem_subdivision, &x, &y, &n);
        for (int i=0; i < n; i++)
            fprintf(f[c], "%g %g\n", x[i], y[i]);
        delete[] x;
        delete[] y;
        printf("Output written to %s.\n", final_filename[c]);
        fclose(f[c]);
    }
}

// Returns pointers to x and y coordinates in **x and **y
// you should free it yourself when you don't need it anymore
// y_prev --- the solution coefficients (all equations)
// comp --- which component you want to process
// plotting_elem_subdivision --- the number of subdivision of the element
// x, y --- the doubles list of x,y
// n --- the number of points

void Linearizer::get_xy(double *y_prev, int comp,
        int plotting_elem_subdivision,
        double **x, double **y, int *n)
{
    int n_eq = this->mesh->get_n_eq();
    int n_active_elem = this->mesh->get_n_active_elem();
    Iterator *I = new Iterator(this->mesh);

    *n = n_active_elem * (plotting_elem_subdivision+1);
    double *x_out = new double[*n];
    double *y_out = new double[*n];

    // FIXME:
    if(n_eq > MAX_EQN_NUM) {
      printf("n_eq = %d\n", n_eq);
        error("number of equations too high in plot_solution().");
    }
    // FIXME
    if(plotting_elem_subdivision > MAX_PTS_NUM)
        error("plotting_elem_subdivision too high in plot_solution().");
    double phys_u_prev[MAX_EQN_NUM][MAX_PTS_NUM];
    double phys_du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM];
        
    Element *e;
    int counter = 0;
    while ((e = I->next_active_element()) != NULL) {
        if (counter >= n_active_elem)
            error("Internal error: wrong n_active_elem");
        // FIXME:
        if(e->p > MAX_POLYORDER)
            error("element degree too hign in plot(solution).");
        double coeffs[MAX_EQN_NUM][MAX_COEFFS_NUM];
        e->get_coeffs(y_prev, coeffs, 
                      this->mesh->bc_left_dir_values,
                      this->mesh->bc_right_dir_values);

        double pts_array[MAX_PTS_NUM];
        double h = 2./plotting_elem_subdivision;

        for (int j=0; j<plotting_elem_subdivision+1; j++)
            pts_array[j] = -1 + j*h;
        e->get_solution(coeffs,
                plotting_elem_subdivision+1, pts_array,
                phys_u_prev, phys_du_prevdx);
        double a = e->x1;
        double b = e->x2;
        for (int j=0; j<plotting_elem_subdivision+1; j++) {
            x_out[counter*(plotting_elem_subdivision+1) + j] =
                (a + b)/2 + pts_array[j] * (b-a)/2;
            y_out[counter*(plotting_elem_subdivision+1) + j] =
                phys_u_prev[comp][j];
        }
        counter++;
    }
    *x = x_out;
    *y = y_out;
    delete I;
}

