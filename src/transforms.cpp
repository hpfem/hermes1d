// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "transforms.h"

#define N_chebyshev (3+2*(MAX_P-1))

double chebyshev_points[N_chebyshev];
double chebyshev_matrix[N_chebyshev][N_chebyshev];
double transformation_matrix[N_chebyshev][MAX_P+1];
int transformation_matrix_initialized=0;

void fill_chebyshev_points()
{
    for (int i=0; i < N_chebyshev; i++)
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
        else {
            printf("XXX: %d %f\n", i, x);
            return lobatto_fn_tab_1d[i/2](map_left(x));
        }
    } 
}

void fill_chebyshev_matrix()
{
    fill_chebyshev_points();
    for (int i=0; i < N_chebyshev; i++)
        for (int j=0; j < N_chebyshev; j++)
            chebyshev_matrix[i][j] = phi(i, chebyshev_points[j]);
}

void fill_transformation_matrix()
{
    if (transformation_matrix_initialized)
        return;
    fill_chebyshev_matrix();

    for (int i=0; i < MAX_P; i++) {
        Matrix *_mat = new DenseMatrix(N_chebyshev);
        _mat->zero();
        for (int _i=0; _i < N_chebyshev; _i++)
            for (int _j=0; _j < N_chebyshev; _j++)
                _mat->add(_i, _j, chebyshev_matrix[_i][_j]);
        double f[N_chebyshev];
        for (int j=0; j < N_chebyshev; j++)
            f[j] = lobatto_fn_tab_1d[i](chebyshev_points[j]);
        solve_linear_system(_mat, f);
        for (int j=0; j < N_chebyshev; j++)
            transformation_matrix[j][i] = f[j];
    }
    transformation_matrix_initialized = 1;
}

void transform_element_refined(int comp, double *y_prev, double *y_prev_ref, Element
        *e, Element *e_ref_left, Element *e_ref_right, Mesh *mesh, Mesh
        *mesh_ref)
{
    double y_prev_loc[MAX_P+1];
    double y_prev_loc_trans[N_chebyshev+1];
    if (e->dof[comp][0] == -1)
        y_prev_loc[0] = mesh->bc_left_dir_values[comp];
    else
        y_prev_loc[0] = y_prev[e->dof[comp][0]];
    if (e->dof[comp][1] == -1)
        y_prev_loc[1] = mesh->bc_right_dir_values[comp];
    else
        y_prev_loc[1] = y_prev[e->dof[comp][1]];
    for (int i=2; i < e->p + 1; i++)
        y_prev_loc[i] = y_prev[e->dof[comp][i]];
    fill_transformation_matrix();
    for (int i=0; i < 3 + 2*(e->p - 1); i++) {
        y_prev_loc_trans[i] = 0.;
        for (int j=0; j < e->p + 1; j++)
            y_prev_loc_trans[i] += transformation_matrix[i][j] * y_prev_loc[j];
    }
    // copying computed coefficients into the elements e_ref_left and
    // e_ref_right
    if (e->dof[comp][0] != -1)
        y_prev_ref[e_ref_left->dof[comp][0]] = y_prev_loc_trans[0];
    y_prev_ref[e_ref_left->dof[comp][1]] = y_prev_loc_trans[1];
    y_prev_ref[e_ref_right->dof[comp][0]] = y_prev_loc_trans[1];
    if (e->dof[comp][1] != -1)
        y_prev_ref[e_ref_right->dof[comp][1]] = y_prev_loc_trans[2];
    if (e_ref_left->p != e_ref_right->p)
        error("internal error in transform_element: the left and right elements must have the same order.");
    int counter = 0;
    for (int p=2; p < e_ref_left->p + 1; p++) {
        y_prev_ref[e_ref_left->dof[comp][p]] = y_prev_loc_trans[3+counter];
        counter++;
        y_prev_ref[e_ref_right->dof[comp][p]] = y_prev_loc_trans[3+counter];
        counter++;
    }
}

void transform_element_unrefined(int comp, double *y_prev, double *y_prev_ref, Element
        *e, Element *e_ref, Mesh *mesh, Mesh *mesh_ref)
{
    for (int p=0; p < e->p + 1; p++) {
        if (e->dof[comp][p] != -1)
            y_prev_ref[e_ref->dof[comp][p]] = y_prev[e->dof[comp][p]];
    }
    for (int p=e->p+1; p < e_ref->p + 1; p++) {
        y_prev_ref[e_ref->dof[comp][p]] = 0.;
    }
}

/* This only works after the dofs are assigned in the reference (and coarse)
 * solution. */
void transfer_solution(Mesh *mesh, Mesh *mesh_ref, double *y_prev, double *y_prev_ref)
{
    Iterator *I = new Iterator(mesh);
    Iterator *I_ref = new Iterator(mesh_ref);
    Element *e, *e_ref, *e_ref_left, *e_ref_right;
    for (int comp=0; comp < mesh->get_n_eq(); comp++) {
        I->reset();
        I_ref->reset();
        while ((e = I->next_active_element()) != NULL) {
            e_ref = I_ref->next_active_element();
            if (e->level == e_ref->level)
                transform_element_unrefined(comp, y_prev, y_prev_ref, e,
                    e_ref, mesh, mesh_ref);
            else if (e->level + 1 == e_ref->level) {
                e_ref_left = e_ref;
                e_ref_right = I_ref->next_active_element();
                transform_element_refined(comp, y_prev, y_prev_ref, e,
                    e_ref_left, e_ref_right, mesh, mesh_ref);
            }
            else
                error("internal error in transfer_solution: element orders mismatch.");
        }
    }
}
