// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "transforms.h"

#define N_chebyshev_max (3+2*(MAX_P-1))

typedef double ChebyshevMatrix[N_chebyshev_max][N_chebyshev_max];
typedef double TransformationMatrix[N_chebyshev_max][MAX_P+1];
// for one particular p:
TransformationMatrix transformation_matrix;
// for all p:
TransformationMatrix transformation_matrix_list[MAX_P+1];
int transformation_matrix_initialized=0;

TransformationMatrix *get_transformation_matrix(int p)
{
    return & (transformation_matrix_list[p]);
}

void fill_chebyshev_points(int n, double *chebyshev_points)
{
    for (int i=0; i < n; i++)
        chebyshev_points[i] = -cos(i*M_PI/(n-1));
    /*
    for (int i=0; i < N_chebyshev_max; i++)
        printf("%f ", chebyshev_points[i]);
    printf("\n");
    printf("done.\n");
    */
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
            return lobatto_fn_tab_1d[0](map_right(x));
        else if (i % 2 == 1)
            return 0;
        else {
            return lobatto_fn_tab_1d[i/2](map_right(x));
        }
    } 
}

void fill_chebyshev_matrix(int n, ChebyshevMatrix *chebyshev_matrix)
{
    double chebyshev_points[N_chebyshev_max];
    fill_chebyshev_points(n, chebyshev_points);
    for (int i=0; i < n; i++) {
        for (int j=0; j < n; j++) {
            //printf("XXX %d %d %f \n", i, j, phi(j, chebyshev_points[i]));
            //printf("XXX %d %d %f \n", i, j, (*chebyshev_matrix)[i][j]);
            (*chebyshev_matrix)[i][j] = phi(j, chebyshev_points[i]);
            //printf("%f ", (*chebyshev_matrix)[i][j]);
        }
        //printf("\n");
    }
    //error("stop.");
}

void fill_transformation_matrix(int p, int p_ref, TransformationMatrix
        transformation_matrix)
{
    double chebyshev_points[N_chebyshev_max];
    ChebyshevMatrix chebyshev_matrix;
    int n = 3+2*(p_ref-1);
    fill_chebyshev_points(n, chebyshev_points);
    fill_chebyshev_matrix(n, &chebyshev_matrix);

    for (int j=0; j < p+1; j++) {
        // FIXME: don't compute the matrix all the time, e.g.  move this out of
        // the cycle and make sure the gaussian elimination doesn't overwrite
        // the _mat.
        Matrix *_mat = new DenseMatrix(n);
        _mat->zero();
        for (int _i=0; _i < n; _i++)
            for (int _j=0; _j < n; _j++)
                _mat->add(_i, _j, chebyshev_matrix[_i][_j]);
        printf("chebyshev stuff:\n");
        for (int _i=0; _i < n; _i++) {
            for (int _j=0; _j < n; _j++) {
                printf("%f ", chebyshev_matrix[_i][_j]);
            }
            printf("\n");
        }
        printf("----- END ---- \n");
        double f[n];
        for (int i=0; i < n; i++)
            f[i] = lobatto_fn_tab_1d[j](chebyshev_points[i]);
        printf("chebyshev_points\n");
        for (int i=0; i < 5; i++)
            printf("%f ", chebyshev_points[i]);
        printf("\n");
        printf("XXXXXX\n");
        for (int i=0; i < n; i++)
            printf("%f ", f[i]);
        printf("\n");
        solve_linear_system(_mat, f);
        for (int i=0; i < n; i++)
            transformation_matrix[i][j] = f[i];
    }
    for (int i=0; i < n; i++) {
        for (int j=0; j < p+1; j++) {
            printf("%f ", transformation_matrix[i][j]);
        }
        printf("\n");
    }
    error("stop.");
}

void transform_element_refined(int comp, double *y_prev, double *y_prev_ref, Element
        *e, Element *e_ref_left, Element *e_ref_right, Mesh *mesh, Mesh
        *mesh_ref)
{
    printf("ELEMENT: %d %f %f\n", e->id, e->x1, e->x2);
    double y_prev_loc[MAX_P+1];
    double y_prev_loc_trans[N_chebyshev_max+1];
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
    for (int i=0; i < e->p + 1; i++)
        printf("y_prev_loc[%d] = %f\n", i, y_prev_loc[i]);
    TransformationMatrix transformation_matrix;
    fill_transformation_matrix(e->p, e_ref_left->p, transformation_matrix);
    //fill_transformation_matrix();
    //double TransformationMatrix *transformation_matrix =
    //    get_transformation_matrix(e_ref_left->p + 1);
    for (int i=0; i < 3 + 2*(e_ref_left->p - 1); i++) {
        y_prev_loc_trans[i] = 0.;
        for (int j=0; j < e->p + 1; j++)
            y_prev_loc_trans[i] += transformation_matrix[i][j] * y_prev_loc[j];
    }
    for (int i=0; i < 3 + 2*(e_ref_left->p - 1); i++)
        printf("y_prev_loc_trans[%d] = %f\n", i, y_prev_loc_trans[i]);
    printf("----------------------\n");
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
