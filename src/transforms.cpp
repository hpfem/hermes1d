// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "transforms.h"

typedef double ProjMatrix[MAX_P+1][MAX_P+1];
typedef double TransMatrix[MAX_P+1][MAX_P+1];

TransMatrix trans_matrix_left;    // transforms coefficients of Lobatto shape functions 
                                  // from (-1, 1) to (-1, 0)
TransMatrix trans_matrix_right;   // transforms coefficients of Lobatto shape functions 
                                  // from (-1, 1) to (0, 1)
int trans_matrices_initialized = 0;

// transform values from (-1, 0) to (-1, 1)
#define map_left(x) (2*x+1)
// transform values from (0, 1) to (-1, 1)
#define map_right(x) (2*x-1)

double lobatto_left(int i, double x) // x \in (-1, 0)
{
  return lobatto_fn_tab_1d[i](map_left(x));
}

double lobatto_right(int i, double x) // x \in (0, 1)
{
  return lobatto_fn_tab_1d[i](map_right(x));
}

double lobatto(int i, double x) // x \in (-1, 1)
{
  return lobatto_fn_tab_1d[i](x);
}

// Fills projection matrix, i.e., the matrix of L2 products 
// of Lobatto shape functions transformed to (-1, 0). The matrix
// is the same for interval (0, 1).
void fill_proj_matrix(int n, ProjMatrix *proj_matrix)
{
  int order = 2*MAX_P; // two times max polyorder on reference element

  double phys_x[MAX_PTS_NUM];                  // quad points
  double phys_weights[MAX_PTS_NUM];            // quad weights
  int    pts_num = 0;
  create_phys_element_quadrature(-1, 0, order, phys_x, phys_weights,
                                 &pts_num); 

  // L2 product of Lobatto shape functions transformed to (-1, 0). 
  // Obviously this is the same as L2 product of Lobatto shape 
  // functions transformed to (0, 1).
  for (int i=0; i < n; i++) {
    for (int j=0; j < n; j++) {
      double result = 0;
      for (int k=0; k < pts_num; k++ ) {
        result += phys_weights[k] * lobatto_left(i, phys_x[k]) * lobatto_left(j, phys_x[k]);
      }
      (*proj_matrix)[i][j] = result;
    }
  }
}

void fill_trans_matrices(TransMatrix trans_matrix_left, 
                         TransMatrix trans_matrix_right)
{
    int order = 2*MAX_P;
    ProjMatrix proj_matrix;
    const int n = MAX_P + 1;
    fill_proj_matrix(n, &proj_matrix);

    // prepare quadrature in (-1, 0) and (0, 1)
    double phys_x_left[MAX_PTS_NUM];                     // quad points
    double phys_x_right[MAX_PTS_NUM];                    // quad points
    double phys_weights_left[MAX_PTS_NUM];               // quad weights
    double phys_weights_right[MAX_PTS_NUM];              // quad weights
    int    pts_num_left = 0;
    int    pts_num_right = 0;
    create_phys_element_quadrature(-1, 0, order, phys_x_left, phys_weights_left,
                                   &pts_num_left); 
    create_phys_element_quadrature(0, 1, order, phys_x_right, phys_weights_right,
                                   &pts_num_right); 

    // loop over shape functions on coarse element
    for (int j=0; j < MAX_P; j++) {
        // backup of projectionmatrix
        Matrix *mat_left = new DenseMatrix(n);
        Matrix *mat_right = new DenseMatrix(n);
        mat_left->zero();
        mat_right->zero();
        for (int r=0; r < n; r++) {
	    for (int s=0; s < n; s++) {
                mat_left->add(r, s, proj_matrix[r][s]);
                mat_right->add(r, s, proj_matrix[r][s]);
            }
        }
        // fill right-hand side vectors f_left and f_right for j-th 
        // Lobatto shape function on (-1, 0) and (0, 1), respectively
        double f_left[n];
        double f_right[n];
        for (int i=0; i < n; i++) {
          f_left[i] = 0;
          f_right[i] = 0;
          for (int k=0; k < pts_num_left; k++) {
            f_left[i] += phys_weights_left[k] * lobatto(j, phys_x_left[k]) *
                                                lobatto_left(i, phys_x_left[k]);
          }
          for (int k=0; k < pts_num_right; k++) {
            f_right[i] += phys_weights_right[k] * lobatto(j, phys_x_right[k]) *
                                                  lobatto_right(i, phys_x_right[k]);
	  }
        }
        // for each 'j' we get a new column in the 
        // transformation matrices
        solve_linear_system(mat_left, f_left);
        solve_linear_system(mat_right, f_right);
        for (int i=0; i < n; i++) {
            trans_matrix_left[i][j] = f_left[i];
            trans_matrix_right[i][j] = f_right[i];
        }
    }
    /* DEBUG
    for (int i=0; i < n; i++) {
        for (int j=0; j < p+1; j++) {
            printf("transf_matrix_left[%d][%d] = %g\n", i, j, transf_matrix_left[i][j]);
            printf("transf_matrix_right[%d][%d] = %g\n", i, j, transf_matrix_right[i][j]);
        }
        printf("\n");
    }
    //error("stop.");
    */
}

// Transfers solution from coarse mesh element 'e' to a pair of fine mesh elements 
// 'e_ref_left' and 'e_ref_right' (obtained via hp-refinement of 'e'). Used is the 
// global coefficient vector y_prev corresponding to the coarse mesh, and the 
// connectivity information on 'e', 'e_ref_left' and 'e_ref_right'. Result are 
// new solution coefficients on 'e_ref_left' and 'e_ref_right' which are stored 
// in the new global coefficient vector y_prev_ref.
// WARNING: For this to work, element DOF must be assigned correctly 
// in all three elements 'e', 'e_ref_left' and 'e_ref_right'!
void transform_element_refined(int comp, double *y_prev, double *y_prev_ref, Element
			       *e, Element *e_ref_left, Element *e_ref_right, 
                               double *bc_left_dir_values, double *bc_right_dir_values)
{
    //printf("ELEMENT: %d %f %f\n", e->id, e->x1, e->x2);
    double y_prev_loc[MAX_P+1];
    double y_prev_loc_trans_left[MAX_P+1];
    double y_prev_loc_trans_right[MAX_P+1];
    if (e->dof[comp][0] == -1)
        y_prev_loc[0] = bc_left_dir_values[comp];
    else
        y_prev_loc[0] = y_prev[e->dof[comp][0]];
    if (e->dof[comp][1] == -1)
        y_prev_loc[1] = bc_right_dir_values[comp];
    else
        y_prev_loc[1] = y_prev[e->dof[comp][1]];
    for (int i=2; i < e->p + 1; i++)
        y_prev_loc[i] = y_prev[e->dof[comp][i]];
    /*
    for (int i=0; i < e->p + 1; i++)
        printf("y_prev_loc[%d] = %f\n", i, y_prev_loc[i]);
        */
    if (trans_matrices_initialized == 0) {
      fill_trans_matrices(trans_matrix_left, trans_matrix_right);
      trans_matrices_initialized = 1;
    }
    // transform coefficients on the left son
    for (int i=0; i < e_ref_left->p; i++) {
        y_prev_loc_trans_left[i] = 0.;
        for (int j=0; j < e->p + 1; j++)
            y_prev_loc_trans_left[i] += trans_matrix_left[i][j] * y_prev_loc[j];
    }
    // transform coefficients on the right son
    for (int i=0; i < e_ref_right->p; i++) {
        y_prev_loc_trans_right[i] = 0.;
        for (int j=0; j < e->p + 1; j++)
            y_prev_loc_trans_right[i] += trans_matrix_right[i][j] * y_prev_loc[j];
    }

    // copying computed coefficients into the elements e_ref_left and
    // e_ref_right
    if (e->dof[comp][0] != -1)
        y_prev_ref[e_ref_left->dof[comp][0]] = y_prev_loc_trans_left[0];
    y_prev_ref[e_ref_left->dof[comp][1]] = y_prev_loc_trans_left[1];
    y_prev_ref[e_ref_right->dof[comp][0]] = y_prev_loc_trans_right[1];
    if (e->dof[comp][1] != -1)
        y_prev_ref[e_ref_right->dof[comp][1]] = y_prev_loc_trans_right[2];
    if (e_ref_left->p != e_ref_right->p)
        error("internal error in transform_element: the left and right elements must have the same order.");
    int counter = 0;
    for (int p=2; p < e_ref_left->p + 1; p++) {
        y_prev_ref[e_ref_left->dof[comp][p]] = y_prev_loc_trans_right[3+counter];
        counter++;
        y_prev_ref[e_ref_right->dof[comp][p]] = y_prev_loc_trans_right[3+counter];
        counter++;
    }
}

// Transfers solution from coarse mesh element 'e' to fine mesh element 'e_ref' 
// (obtained via p-refinement of 'e'). Used is the global coefficient vector 
// y_prev corresponding to the coarse mesh, and the connectivity information 
// on 'e' and 'e_ref'. Result are new solution coefficients on 'e_ref' which 
// are stored in the new global coefficient vector y_prev_ref.
// WARNING: For this to work, element DOF must be assigned correctly 
// in both elements 'e' and 'e_ref'!
void transform_element_unrefined(int comp, double *y_prev, double *y_prev_ref, 
        Element *e, Element *e_ref)
{
    for (int p=0; p < e->p + 1; p++) {
        if (e->dof[comp][p] != -1)
            y_prev_ref[e_ref->dof[comp][p]] = y_prev[e->dof[comp][p]];
    }
    for (int p=e->p+1; p < e_ref->p + 1; p++) {
        y_prev_ref[e_ref->dof[comp][p]] = 0.;
    }
}

// Transfers solution from the coarse mesh to the reference one. The
// solution will remain identical, but a new coefficient vector 
// y_prev_ref will be constructed/
// WARNING: For this to work, element DOF must be assigned correctly 
// in both the coarse and fine meshes!
void transfer_solution(Mesh *mesh, Mesh *mesh_ref, double *y_prev, 
                       double *y_prev_ref)
{
    Iterator *I = new Iterator(mesh);
    Iterator *I_ref = new Iterator(mesh_ref);

    // simultaneous traversal of 'mesh' and 'mesh_ref'
    Element *e, *e_ref, *e_ref_left, *e_ref_right;
    for (int comp=0; comp < mesh->get_n_eq(); comp++) {
        I->reset();
        I_ref->reset();
        while ((e = I->next_active_element()) != NULL) {
            e_ref = I_ref->next_active_element();
            if (e->level == e_ref->level)
                transform_element_unrefined(comp, y_prev, 
                                            y_prev_ref, e, e_ref);
            else {
                e_ref_left = e_ref;
                e_ref_right = I_ref->next_active_element();
                transform_element_refined(comp, y_prev, y_prev_ref, e,
					  e_ref_left, e_ref_right, 
                                          mesh->bc_left_dir_values, 
                                          mesh->bc_left_dir_values);
            }
        }
    }
}
