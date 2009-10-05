#ifndef _WEAKFORM_H_
#define _WEAKFORM_H_

#include "mesh.h"
#include "quad_std.h"
#include "lobatto.h"
#include "matrix.h"

typedef double (*matrix_form) (int pts_num, double *pts, double *weights,
        double *u, double *dudx, double *v, double *dvdx, double *u_prev,
        double *du_prevdx, void *user_data);

typedef double (*vector_form) (int pts_num, double *pts, double *weights,
        double *u_prev, double *du_prevdx, double *v, double *dvdx,
        void *user_data);

class DiscreteProblem {

public:
    DiscreteProblem(int neq, Mesh *mesh);

    void add_matrix_form(int i, int j, matrix_form fn);
    void add_vector_form(int i, vector_form fn);
    void assemble(Matrix *mat, double *res, double *y_prev, int matrix_flag);
    void assemble_matrix_and_vector(Matrix *mat, double *res, double *y_prev); 
    void assemble_matrix(Matrix *mat, double *y_prev);
    void assemble_vector(double *res, double *y_prev);

private:
    int neq;
    Mesh *mesh;
    matrix_form _matrix_form;
    vector_form _vector_form;

};

void element_quadrature(double a, double b, 
                        int order, double *pts, double *weights, int &num);
void element_shapefn(double a, double b, 
		     int k, int order, double *val, double *der);

#endif
