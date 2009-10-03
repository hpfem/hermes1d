#ifndef _WEAKFORM_H_
#define _WEAKFORM_H_

typedef double (*matrix_form) (int pts_num, double *pts, double *weights,
        double *u, double *dudx, double *v, double *dvdx, double *u_prev,
        double *du_prevdx);

typedef double (*vector_form) (int pts_num, double *pts, double *weights,
        double *u_prev, double *du_prevdx, double *v, double *dvdx);

class DiscreteProblem {


public:
    DiscreteProblem(int neq);

    void add_matrix_form(int i, int j, matrix_form fn);
    void add_vector_form(int i, int j, vector_form fn);

    int neq;
};

#endif
