#ifndef _WEAKFORM_H_
#define _WEAKFORM_H_

typedef double (*biform_val_t) (int pts_num, double *pts, double *weights,
        double *u, double *dudx, double *v, double *dvdx, double *u_prev,
        double *du_prevdx);

class WeakForm {


public:
    WeakForm(int neq);

    void add_biform(int i, int j, biform_val_t fn);

    int neq;
};

#endif
