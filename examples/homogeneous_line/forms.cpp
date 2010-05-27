double jacobian_1_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_1_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (R + L/dt)*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_2_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (G + L/dt)*u[i]*v[i]*weights[i];
    }
    return val;
};


double jacobian_2_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double residual_1(int num, double *x, double *weights,
                  double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[0][0][i] + R*u_prev[0][1][i] + L*(u_prev[0][1][i] - u_prev[1][1][i])/dt)*v[i]*weights[i];
    }
    return val;
};

double residual_2(int num, double *x, double *weights,
                  double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[0][1][i] + G*u_prev[0][0][i] + C*(u_prev[0][0][i] - u_prev[1][0][i])/dt)*v[i]*weights[i];
    }
    return val;
};


double jacobian_surf_left_U(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM],
                                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{     
    double U = sin(2*M_PI*1e6*t);
    U = 10.0;
    // cout << "time = " << t << ", U = " << U << endl;
    return U;
}

double jacobian_surf_right_I(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM],
                                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
    return 0;
}
