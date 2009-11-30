#include "hermes1d.h"
#include "solver_umfpack.h"

// ********************************************************************
// This example shows solution of harmonic steady state on the homoegenous transmision line.
// Wave propagation is described by two linear differencial equation with complex coefficients.
// dU(x)/dx = -(R+j\omega L)I(x)
// dI(x)/dx = -(G+j\omega C)U(x)
// These equations are rewrited into four equations with real coefficients
// Ur' - R*Ir + \omega*L*Ii = 0
// Ui' - R*Ii - \omega*L*Ir = 0
// Ir' - G*Ur + \omega*C*Ui = 0
// Ii' - G*Ui - \omega*C*Ur = 0

// in an interval (0, 10) equipped with Dirichlet bdy conditions
// x1(0) = 1, x2(0) = 0, x3(0) = 0, x4(0) = 0

double L=25e-9;               // induktance [H/m]
double C=10e-12;              // capacitance [F/m]
double G=1e-9;                // conductance [S/m]
double R=1e-3;                // resistance [Ohm/m]
double l=10 ;                  // length of the line [m]
double omega=2*M_PI*3e8;        //
double Zl=60;                  // load impedance[Ohm]

//double Val_newton_alpha_U_Re=-R/Zl;
//double Val_newton_alpha_U_Im=-omega*L/Zl;
//double Val_newton_alpha_I_Re=-G*Zl;
//double Val_newton_alpha_I_Im=-omega*C*Zl;

// General input:
static int N_eq = 4;
int N_elem = 1000;          // number of elements
double A = 0, B = l;        // domain end points
int P_init = 2;             // initial polynomal degree

// Error tolerance
double TOL_NEWTON = 1e-2;   // tolerance for the Newton's method on basic mesh

// Boundary conditions
double Val_dir_left_1 = 1;  // real part of the voltage at the beginnig of the line
double Val_dir_left_2 = 0;  // imaginary part of the voltage at the beginnig of the line
double Val_dir_left_3 = 0;  // real part of the voltage at the beginnig of the line
double Val_dir_left_4 = 0;  // imaginary part of the voltage at the beginnig of the line

//At the end of the line is an indirect boundary condition U(l) = I(l)*Zl see below


// ********************************************************************

double jacobian_1_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_1_3(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -R*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_1_4(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += +omega*L*u[i]*v[i]*weights[i];
    }
    return val;
};


double jacobian_2_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_2_3(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -omega*L*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_2_4(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -R*u[i]*v[i]*weights[i];
    }
    return val;
};


double jacobian_3_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -G*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_3_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += +omega*C*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_3_3(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};



double jacobian_4_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -omega*C*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_4_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -G*weights[i];
    }
    return val;
};

double jacobian_4_4(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                    double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
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
                  double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                  double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[0][i] - R*u_prev[2][i] +L*omega*u_prev[3][i])*v[i]*weights[i];
    }
    return val;
};

double residual_2(int num, double *x, double *weights,
                  double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                  double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[1][i] - R*u_prev[3][i]-omega*L*u_prev[2][i])*v[i]*weights[i];
    }
    return val;
};

double residual_3(int num, double *x, double *weights,
                  double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                  double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[2][i] - G*u_prev[0][i] + omega*C*u_prev[1][i])*v[i]*weights[i];
    }
    return val;
};


double residual_4(int num, double *x, double *weights,
                  double u_prev[MAX_EQN_NUM][MAX_PTS_NUM],
                  double du_prevdx[MAX_EQN_NUM][MAX_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[3][i] - G*u_prev[1][i]-omega*C*u_prev[0][i])*v[i]*weights[i];
    }
    return val;
};


double jacobian_surf_right_U_Re(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_EQN_NUM],
                                double du_prevdx[MAX_EQN_NUM], void *user_data)
{     
    return 1*u*v;
}


double jacobian_surf_right_U_Im(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_EQN_NUM],
                                double du_prevdx[MAX_EQN_NUM], void *user_data)
{
    return -Zl*u*v;
}

double jacobian_surf_right_I_Re(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_EQN_NUM],
                                double du_prevdx[MAX_EQN_NUM], void *user_data)
{
    return 1*u*v;
}


double jacobian_surf_right_I_Im(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_EQN_NUM],
                                double du_prevdx[MAX_EQN_NUM], void *user_data)
{
    return -Zl*u*v;
}


/******************************************************************************/
int main() {
    // create mesh
    Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
    mesh->set_bc_left_dirichlet(0, Val_dir_left_1);
    mesh->set_bc_left_dirichlet(1, Val_dir_left_2);
    mesh->set_bc_left_dirichlet(2, Val_dir_left_3);
    mesh->set_bc_left_dirichlet(3, Val_dir_left_4);
    int N_dof = mesh->assign_dofs();
    printf("N_dof = %d\n", N_dof);

    // register weak forms
    DiscreteProblem *dp = new DiscreteProblem();
    dp->add_matrix_form(0, 0, jacobian_1_1);
    dp->add_matrix_form(0, 2, jacobian_1_3);
    dp->add_matrix_form(0, 3, jacobian_1_4);
    dp->add_matrix_form(1, 1, jacobian_2_2);
    dp->add_matrix_form(1, 2, jacobian_2_3);
    dp->add_matrix_form(1, 3, jacobian_2_4);
    dp->add_matrix_form(2, 0, jacobian_3_1);
    dp->add_matrix_form(2, 1, jacobian_3_2);
    dp->add_matrix_form(2, 2, jacobian_3_3);
    dp->add_matrix_form(3, 0, jacobian_4_1);
    dp->add_matrix_form(3, 1, jacobian_4_2);
    dp->add_matrix_form(3, 3, jacobian_4_4);
    dp->add_vector_form(0, residual_1);
    dp->add_vector_form(1, residual_2);
    dp->add_vector_form(2, residual_3);
    dp->add_vector_form(3, residual_4);
    dp->add_matrix_form_surf(0, 0, jacobian_surf_right_U_Re, BOUNDARY_RIGHT);
    dp->add_matrix_form_surf(0, 2, jacobian_surf_right_U_Im, BOUNDARY_RIGHT);
    dp->add_matrix_form_surf(1, 1, jacobian_surf_right_I_Re, BOUNDARY_RIGHT);
    dp->add_matrix_form_surf(1, 3, jacobian_surf_right_I_Im, BOUNDARY_RIGHT);

    // Allocate vectors res and y_prev
    double *res = new double[N_dof];
    double *y_prev = new double[N_dof];
    if (res == NULL || y_prev == NULL)
      error("res or y_prev could not be allocated in main().");

    // Set zero initial condition for the Newton's method
    for(int i=0; i<N_dof; i++) y_prev[i] = 0;

    // Newton's loop
    int newton_iterations = 1;
    CooMatrix *mat = NULL;
    while (1) {
        // Reset the matrix:
        if (mat != NULL) delete mat;
        mat = new CooMatrix();

        // Construct residual vector
        dp->assemble_matrix_and_vector(mesh, mat, res, y_prev);

        // Calculate L2 norm of residual vector
        double res_norm = 0;
        for(int i=0; i<N_dof; i++) res_norm += res[i]*res[i];
        res_norm = sqrt(res_norm);

        // If residual norm less than TOL_NEWTON
        // latest solution is in y_prev
        printf("Residual norm: %.15f\n", res_norm);
        if(res_norm < TOL_NEWTON) break;

        // Change sign of vector res
        for(int i=0; i<N_dof; i++) res[i]*= -1;

        // Solve the matrix system
        solve_linear_system_umfpack((CooMatrix*)mat, res);

        // Update y_prev by new solution which is in res
        for(int i=0; i<N_dof; i++) y_prev[i] += res[i];

        newton_iterations++;
    }
    printf("Finished Newton loop (%d iter).\n", newton_iterations);

    // plotting the solution
    Linearizer l(mesh);
    const char *out_filename = "solution.gp";
    l.plot_solution(out_filename, y_prev);

    printf("Done.\n");
    if (y_prev != NULL) delete[] y_prev;
    if (res != NULL) delete[] res;
    return 1;
}
