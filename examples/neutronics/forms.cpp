// bilinear form for the Jacobi matrix 
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...basis function
// v...test function
// u_prev...previous solution
double jacobian_vol_inner(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  int m = 0;
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[0][0]*dudx[i]*dvdx[i] + Sa[0][m]*u[i]*v[i]) * weights[i]; // inner core
  }
  return val;
}
double jacobian_vol_outer(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  int m = 1;
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[0][m]*dudx[i]*dvdx[i] + Sa[0][m]*u[i]*v[i]) * weights[i]; // outer core
  }
  return val;
}
double jacobian_vol_vacuum(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  int m = 2;
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[0][m]*dudx[i]*dvdx[i] + Sa[0][m]*u[i]*v[i]) * weights[i]; // vacuum
  }
  return val;
}

double residual_vol_inner(int num, double *x, double *weights, 
                double u_prev[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  int m = 0;
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[0][m] * du_prevdx[0][i] * dvdx[i] + Sa[0][m] * u_prev[0][i] * v[i] // inner core
    				- fs[x[i]]*v[i]) * weights[i];
  }
  return val;
}
double residual_vol_outer(int num, double *x, double *weights, 
                double u_prev[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  int m = 1;
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[0][m] * du_prevdx[0][i] * dvdx[i] + Sa[0][m] * u_prev[0][i] * v[i] // outer core
    				- fs[x[i]]*v[i]) * weights[i];
  }
  return val;
}
double residual_vol_vacuum(int num, double *x, double *weights, 
                double u_prev[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  int m = 2;
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[0][m] * du_prevdx[0][i] * dvdx[i] + Sa[0][m] * u_prev[0][i] * v[i] // vacuum
    				- fs[x[i]]*v[i]) * weights[i];
  }
  return val;
}

double residual_surf_left(double x, double u_prev[MAX_EQN_NUM], 
        double du_prevdx[MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
  return D[0][0] * Val_neumann_left * v; 
}

double jacobian_surf_right(double x, double u, double dudx,
        double v, double dvdx, double u_prev[MAX_EQN_NUM], 
        double du_prevdx[MAX_EQN_NUM], void *user_data)
{
  return Val_albedo_right*u*v;
}

double residual_surf_right(double x, double u_prev[MAX_EQN_NUM], 
        double du_prevdx[MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
  return Val_albedo_right * u_prev[0] * v; 
}
