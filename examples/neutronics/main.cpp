//#define NDEBUG

#include "hermes1d.h"
#include "solver_umfpack.h"

#include <map>
#include <cassert>

using std::map;
using std::cout;
using std::endl;

// ********************************************************************

// This example solves the eigenvalue problem for the neutron diffusion equation 
//	-(D.u')' + Sa.u = 1/k.nSf.u 
// in an environment composed of three slabs - inner core, outer core and a
// reflector. Reflective condition is prescribed on the left (homogeneous 
// Neumann BC) and there is vacuum on the outside of the reflector on the right
// (modelled by a Newton BC "albedo.u + D.u' = 0").

bool flag = false;				// flag for debugging purposes

// General input:
int N_elem = 60;          // number of elements
int P_init = 3;           // initial polynomal degree
int max_SI = 1000;        // Max. number of the eigenvalue iteration

// Geometry and materials
const int N_MAT = 3;			           // number of materials
const int N_GRP = 1;			           // number of energy groups in multigroup approximation
double interfaces[N_MAT+1] = { 0, 50, 100, 125 };  // Coordinates of material regions interfaces [cm]


// Matrix solver
const int MATRIX_SOLVER = 1;            // 0... default (LU decomposition)
                                        // 1... UMFPACK
                                        // 2... CG (no preconditioning)
                                        // Only relevant for iterative matrix solvers:
const double MATRIX_SOLVER_TOL = 1e-7;  // Tolerance for residual in L2 norm
const int MATRIX_SOLVER_MAXITER = 150;  // Max. number of iterations

// Newton's method
double NEWTON_TOL = 1e-5;               // tolerance for the Newton's method
int NEWTON_MAXITER = 150;               // max. number of Newton iterations
double TOL_SI = 1e-8;                   // tol. for the source (eigenvalue) iteration

// Boundary conditions
double Val_neumann_left = 0.0;		// total reflection on the left (zero Neumann)
double Val_albedo_right = 0.5; 		// vacuum on the right

// Physical properties of each material type
static double D[N_GRP][N_MAT] = 	
	{ { 0.650, 0.750, 1.150 } };		// diffusion coefficient
static double Sa[N_GRP][N_MAT] = 	
	{ { 0.120, 0.100, 0.010 } };		// absorption cross-section
static double nSf[N_GRP][N_MAT] = 
	{ { 0.185, 0.150, 0.000 } };		// fission-yield cross section (\nu \Sigma_f)
static double chi[N_GRP] = 
	{ 1.0 };				// fission spectrum (for multigroup calc.)

// Other physical properties
static double nu = 2.43; 		// mean number of neutrons released by fission
static double eps = 3.204e-11;		// mean energy release of each fission evt [J]

// Table defining the fission sources distribution in each step of the source
// iteration; in the case of multigroup approximation, one such table for each 
// group with chi^g > 0 would be needed
map<double, double> fs;

/******************************************************************************/

// Get material number at Gauss point x
int get_material_number(double x)
{
  for (int i = 0; i < N_MAT; i++) {
    if (x >= interfaces[i] && x < interfaces[i+1]) return i;
  }
  return -1;
}

// Calculate \int_e \nu\Sigma_f(x) u(x) \,\mathrm{d}x over the specified element
double calc_elem_fission_yield(Element *elem)
{
  double val_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM];
  double der_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM];
  int pts_num;
  int order = 2*elem->p;
  int m;
  elem->get_solution_quad(0,  order, val_phys, der_phys);

  double phys_x[MAX_QUAD_PTS_NUM];
  double phys_weights[MAX_QUAD_PTS_NUM];
  create_phys_element_quadrature(elem->x1, elem->x2, order, phys_x, phys_weights, &pts_num); 

  int n_grp = elem->dof_size;
  double yield, val;
  for (int i=0; i<pts_num; i++) {
    val = 0;
    m = get_material_number(phys_x[i]);
    //if (flag) cout << "\t" << phys_x[i] << "\t" << m << endl; 
   	
    for (int g=0; g<n_grp; g++) val += nSf[g][m]*val_phys[g][i];
    yield += val * phys_weights[i];
  }

  return yield;
}

// Calculate \int_\Omega \nu\Sigma_f(x) u(x) \,\mathrm{d}x over the whole mesh
double calc_fission_yield(Mesh* mesh)
{
  double fis_yield = 0;
  Iterator *I = new Iterator(mesh);

  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    //if (flag) printf("%d : (%f,%f) \n", e->id, e->x1, e->x2);
    if (flag) cout << e->id << " : (" << e->x1 << "," << e->x2 << ")" << endl;
    //if (flag) cout << "hi";
    fis_yield += calc_elem_fission_yield(e);
  }
  
  delete I;
  return fis_yield;
}

// Create table [x, f(x)], where x are the quadrature points in every element
// and f(x) are the corresponding fission source strengths in these points:
// f(x) = 1/keff \nu\Sigma_f(x) u(x)
void fis_src_distribution(Mesh* mesh, double keff, map<double, double>* fs)
{  
  // traverse 'mesh' and calculate fission source distribution over all elements
  Iterator *I = new Iterator(mesh);
  Element *e;
	
  while ((e = I->next_active_element()) != NULL) 
  {
    // quadrature order made to match the setting in 
    // DiscreteProblem::process_vol_forms(...) to produce the same number of 
    // quad points
    int order = 4*e->p; 
	
    // obtain quadrature points in element 'e'
    int pts_num;                                // num of quad points
    double phys_pts[MAX_QUAD_PTS_NUM];          // quad points
    double phys_weights[MAX_QUAD_PTS_NUM];      // quad weights

    create_phys_element_quadrature( e->x1, e->x2, order, phys_pts, phys_weights,
    															  &pts_num ); 
 
    // append contribution from this element to the fission sources table                                                	
    double u[MAX_EQN_NUM];              // neutron flux
    double dudx[MAX_EQN_NUM];           // gradient of flux
    
    for (int i = 0; i < pts_num; i++) {
      double x = phys_pts[i];
      int m = get_material_number(x);    	
      e->get_solution_point(x, u, dudx);
      (*fs)[x] = chi[0]/keff * nSf[0][m] * u[0];
    }
  }
  
  delete I;
}

// Set the solution to a constant value 'val' everywhere
void set_const_flux(Mesh* mesh, double val)
{
  Iterator *I = new Iterator(mesh);

  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    e->coeffs[0][0] = val;
    e->coeffs[0][1] = val;
  }
  
  delete I;
}

// Multiply solution at all points by 'val'
void multiply_solution(Mesh* mesh, double val)
{
  int n_dof = mesh->get_n_dof();
  double *y = new double[n_dof];
  
  Iterator *I = new Iterator(mesh);
  Element *e;
  
  while ((e = I->next_active_element()) != NULL) e->copy_coeffs_to_vector(y);

  for (int i = 0; i < n_dof; i++) y[i] *= val;
  
  I->reset();	
  while ((e = I->next_active_element()) != NULL)
    e->get_coeffs_from_vector(y);  
  	
  delete I; 
  delete [] y;
}

// Normalize the eigenfunction representing the neutron flux so that the total
// power it generates equals to 'desired_power' [W]
void normalize_to_power(Mesh* mesh, double desired_power)
{
  // Calculate total power generated by the computed flux 'u': 
  // P(u) = \eps \int_\Omega \Sigma_f(x) u(x) \,\mathrm{d}x 
  double P = eps * calc_fission_yield(mesh) / nu;
	
  // Calculate normalization constant 'c', so that P(c u) = 'desired_power' 
  double c = desired_power / P;
	
  // Multiply the computed flux by the normalization constant
  multiply_solution(mesh, c);
}

  
/******************************************************************************/
// Weak forms for Jacobi matrix and residual

// bilinear form for the Jacobi matrix 
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...basis function
// v...test function
// u_prev...previous solution
double jacobian_vol(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  double val = 0;
  int m;
	
  for(int i = 0; i<num; i++) {
	  m = get_material_number(x[i]);
	  assert(m >= 0);
    val += ( 	D[0][m]*dudx[i]*dvdx[i] + Sa[0][m]*u[i]*v[i]  ) * weights[i];
  }
  
  return val;
}

double residual_vol(int num, double *x, double *weights, 
                double u_prev[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  double val = 0;
  int m;
 	
  for(int i = 0; i<num; i++) {
	  m = get_material_number(x[i]); 
	  assert(m >= 0);
    val += (	D[0][m] * du_prevdx[0][i]*dvdx[i] + Sa[0][m] * u_prev[0][i]*v[i] 
    				- fs[x[i]]*v[i] )*weights[i];
  }

  return val;
}

double residual_surf_left(double x, double u_prev[MAX_EQN_NUM], 
        double du_prevdx[MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
  return Val_neumann_left * v; 
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


/******************************************************************************/

int main() {
  // Create coarse mesh, enumerate basis functions
  Mesh *mesh = new Mesh(interfaces[0],interfaces[N_MAT],N_elem,P_init,N_GRP);
  printf("N_dof = %d\n", mesh->assign_dofs());

	// Initial approximation: keff = 1, u = 1
	double keff = 1.0, keff_old;
  set_const_flux(mesh, 1.0);
  
  // Register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_vol);
  dp->add_vector_form(0, residual_vol);
  dp->add_vector_form_surf(0, residual_surf_left, BOUNDARY_LEFT);
  dp->add_matrix_form_surf(0, 0, jacobian_surf_right, BOUNDARY_RIGHT);
  dp->add_vector_form_surf(0, residual_surf_right, BOUNDARY_RIGHT);

  // Source iteration (power method)
  for (int i = 0; i < max_SI; i++)
  {	
    // Obtain fission source
    fis_src_distribution(mesh, keff, &fs);
		
    // Newton's loop		
    newton(dp, mesh, MATRIX_SOLVER, MATRIX_SOLVER_TOL, MATRIX_SOLVER_MAXITER,
           NEWTON_TOL, NEWTON_MAXITER);
	
    // Update the eigenvalue
    keff_old = keff;
    keff = calc_fission_yield(mesh);		
    printf("keff_%d = %f\n", i, keff);
		
    if (fabs(keff - keff_old)/keff < TOL_SI) break;
  }
	
  flag = true;
  
  // Plot the critical (i.e. steady-state) neutron flux
  Linearizer l(mesh);
  l.plot_solution("solution.gp");
  
  // Normalize so that the absolute neutron flux generates 320 Watts of energy
  // (note that, using the symmetry condition at the origin, we've solved for  
  // flux only in the right half of the reactor) 
  normalize_to_power(mesh, 320/2);	

  // Plot the solution and mesh
  l.plot_solution("solution_320W.gp");	
  mesh->plot("mesh.gp");

  printf("Done.\n");
  return 1;
}
