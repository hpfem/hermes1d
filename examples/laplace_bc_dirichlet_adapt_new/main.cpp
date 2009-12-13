#include "hermes1d.h"

// ********************************************************************

// This example uses hp-adaptivity to solve the Poisson equation 
// -u'' - f = 0 in an interval (A, B), equipped with Dirichlet 
// boundary conditions on both end points. A series of small reference 
// solutions (we call tehm fast trial refinements, FTR) is used both 
// to decide what elements will be refined, and how they will be refined. 

// General input:
static int N_eq = 1;
int N_elem = 2;                         // Number of elements
double A = -M_PI, B = M_PI;             // Domain end points
int P_init = 1;                         // Initial polynomal degree

// Stopping criteria for Newton
double TOL_NEWTON_COARSE = 1e-10;        // Coarse mesh
double TOL_NEWTON_REF = 1e-10;           // Reference mesh

// Adaptivity
const int ADAPT_TYPE = 0;               // 0... hp-adaptivity
                                        // 1... h-adaptivity
                                        // 2... p-adaptivity
const double THRESHOLD = 0.7;           // Refined will be all elements whose error
                                        // is greater than THRESHOLD*max_elem_error
const double TOL_ERR_REL = 1e-8;        // Tolerance for the relative error between 
                                        // the coarse mesh and reference solutions
const int NORM = 1;                     // To measure errors:
                                        // 1... H1 norm
                                        // 0... L2 norm

// Boundary conditions
double Val_dir_left = 0;                // Dirichlet condition left
double Val_dir_right = 0;               // Dirichlet condition right

// Function f(x)
double f(double x) {
  return sin(x);
  //return 2;
}

// Exact solution:
// When changing exact solution, do not 
// forget to update interval accordingly
const int EXACT_SOL_PROVIDED = 1;
double exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
  u[0] = sin(x);
  dudx[0] = cos(x);
  //u[0] = 1. - x*x;
  //dudx[0] = -2.*x;
}

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  printf("here!\n"); exit(0);

  // Create coarse mesh, set Dirichlet BC, enumerate 
  // basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left);
  mesh->set_bc_right_dirichlet(0, Val_dir_right);
  printf("N_dof = %d\n", mesh->assign_dofs());

  // Create discrete problem on coarse mesh
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian);
  dp->add_vector_form(0, residual);

  // Convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Convergence History", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");

  // Newton's loop on coarse mesh
  int success, iter_num;
  success = newton(dp, mesh, TOL_NEWTON_COARSE, iter_num);
  if (!success) error("Newton's method did not converge."); 
  printf("Finished initial coarse mesh Newton's iteration (%d iter).\n", 
         iter_num);

  // Main adaptivity loop
  int adapt_iterations = 1;
  double elem_errors[MAX_ELEM_NUM];      // This array decides what 
                                         // elements will be refined.
  ElemPtr2 ref_elem_pairs[MAX_ELEM_NUM]; // To store element pairs from the 
                                         // FTR solution. Decides how 
                                         // elements will be hp-refined. 

  while(1) {
    printf("============ Adaptivity step %d ============\n", adapt_iterations); 
    
    // For every element perform its fast trial refinement (FTR),
    // calculate the norm of the difference between the FTR
    // solution and the coarse mesh solution, and store this 
    // error in the elem_errors[] array.
    int n_elem = mesh->get_n_active_elem();
    for (int i=0; i < n_elem; i++) {

      // replicate coarse mesh incl. dof[] and coeffs[] arrays
      Mesh *mesh_ref_local = mesh->replicate();

      // refine one element starting with element 'i'
      mesh_ref_local->reference_refinement(i, 1);
      printf("Elem [%d]: fine mesh created (%d DOF).\n", 
             i, mesh_ref_local->assign_dofs());

      // Transfer solution from coarse mesh to the FTR mesh
      // FIXME: this should be done on one element only,
      //        in all others the solution remains the same !!! 
      transfer_solution_forward(mesh, mesh_ref_local);
      printf("Elem [%d]: coarse mesh solution copied to fine mesh.\n", i);

      // Newton's loop on the FTR mesh
      success = newton(dp, mesh_ref_local, TOL_NEWTON_REF, iter_num);
      if (!success) error("Newton's method did not converge."); 
      printf("Elem [%d]: finished fine mesh Newton's iteration (%d iter).\n", 
             i, iter_num);

      // Print FTR solutions (enumerated) 
      Linearizer *lxx = new Linearizer(mesh_ref_local);
      char out_filename[255];
      sprintf(out_filename, "solution_ref_%d.gp", i);
      lxx->plot_solution(out_filename);
      delete lxx;

      // Calculate norm of the difference between the locally refined 
      // and coarse mesh solutions.
      double err_est_array[MAX_ELEM_NUM];
      elem_errors[i] = calc_error_estimate(NORM, mesh, mesh_ref_local, 
                       err_est_array);
      printf("Elem [%d]: absolute error (est) = %g\n", i, elem_errors[i]);

      // Store pointers to the reference element pair for element 'i'
      // in the ref_elem_pairs[] array
      Iterator *I = new Iterator(mesh);
      Iterator *I_ref = new Iterator(mesh_ref_local);
      Element *e;
      while (I->next_active_element()->id <= i) {
  	  I_ref->next_active_element()->copy_into(ref_elem_pairs[e->id][0]);
          // coarse element 'e' was split in space
          if (e->level != ref_elem_pairs[e->id][0]->level) {
    	    I_ref->next_active_element()->copy_into(ref_elem_pairs[e->id][1]);
          }
      }
      delete I;
      delete I_ref;
      delete mesh_ref_local;
    }  

    // TODO: calculate the global error estimate based on the 
    //       difference between the coarse mesh and the 
    //       ref_elem_pairs[] array
    double err_est_total = calc_error_estimate(NORM, mesh, ref_elem_pairs);

    // TODO: use the ref_elem_pairs[] array
    double ref_sol_norm = calc_solution_norm(NORM, mesh, ref_elem_pairs);

    // Calculate an estimate of the global relative error
    double err_est_rel = err_est_total/ref_sol_norm;
    printf("Relative error (est) = %g %%\n", 100.*err_est_rel);

    // If exact solution available, also calculate exact error
    if (EXACT_SOL_PROVIDED) {
      // Calculate element errors wrt. exact solution
      double err_exact_total = calc_error_exact(NORM, mesh, exact_sol);
     
      // Calculate the norm of the exact solution
      // (using a fine subdivision and high-order quadrature)
      int subdivision = 500; // heuristic parameter
      int order = 20;        // heuristic parameter
      double exact_sol_norm = calc_solution_norm(NORM, exact_sol, N_eq, A, B,
                                                  subdivision, order);
      // Calculate an estimate of the global relative error
      double err_exact_rel = err_exact_total/exact_sol_norm;
      printf("Relative error (exact) = %g %%\n", 100.*err_exact_rel);
      graph.add_values(0, mesh->get_n_dof(), 100 * err_exact_rel);
    }

    // add entry to DOF convergence graph
    graph.add_values(1, mesh->get_n_dof(), 100 * err_est_rel);

     // Decide whether the relative error is sufficiently small
    if(err_est_rel*100 < TOL_ERR_REL) break;

    // debug
    //if (adapt_iterations == 1) break;

    /***** Perform the refinements *****/

    // Returns updated coarse mesh with the last solution on it. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, elem_errors,
          mesh, ref_elem_pairs);

    adapt_iterations++;
  }

  // TODO: replace mesh_ref with ref_element_pairs[]
  // Plot meshes, results, and errors
  /*
  adapt_plotting(mesh, mesh_ref,
           NORM, EXACT_SOL_PROVIDED, exact_sol);
  */

  // Save convergence graph
  graph.save("conv_dof.gp");

  printf("Done.\n");
  return 1;
}

















