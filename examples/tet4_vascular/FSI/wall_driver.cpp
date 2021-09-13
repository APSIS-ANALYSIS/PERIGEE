// ==================================================================
// wall_driver.cpp
//
// wall mechanics solver.
//
// Date: Sep. 13 2021
// ==================================================================

int main( int argc, char *argv[] )
{
  // Clean potentially pre-existing hdf5 files of prestress
  SYS_T::execute("rm -rf prestress_p*.h5");

  // (Pseudo-) time integration parameters
  double genA_rho_inf = 0.0;
  bool is_backward_Euler = true;

  // Estimate of num nonzeros per row for the sparse tangent matrix
  int nz_estimate = 300;

  // Prestress tolerance
  double prestress_disp_tol = 1.0e-6;

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-3;           // convergence criterion relative tolerance
  double nl_atol = 1.0e-6;           // convergence criterion absolute tolerance
  double nl_dtol = 1.0e3;            // divergence criterion
  int    nl_maxits = 20;             // maximum number if nonlinear iterations
  int    nl_refreq = 4;              // frequency of tangent matrix renewal
  int    nl_threshold = 4;           // threshold of tangent matrix renewal

  // Time stepping parameters
  double initial_time = 0.0;         // time of initial condition
  double initial_step = 0.1;         // time step size
  int    initial_index = 0;          // index of initial condition
  double final_time = 1.0;           // end time of simulation
  bool   is_record_sol = false;      // bool flag to decide if one wants to record the solution
  std::string sol_bName("PS_");      // base name of the solution file
  int    ttan_renew_freq = 1;        // frequency of tangent matrix renewal
  int    sol_record_freq = 1;        // frequency for recording the solution



  return EXIT_SUCCESS;
}

// EOF
