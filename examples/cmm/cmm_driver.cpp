// ==================================================================
// cmm_driver.cpp
//
// Tetrahedral element based finite element code for 3D Coupled-Momentum
// Method using the Variational Multiscale Formulation and Generalize-
// alpha time stepping.
//
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"

int main( int argc, char *argv[] )
{
  // Number of quadrature points for tets and triangles
  // Suggested values: 5 / 4 for linear, 29 / 13 for quadratic
  int nqp_tet = 5, nqp_tri = 4;

  // Estimate of num nonzeros per row for the sparse tangent matrix 
  int nz_estimate = 300;

  // Fluid properties
  double fluid_density = 1.065;
  double fluid_mu = 3.5e-2;
  double c_tauc = 1.0;               // scaling factor for tau_c: 0.0, 0.125, or 1.0

  // Wall properties: density, Poisson ratio, shear correction factor (kappa)
  double wall_density = 1.0;
  double wall_poisson = 0.5;
  double wall_kappa   = 5.0 / 6.0;

  // Inflow file
  std::string inflow_file("inflow_fourier_series.txt");

  double inflow_thd_time = 1.0;      // prescribed time for inflow to reach steadness
  double inflow_tgt_rate = 1.0;      // prescribed flow rate at steady state

  // LPN file
  std::string lpn_file("lpn_rcr_input.txt");

  // Backflow stabilization
  double bs_beta = 0.2;

  // Generalized-alpha rho_inf
  double genA_rho_inf = 0.5;

  // Partition filename prefix
  std::string part_file("part");

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-3;           // convergence criterion relative tolerance
  double nl_atol = 1.0e-6;           // convergence criterion absolute tolerance
  double nl_dtol = 10.0;             // divergence criterion
  int    nl_maxits = 20;             // maximum number if nonlinear iterations
  int    nl_refreq = 4;              // frequency of tangent matrix renewal
  int    nl_threshold = 4;           // threshold of tangent matrix renewal

  // Time stepping parameters
  double initial_time = 0.0;         // time of initial condition
  double initial_step = 0.1;         // time step size
  int    initial_index = 0;          // index of initial condition
  double final_time = 1.0;           // final time
  std::string sol_bName("SOL_");     // base name of the solution file
  int    ttan_renew_freq = 1;        // frequency of tangent matrix renewal
  int    sol_record_freq = 1;        // frequency of recording the solution

  // Restart options
  bool   is_restart = false;
  int    restart_index = 0;          // restart solution time index
  double restart_time = 0.0;         // restart time
  double restart_step = 1.0e-3;      // restart simulation time step size
  std::string restart_name = "SOL_"; // restart solution base name

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::print_perigee_art(); 

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading command line arguments... \n");

  SYS_T::GetOptionInt(   "-nqp_tet",         nqp_tet);
  SYS_T::GetOptionInt(   "-nqp_tri",         nqp_tri);
  SYS_T::GetOptionInt(   "-nz_estimate",     nz_estimate);
  SYS_T::GetOptionReal(  "-bs_beta",         bs_beta);
  SYS_T::GetOptionReal(  "-rho_inf",         genA_rho_inf);
  SYS_T::GetOptionReal(  "-fl_density",      fluid_density);
  SYS_T::GetOptionReal(  "-fl_mu",           fluid_mu);
  SYS_T::GetOptionReal(  "-c_tauc",          c_tauc);
  SYS_T::GetOptionReal(  "-wall_density",    wall_density);
  SYS_T::GetOptionReal(  "-wall_poisson",    wall_poisson);
  SYS_T::GetOptionReal(  "-wall_kappa",      wall_kappa);
  SYS_T::GetOptionString("-inflow_file",     inflow_file);
  SYS_T::GetOptionReal(  "-inflow_thd_time", inflow_thd_time);
  SYS_T::GetOptionReal(  "-inflow_tgt_rate", inflow_tgt_rate);
  SYS_T::GetOptionString("-lpn_file",        lpn_file);
  SYS_T::GetOptionString("-part_file",       part_file);
  SYS_T::GetOptionReal(  "-nl_rtol",         nl_rtol);
  SYS_T::GetOptionReal(  "-nl_atol",         nl_atol);
  SYS_T::GetOptionReal(  "-nl_dtol",         nl_dtol);
  SYS_T::GetOptionInt(   "-nl_maxits",       nl_maxits);
  SYS_T::GetOptionInt(   "-nl_refreq",       nl_refreq);
  SYS_T::GetOptionInt(   "-nl_threshold",    nl_threshold);
  SYS_T::GetOptionReal(  "-init_time",       initial_time);
  SYS_T::GetOptionReal(  "-fina_time",       final_time);
  SYS_T::GetOptionReal(  "-init_step",       initial_step);
  SYS_T::GetOptionInt(   "-init_index",      initial_index);
  SYS_T::GetOptionInt(   "-ttan_freq",       ttan_renew_freq);
  SYS_T::GetOptionInt(   "-sol_rec_freq",    sol_record_freq);
  SYS_T::GetOptionString("-sol_name",        sol_bName);
  SYS_T::GetOptionBool(  "-is_restart",      is_restart);
  SYS_T::GetOptionInt(   "-restart_index",   restart_index);
  SYS_T::GetOptionReal(  "-restart_time",    restart_time);
  SYS_T::GetOptionReal(  "-restart_step",    restart_step);
  SYS_T::GetOptionString("-restart_name",    restart_name);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint(      "-nqp_tet:",         nqp_tet);
  SYS_T::cmdPrint(      "-nqp_tri:",         nqp_tri);
  SYS_T::cmdPrint(      "-nz_estimate:",     nz_estimate);
  SYS_T::cmdPrint(      "-bs_beta:",         bs_beta);
  SYS_T::cmdPrint(      "-rho_inf:",         genA_rho_inf);
  SYS_T::cmdPrint(      "-fl_density:",      fluid_density);
  SYS_T::cmdPrint(      "-fl_mu:",           fluid_mu);
  SYS_T::cmdPrint(      "-c_tauc:",          c_tauc);
  SYS_T::cmdPrint(      "-wall_density:",    wall_density);
  SYS_T::cmdPrint(      "-wall_poisson:",    wall_poisson);
  SYS_T::cmdPrint(      "-wall_kappa:",      wall_kappa);

  // If the inflow file exists, print its filename.
  // Otherwise, print parameters for linear2steady inflow setting 
  if( SYS_T::file_exist( inflow_file ) )
    SYS_T::cmdPrint(    "-inflow_file:",     inflow_file);
  else
  {
    SYS_T::cmdPrint(    "-inflow_thd_time:", inflow_thd_time);
    SYS_T::cmdPrint(    "-inflow_tgt_rate:", inflow_tgt_rate);
  }

  SYS_T::cmdPrint(      "-lpn_file:",        lpn_file);
  SYS_T::cmdPrint(      "-part_file:",       part_file);
  SYS_T::cmdPrint(      "-nl_rtol:",         nl_rtol);
  SYS_T::cmdPrint(      "-nl_atol:",         nl_atol);
  SYS_T::cmdPrint(      "-nl_dtol:",         nl_dtol);
  SYS_T::cmdPrint(      "-nl_maxits:",       nl_maxits);
  SYS_T::cmdPrint(      "-nl_refreq:",       nl_refreq);
  SYS_T::cmdPrint(      "-nl_threshold:",    nl_threshold);
  SYS_T::cmdPrint(      "-init_time:",       initial_time);
  SYS_T::cmdPrint(      "-init_step:",       initial_step);
  SYS_T::cmdPrint(      "-init_index:",      initial_index);
  SYS_T::cmdPrint(      "-fina_time:",       final_time);
  SYS_T::cmdPrint(      "-ttan_freq:",       ttan_renew_freq);
  SYS_T::cmdPrint(      "-sol_rec_freq:",    sol_record_freq);
  SYS_T::cmdPrint(      "-sol_name:",        sol_bName);

  if(is_restart)
  {
    SYS_T::commPrint(   "-is_restart: true \n");
    SYS_T::cmdPrint(    "-restart_index:",   restart_index);
    SYS_T::cmdPrint(    "-restart_time:",    restart_time);
    SYS_T::cmdPrint(    "-restart_step:",    restart_step);
    SYS_T::cmdPrint(    "-restart_name:",    restart_name);
  }
  else SYS_T::commPrint("-is_restart: false \n");

  // ===== Record important solver options =====
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("solver_cmd.h5",
        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

    cmdh5w->write_doubleScalar(  "fl_density",      fluid_density);
    cmdh5w->write_doubleScalar(  "fl_mu",           fluid_mu);
    cmdh5w->write_doubleScalar(  "wall_density",    wall_density);
    cmdh5w->write_doubleScalar(  "wall_poisson",    wall_poisson);
    cmdh5w->write_doubleScalar(  "wall_kappa",      wall_kappa);
    cmdh5w->write_doubleScalar(  "init_step",       initial_step);
    cmdh5w->write_intScalar(     "sol_record_freq", sol_record_freq);
    cmdh5w->write_string(        "lpn_file",        lpn_file);

    if( SYS_T::file_exist( inflow_file ) )
      cmdh5w->write_string(      "inflow_file",     inflow_file);
    else
    {
      cmdh5w->write_doubleScalar("inflow_thd_time", inflow_thd_time );
      cmdh5w->write_doubleScalar("inflow_tgt_rate", inflow_tgt_rate );
    }
    delete cmdh5w; H5Fclose(cmd_file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);



  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
