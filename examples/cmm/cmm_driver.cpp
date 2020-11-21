// ==================================================================
// cmm_driver.cpp
//
// Tetrahedral element based finite element code for the 3D Coupled-Momentum
// Method using the Variational Multiscale Formulation and Generalized-
// alpha time stepping.
//
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_EBC_wall.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "CVFlowRate_Unsteady.hpp"
#include "CVFlowRate_Linear2Steady.hpp"
// #include "GenBC_Resistance.hpp"
// #include "GenBC_RCR.hpp"
// #include "GenBC_Inductance.hpp"
#include "PLocAssem_Tet_VMS_NS_GenAlpha.hpp"
// #include "PGAssem_NS_FEM.hpp"
#include "PTime_NS_Solver.hpp"

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

  double inflow_thd_time = 1.0;      // time for linearly increasing inflow to reach steady state
  double inflow_tgt_rate = 1.0;      // inflow upon reaching steady state

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
  double final_time = 1.0;           // end time of simulation
  std::string sol_bName("SOL_");     // base name of the solution file
  int    ttan_renew_freq = 1;        // frequency of tangent matrix renewal
  int    sol_record_freq = 1;        // frequency for recording the solution

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
  // Otherwise, print parameters for linear2steady inflow setting. 
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

  // ===== Read data from partition files =====
  // Control points' xyz coordinates
  FEANode * fNode = new FEANode(part_file, rank);

  // Local sub-domain's IEN array
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // Global mesh info
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);

  // Mesh partition info
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  // Local sub-domain's element indices
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // Local sub-domain's nodal (Dirichlet) BC
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);

  // Local sub-domain's inflow (Dirichlet) BC
  ALocal_Inflow_NodalBC * locinfnbc = new ALocal_Inflow_NodalBC(part_file, rank);

  // Local sub-domain's outflow elemental (Neumann) BC
  ALocal_EBC * locebc = new ALocal_EBC_outflow(part_file, rank);

  // Local sub-domain's wall elemental (Neumann) BC for CMM
  ALocal_EBC * locebc_wall = new ALocal_EBC_wall(part_file, rank, "ebc_wall");

  // Local sub-domain's nodal indices
  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::commPrint("===> Data from HDF5 files read from disk.\n");

  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) assigned for FEM analysis. \n", size);

  // ===== Prescribed inflow =====
  SYS_T::commPrint("===> Set up inflow. \n");

  ICVFlowRate * inflow_rate_ptr = nullptr;

  // If inflow file exists, prescribe it. Otherwise, prescribe an inflow that 
  // linearly increases until a steady flow rate.
  if( SYS_T::file_exist( inflow_file ) )
    inflow_rate_ptr = new CVFlowRate_Unsteady( inflow_file.c_str() );
  else
    inflow_rate_ptr = new CVFlowRate_Linear2Steady( inflow_thd_time, inflow_tgt_rate );

  inflow_rate_ptr->print_info();

  // ===== Quadrature rules =====
  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  // ===== Finite element containers =====
  SYS_T::commPrint("===> Set up volumetric and surface element containers. \n");
  FEAElement * elementv = nullptr;
  FEAElement * elements = nullptr;

  if( GMIptr->get_elemType() == 501 )          // linear tet
  {
    if( nqp_tet > 5 ) SYS_T::commPrint("Warning: Requested > 5 volumetric quadrature points for a linear tet element.\n");
    if( nqp_tri > 4 ) SYS_T::commPrint("Warning: Reqiested > 4 surface quadrature points for a linear tri element.\n");

    elementv = new FEAElement_Tet4( nqp_tet );
    elements = new FEAElement_Triangle3_3D_der0( nqp_tri );
  }
  else if( GMIptr->get_elemType() == 502 )     // quadratic tet
  {
    SYS_T::print_fatal_if( nqp_tet < 29, "Error: not enough quadrature points for quadratic tet element.\n" );
    SYS_T::print_fatal_if( nqp_tri < 13, "Error: not enough quadrature points for quadratic tri element.\n" );

    elementv = new FEAElement_Tet10_v2( nqp_tet );
    elements = new FEAElement_Triangle6_3D_der0( nqp_tri );
  }
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  // ===== Generate a sparse matrix for enforcing nodal BCs ====
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);

  pmat->gen_perm_bc(pNode, locnbc);

  // ===== Generalized-alpha =====
  SYS_T::commPrint("===> Set up the generalized-alpha time integration scheme.\n");

  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(
      genA_rho_inf, false );

  tm_galpha_ptr->print_info();

  // ===== Local assembly routine =====
  IPLocAssem * locAssem_ptr = new PLocAssem_Tet_VMS_NS_GenAlpha(
      tm_galpha_ptr, GMIptr->get_nLocBas(),
      quadv->get_num_quadPts(), elements->get_nLocBas(),
      fluid_density, fluid_mu, bs_beta, c_tauc, GMIptr->get_elemType() );


  // ===== Deallocate memory =====
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  delete locElem; delete locnbc; delete locinfnbc;
  delete locebc; delete locebc_wall; delete pNode;
  delete inflow_rate_ptr; delete quadv; delete quads;
  delete elementv; delete elements; delete pmat;
  delete tm_galpha_ptr; delete locAssem_ptr; 

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
