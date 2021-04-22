// ==================================================================
// cmm_driver.cpp
//
// Tetrahedral element based finite element code for the 3D Coupled-Momentum
// Method using the Variational Multiscale Formulation and Generalized-
// alpha time stepping.
//
// ==================================================================
#include "HDF5_Writer.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_EBC_wall.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "ALocal_Ring_NodalBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Triangle3_membrane.hpp"
#include "FEAElement_Triangle6_membrane.hpp"
#include "CVFlowRate_Unsteady.hpp"
#include "CVFlowRate_Linear2Steady.hpp"
#include "GenBC_Resistance.hpp"
#include "GenBC_RCR.hpp"
#include "GenBC_Inductance.hpp"
#include "GenBC_Coronary.hpp"
#include "PLocAssem_Tet_CMM_GenAlpha.hpp"
#include "PGAssem_Tet_CMM_GenAlpha.hpp"
#include "PTime_CMM_Solver.hpp"

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

  // Prestress solver parameters
  bool   prestress_flag = false;

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

  if(prestress_flag) SYS_T::commPrint(   "-prestress_flag: true \n");
  else SYS_T::commPrint("-prestress_flag: false \n");

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

  // ===== Quadrature rules =====
  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  // ===== Finite element containers =====
  SYS_T::commPrint("===> Set up volumetric and surface element containers. \n");
  FEAElement * elementv = nullptr;
  FEAElement * elements = nullptr;
  FEAElement * elementw = nullptr;

  if( GMIptr->get_elemType() == 501 )          // linear tet
  {
    if( nqp_tet > 5 ) SYS_T::commPrint("Warning: Requested > 5 volumetric quadrature points for a linear tet element.\n");
    if( nqp_tri > 4 ) SYS_T::commPrint("Warning: Requested > 4 surface quadrature points for a linear tri element.\n");

    elementv = new FEAElement_Tet4( nqp_tet );
    elements = new FEAElement_Triangle3_3D_der0( nqp_tri );
    elementw = new FEAElement_Triangle3_membrane( nqp_tri );
  }
  else if( GMIptr->get_elemType() == 502 )     // quadratic tet
  {
    SYS_T::print_fatal_if( nqp_tet < 29, "Error: not enough quadrature points for quadratic tet element.\n" );
    SYS_T::print_fatal_if( nqp_tri < 13, "Error: not enough quadrature points for quadratic tri element.\n" );

    elementv = new FEAElement_Tet10_v2( nqp_tet );
    elements = new FEAElement_Triangle6_3D_der0( nqp_tri );
    elementw = new FEAElement_Triangle6_membrane( nqp_tri );
  }
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  // Local sub-domain's nodal (Dirichlet) BC
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);

  // Local sub-domain's inflow (Dirichlet) BC
  ALocal_Inflow_NodalBC * locinfnbc = new ALocal_Inflow_NodalBC(part_file, rank);

  // Local sub-domain's ring (Dirichlet) in-plane motion BC
  ALocal_Ring_NodalBC * locringnbc = new ALocal_Ring_NodalBC(part_file, rank);

  // Local sub-domain's outflow elemental (Neumann) BC
  ALocal_EBC * locebc = new ALocal_EBC_outflow(part_file, rank);

  // Local sub-domain's wall elemental (Neumann) BC for CMM
  ALocal_EBC * locebc_wall = new ALocal_EBC_wall(part_file, rank, quads, "ebc_wall", prestress_flag);

  // Cross check fluid densities specified for the solver vs. wall youngsmod calculation
  if( locebc_wall -> get_fluid_density() != fluid_density )
    SYS_T::commPrint("Warning: Assigned fluid density does not match that used to compute "
                     "wall youngsmod in the preprocessor.\n");

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
    inflow_rate_ptr = new CVFlowRate_Unsteady( inflow_file.c_str(), prestress_flag );
  else
    inflow_rate_ptr = new CVFlowRate_Linear2Steady( inflow_thd_time, inflow_tgt_rate, prestress_flag );

  inflow_rate_ptr->print_info();

  // ===== Generate a sparse matrix for enforcing nodal BCs ====
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);

  pmat->gen_perm_bc(pNode, locnbc);

  // ===== Generalized-alpha =====
  SYS_T::commPrint("===> Set up the generalized-alpha time integration scheme.\n");

  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha( genA_rho_inf, false );

  tm_galpha_ptr->print_info();

  // ===== Local assembly routine =====
  IPLocAssem * locAssem_ptr = new PLocAssem_Tet_CMM_GenAlpha(
      tm_galpha_ptr, GMIptr->get_nLocBas(),
      quadv->get_num_quadPts(), elements->get_nLocBas(),
      fluid_density, fluid_mu, bs_beta,
      wall_density, wall_poisson, wall_kappa,
      c_tauc, GMIptr->get_elemType() );

  // ===== Initial condition =====
  // base generates a parabolic velocity profile at the inlet with unit flow rate
  PDNSolution * base = new PDNSolution_NS( pNode, fNode, locinfnbc, 1 );
  
  PDNSolution * sol = new PDNSolution_NS( pNode, 0 );

  PDNSolution * dot_sol = new PDNSolution_NS( pNode, 0 );

  PDNSolution * sol_wall_disp = new PDNSolution_Wall_Disp( pNode, fNode, 0 );
  
  PDNSolution * dot_sol_wall_disp = new PDNSolution_Wall_Disp( pNode, fNode, 0 );

  if( is_restart )
  {
    initial_index = restart_index;
    initial_time  = restart_time;
    initial_step  = restart_step;

    // Read in pres, velo
    SYS_T::file_check(restart_name.c_str());
    sol->ReadBinary(restart_name.c_str());

    // Read in dot pres, dot velo
    std::string restart_dot_name = "dot_" + restart_name;
    SYS_T::file_check(restart_dot_name.c_str());
    dot_sol->ReadBinary(restart_dot_name.c_str());

    // Read in wall disp: may not exist if restart from rigid soln
    std::string restart_disp_name = restart_name + "_disp";
    std::string restart_dot_disp_name = "dot_" + restart_name + "_disp";
    if( SYS_T::file_exist( restart_disp_name ) )
    {
      sol_wall_disp->ReadBinary(restart_disp_name.c_str());

      // Read in dot wall disp if wall disp exists
      SYS_T::file_check(restart_dot_disp_name.c_str());
      dot_sol_wall_disp->ReadBinary(restart_dot_disp_name.c_str());
    }

    SYS_T::commPrint("===> Read sol from disk as a restart run... \n");
    SYS_T::commPrint("     restart_name: %s \n", restart_name.c_str());
    SYS_T::commPrint("     restart_dot_name: %s \n", restart_dot_name.c_str());
    if( SYS_T::file_exist( restart_disp_name ) )
    {
      SYS_T::commPrint("     restart_disp_name: %s \n", restart_disp_name.c_str());
      SYS_T::commPrint("     restart_dot_disp_name: %s \n", restart_dot_disp_name.c_str());
    }
    SYS_T::commPrint("     restart_time: %e \n", restart_time);
    SYS_T::commPrint("     restart_index: %d \n", restart_index);
    SYS_T::commPrint("     restart_step: %e \n", restart_step);
  }

  // ===== Time step info =====
  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // ===== LPN models =====
  IGenBC * gbc = nullptr;

  if( GENBC_T::get_genbc_file_type( lpn_file.c_str() ) == 1  )
    gbc = new GenBC_Resistance( lpn_file.c_str() );
  else if( GENBC_T::get_genbc_file_type( lpn_file.c_str() ) == 2  )
    gbc = new GenBC_RCR( lpn_file.c_str(), 1000, initial_step );
  else if( GENBC_T::get_genbc_file_type( lpn_file.c_str() ) == 3  )
    gbc = new GenBC_Inductance( lpn_file.c_str() );
  else if( GENBC_T::get_genbc_file_type( lpn_file.c_str() ) == 4  )
    gbc = new GenBC_Coronary( lpn_file.c_str(), 1000, initial_step );
  else
    SYS_T::print_fatal( "Error: GenBC input file %s format cannot be recongnized.\n", lpn_file.c_str() );

  gbc -> print_info();

  // Make sure the gbc number of faces matches that of ALocal_EBC
  SYS_T::print_fatal_if(gbc->get_num_ebc() != locebc->get_num_ebc(),
      "Error: GenBC number of faces does not match with that in ALocal_EBC.\n");

  // ===== Global assembly =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");

  IPGAssem * gloAssem_ptr = new PGAssem_Tet_CMM_GenAlpha( locAssem_ptr, elements, quads,
      GMIptr, locElem, locIEN, pNode, locnbc, locringnbc, locebc, gbc, nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr,
      elements, quads, locIEN, pNode, locnbc, locringnbc, locebc, gbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // ===== Initialize the dot_sol vector by solving the mass matrix =====
  if( is_restart == false )
  {
    SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n");
    PLinear_Solver_PETSc * lsolver_acce = new PLinear_Solver_PETSc(
        1.0e-14, 1.0e-85, 1.0e30, 1000, "mass_", "mass_" );

    KSPSetType(lsolver_acce->ksp, KSPGMRES);
    KSPGMRESSetOrthogonalization(lsolver_acce->ksp,
        KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPGMRESSetRestart(lsolver_acce->ksp, 500);

    PC preproc; lsolver_acce->GetPC(&preproc);
    PCSetType( preproc, PCHYPRE );
    PCHYPRESetType( preproc, "boomeramg" );

    gloAssem_ptr->Assem_mass_residual( sol, locElem, locAssem_ptr, elementv,
        elements, quadv, quads, locIEN, pNode, fNode, locnbc, locringnbc, locebc );

    lsolver_acce->Solve( gloAssem_ptr->K, gloAssem_ptr->G, dot_sol );

    dot_sol -> ScaleValue(-1.0);

    SYS_T::commPrint("\n===> Consistent initial acceleration is obtained. \n");
    lsolver_acce -> print_info();
    delete lsolver_acce;
    SYS_T::commPrint(" The mass matrix lsolver is destroyed. \n\n");
  }

  // ===== Linear solver context =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  PC upc; lsolver->GetPC(&upc);
  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields);
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);

  // ===== Nonlinear solver context =====
  PNonlinear_CMM_Solver * nsolver = new PNonlinear_CMM_Solver( pNode, fNode,
      nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq, nl_threshold, prestress_flag );

  nsolver->print_info();

  // ===== Temporal solver context =====
  PTime_CMM_Solver * tsolver = new PTime_CMM_Solver( sol_bName,
      sol_record_freq, ttan_renew_freq, final_time, prestress_flag );

  tsolver->print_info();

  // ===== Outlet data recording files =====
  for(int ff=0; ff<locebc->get_num_ebc(); ++ff)
  {
    const double dot_face_flrate = gloAssem_ptr -> Assem_surface_flowrate(
        dot_sol, locAssem_ptr, elements, quads, locebc, ff );

    const double face_flrate = gloAssem_ptr -> Assem_surface_flowrate(
        sol, locAssem_ptr, elements, quads, locebc, ff );

    const double face_avepre = gloAssem_ptr -> Assem_surface_ave_pressure(
        sol, locAssem_ptr, elements, quads, locebc, ff );

    // set the gbc initial conditions using the 3D data
    gbc -> reset_initial_sol( ff, face_flrate, face_avepre, timeinfo->get_time() );

    const double dot_lpn_flowrate = dot_face_flrate;
    const double lpn_flowrate = face_flrate;
    const double lpn_pressure = gbc -> get_P( ff, dot_lpn_flowrate, lpn_flowrate );

    // Create the txt files and write the initial flow rates
    if(rank == 0)
    {
      std::ofstream ofile;

      // If this is NOT a restart run, generate a new file, otherwise append to
      // existing file
      if( !is_restart )
        ofile.open( locebc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::trunc );
      else
        ofile.open( locebc->gen_flowfile_name(ff).c_str(), std::ofstream::out | std::ofstream::app );

      // If this is NOT a restart, then record the initial values
      if( !is_restart )
      {
        ofile<<"Time index"<<'\t'<<"Time"<<'\t'<<"dot Flow rate"<<'\t'<<"Flow rate"<<'\t'<<"Face averaged pressure"<<'\t'<<"Reduced model pressure"<<'\n';
        ofile<<timeinfo->get_index()<<'\t'<<timeinfo->get_time()<<'\t'<<dot_face_flrate<<'\t'<<face_flrate<<'\t'<<face_avepre<<'\t'<<lpn_pressure<<'\n';
      }

      ofile.close();
    }
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== Inlet data recording files =====
  const double inlet_face_flrate = gloAssem_ptr -> Assem_surface_flowrate(
      sol, locAssem_ptr, elements, quads, locinfnbc );

  const double inlet_face_avepre = gloAssem_ptr -> Assem_surface_ave_pressure(
      sol, locAssem_ptr, elements, quads, locinfnbc );

  if( rank == 0 )
  {
    std::ofstream ofile;
    if( !is_restart )
      ofile.open( locinfnbc->gen_flowfile_name().c_str(), std::ofstream::out | std::ofstream::trunc );
    else
      ofile.open( locinfnbc->gen_flowfile_name().c_str(), std::ofstream::out | std::ofstream::app );

    if( !is_restart )
    {
      ofile<<"Time index"<<'\t'<<"Time"<<'\t'<<"Flow rate"<<'\t'<<"Face averaged pressure"<<'\n';
      ofile<<timeinfo->get_index()<<'\t'<<timeinfo->get_time()<<'\t'<<inlet_face_flrate<<'\t'<<inlet_face_avepre<<'\n';
    }

    ofile.close();
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");

  tsolver->TM_CMM_GenAlpha(is_restart, base, dot_sol, sol, dot_sol_wall_disp, sol_wall_disp,
      tm_galpha_ptr, timeinfo, inflow_rate_ptr, locElem, locIEN, pNode, fNode,
      locnbc, locinfnbc, locringnbc, locebc, locebc_wall, gbc, pmat, elementv, elements, elementw,
      quadv, quads, locAssem_ptr, gloAssem_ptr, lsolver, nsolver);

  // ===== Print complete solver info =====
  lsolver -> print_info();

  // ===== Deallocate memory =====
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic; delete locElem;
  delete locnbc; delete locinfnbc; delete locringnbc; delete locebc; delete locebc_wall;
  delete pNode; delete inflow_rate_ptr; delete quadv; delete quads; delete elementv;
  delete elements; delete elementw; delete pmat; delete tm_galpha_ptr;
  delete locAssem_ptr; delete base; delete sol; delete dot_sol;
  delete sol_wall_disp; delete dot_sol_wall_disp; delete timeinfo; delete gbc;
  delete gloAssem_ptr; delete lsolver; delete nsolver; delete tsolver;

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
