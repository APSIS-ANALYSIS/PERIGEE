// ==================================================================
// ns code driver.cpp
//
// Finite element code for 3D Navier-Stokes equations using the 
// Variational Multiscale Formulation and Half-explicit RK time 
// stepping.
//
// Author: Yujie Sun
// Date: Apr. 4 2025
// ==================================================================
#include "HDF5_Writer.hpp"
#include "ANL_Tools.hpp"
#include "FlowRateFactory.hpp"
#include "PGAssem_Block_NS_FEM_HERK.hpp"
#include "PTime_NS_HERK_Solver_AccurateA.hpp"
#include "ExplicitRK_FERK1p2s.hpp"
#include "ExplicitRK_EMRK2p2s.hpp"
#include "ExplicitRK_HeunRK2p2s.hpp"
#include "ExplicitRK_RalstonRK2p2s.hpp"
#include "ExplicitRK_SSPRK3p3s.hpp"
#include "ExplicitRK_SSPRK3p4s.hpp"
#include "ExplicitRK_RalstonRK3p3s.hpp"
#include "ExplicitRK_PseudoSymplecticRK3p5q4s.hpp"
#include "ExplicitRK_38RuleRK4p4s.hpp"
#include "ExplicitRK_ClassicRK4p4s.hpp"
#include "ExplicitRK_RalstonRK4p4s.hpp"
#include "Matrix_Free_Tools.hpp"
#include "Matrix_Free_Tools_AccurateA.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature points for tets and triangles
  // Suggested values: 5 / 4 for linear, 17 / 13 for quadratic
  // Number of quadrature points for hexs and quadrangles
  // Suggested values: 8 / 4 for linear, 64 / 16 for quadratic
  int nqp_vol = 5, nqp_sur = 4;

  // Estimate of the nonzero per row for the sparse matrix
  int nz_estimate = 300;

  // fluid properties
  double fluid_density = 1.065;
  double fluid_mu = 3.5e-2;
  
  // Stabilization para for Darcy problem
  double L0 = 0.1;
  double cu = 2.0;
  double cp = 2.0;

  // inflow file & dot_inflow_file
  std::string inflow_file("inflow_fourier_series.txt");
  std::string dot_inflow_file("dot_inflow_fourier_series.txt");

  // LPN file
//   std::string lpn_file("lpn_rcr_input.txt");

  // part file location
  std::string part_file("part");

  // time stepping parameters
  double initial_time = 0.0; // time of the initial condition
  double initial_step = 0.1; // time step size
  int initial_index = 0;     // indiex of the initial condition
  double final_time = 1.0;   // final time
  std::string sol_bName("SOL_"); // base name of the solution file
  int sol_record_freq = 1;   // frequency of recording the solution

  // Restart options
  bool is_restart = false;
  int restart_index = 0;             // restart solution time index
  double restart_time = 0.0;         // restart time
  double restart_step = 1.0e-3;      // restart simulation time step size
  std::string restart_name = "SOL_"; // restart solution base name

  // Yaml options
  bool is_loadYaml = true;
  std::string yaml_file("./runscript.yml");

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::print_perigee_art();

  // ===== Yaml Arguments =====
  SYS_T::GetOptionBool("-is_loadYaml", is_loadYaml);
  SYS_T::GetOptionString("-yaml_file", yaml_file);

  if (is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionInt("-nqp_vol", nqp_vol);
  SYS_T::GetOptionInt("-nqp_sur", nqp_sur);
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);
  SYS_T::GetOptionReal("-fl_density", fluid_density);
  SYS_T::GetOptionReal("-fl_mu", fluid_mu);
  SYS_T::GetOptionReal("-L0", L0);
  SYS_T::GetOptionReal("-cu", cu);
  SYS_T::GetOptionReal("-cp", cp);
  SYS_T::GetOptionString("-inflow_file", inflow_file);
  SYS_T::GetOptionString("-dot_inflow_file", dot_inflow_file);
//   SYS_T::GetOptionString("-lpn_file", lpn_file);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionReal("-init_time", initial_time);
  SYS_T::GetOptionReal("-fina_time", final_time);
  SYS_T::GetOptionReal("-init_step", initial_step);
  SYS_T::GetOptionInt("-init_index", initial_index);
  SYS_T::GetOptionInt("-sol_rec_freq", sol_record_freq);
  SYS_T::GetOptionString("-sol_name", sol_bName);
  SYS_T::GetOptionBool("-is_restart", is_restart);
  SYS_T::GetOptionInt("-restart_index", restart_index);
  SYS_T::GetOptionReal("-restart_time", restart_time);
  SYS_T::GetOptionReal("-restart_step", restart_step);
  SYS_T::GetOptionString("-restart_name", restart_name);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint("-nqp_vol:", nqp_vol);
  SYS_T::cmdPrint("-nqp_sur:", nqp_sur);
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-fl_density:", fluid_density);
  SYS_T::cmdPrint("-fl_mu:", fluid_mu);
  SYS_T::cmdPrint("-L0", L0);
  SYS_T::cmdPrint("-cu", cu);
  SYS_T::cmdPrint("-cp", cp);
  SYS_T::cmdPrint("-inflow_file:", inflow_file);
  SYS_T::cmdPrint("-dot_inflow_file:", dot_inflow_file);
//   SYS_T::cmdPrint("-lpn_file:", lpn_file);
  SYS_T::cmdPrint("-part_file:", part_file);;
  SYS_T::cmdPrint("-init_time:", initial_time);
  SYS_T::cmdPrint("-init_step:", initial_step);
  SYS_T::cmdPrint("-init_index:", initial_index);
  SYS_T::cmdPrint("-fina_time:", final_time);
  SYS_T::cmdPrint("-sol_rec_freq:", sol_record_freq);
  SYS_T::cmdPrint("-sol_name:", sol_bName);
  if(is_restart)
  {
    SYS_T::commPrint("-is_restart: true \n");
    SYS_T::cmdPrint("-restart_index:", restart_index);
    SYS_T::cmdPrint("-restart_time:", restart_time);
    SYS_T::cmdPrint("-restart_step:", restart_step);
    SYS_T::cmdPrint("-restart_name:", restart_name);
  }
  else SYS_T::commPrint("-is_restart: false \n");

  // ===== Record important solver options =====
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("solver_cmd.h5",
        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

    cmdh5w->write_doubleScalar("fl_density", fluid_density);
    cmdh5w->write_doubleScalar("fl_mu", fluid_mu);
    cmdh5w->write_doubleScalar("init_step", initial_step);
    cmdh5w->write_intScalar("sol_record_freq", sol_record_freq);
    cmdh5w->write_intScalar("nqp_vol", nqp_vol);
    cmdh5w->write_intScalar("nqp_sur", nqp_sur);
    // cmdh5w->write_string("lpn_file", lpn_file);
    cmdh5w->write_string("inflow_file", inflow_file);
    cmdh5w->write_string("dot_inflow_file", dot_inflow_file);
    cmdh5w->write_string("sol_bName", sol_bName);
    delete cmdh5w; H5Fclose(cmd_file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== Data from Files =====
  // Control points' xyz coordinates
  auto fNode = SYS_T::make_unique<FEANode>(part_file, rank);

  // Local sub-domain's IEN array
  auto locIEN = SYS_T::make_unique<ALocal_IEN>(part_file, rank);
  
  // Local sub-domain's element indices
  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_file, rank);

  // Local sub-domain's nodal bc
  auto locnbc = SYS_T::make_unique<ALocal_NBC>(part_file, rank);

  // Local sub-domain's inflow bc
  auto locinfnbc = SYS_T::make_unique<ALocal_InflowBC>(part_file, rank);

  // Local sub-domain's elemental bc
  auto locebc = SYS_T::make_unique<ALocal_EBC>(part_file, rank);

  // Local sub-domain's nodal indices
  auto pNode = SYS_T::make_unique<APart_Node>(part_file, rank);

  const int nlocalnode = pNode->get_nlocalnode();

  PetscInt local_row_size, local_col_size;
  local_row_size = 4 * nlocalnode;
  local_col_size = local_row_size;

  SYS_T::commPrint("===> Data from HDF5 files are read from disk.\n");

  SYS_T::print_fatal_if( size!= ANL_T::get_cpu_size(part_file, rank),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis. \n", size);

  // ===== Inflow flow rate =====
  SYS_T::commPrint("===> Setup inflow flow rate. \n");

  auto inflow_rate = FlowRateFactory::createFlowRate(inflow_file);

  auto dot_inflow_rate = FlowRateFactory::createFlowRate(dot_inflow_file);

  inflow_rate->print_info();

  dot_inflow_rate->print_info();

  // ===== LPN models =====
//   auto gbc = GenBCFactory::createGenBC(lpn_file, initial_time, initial_step, 
//       initial_index, 1000);

//   gbc -> print_info();

  // Make sure the gbc number of faces matches that of ALocal_EBC
//   SYS_T::print_fatal_if(gbc->get_num_ebc() != locebc->get_num_ebc(),
//       "Error: GenBC number of faces does not match with that in ALocal_EBC.\n");

  // ===== Generate a sparse matrix for the enforcement of essential BCs
  auto pmat = SYS_T::make_unique<Matrix_PETSc>(pNode.get(), locnbc.get());

  pmat->gen_perm_bc(pNode.get(), locnbc.get());

  // ===== Half Explicit Runge Kutta scheme =====
  SYS_T::commPrint("===> Setup the Runge Kutta time scheme.\n");

  std::unique_ptr<ITimeMethod_RungeKutta> tm_RK = SYS_T::make_unique<ExplicitRK_PseudoSymplecticRK3p5q4s>();

  tm_RK->print_coefficients();
 
    // ===== HERK Local Assembly routine =====
  auto locAssem = SYS_T::make_unique<PLocAssem_Block_VMS_NS_HERK>(
        ANL_T::get_elemType(part_file, rank), nqp_vol, nqp_sur, tm_RK.get(),
        fluid_density, fluid_mu, L0, cu, cp );

  // ===== Initial condition =====
  std::unique_ptr<PDNSolution> base =
    SYS_T::make_unique<PDNSolution_NS>( pNode.get(), fNode.get(), locinfnbc.get(), 1 );    

  std::unique_ptr<PDNSolution> sol =
    SYS_T::make_unique<PDNSolution_NS>( pNode.get(), 0 );

  std::unique_ptr<PDNSolution> velo =
    SYS_T::make_unique<PDNSolution_V>( pNode.get(), 0, true, "velo" );

  std::unique_ptr<PDNSolution> pres =
    SYS_T::make_unique<PDNSolution_P>( pNode.get(), 0, true, "pres" );

  std::unique_ptr<PDNSolution> dot_velo =
    SYS_T::make_unique<PDNSolution_V>( pNode.get(), 0, true, "dot_velo" );

  if( is_restart )
  {
    initial_index = restart_index;
    initial_time  = restart_time;
    initial_step  = restart_step;

    // Read sol file
    SYS_T::file_check(restart_name);
    sol->ReadBinary(restart_name);

    SYS_T::commPrint("===> Read sol from disk as a restart run... \n");
    SYS_T::commPrint("     restart_name: %s \n", restart_name.c_str());
    SYS_T::commPrint("     restart_time: %e \n", restart_time);
    SYS_T::commPrint("     restart_index: %d \n", restart_index);
    SYS_T::commPrint("     restart_step: %e \n", restart_step);
  }

  // ===== Global assembly =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  auto gloAssem = SYS_T::make_unique<PGAssem_Block_NS_FEM_HERK>( 
      std::move(locIEN), std::move(locElem), std::move(fNode), 
      std::move(pNode), std::move(locnbc), std::move(locebc), 
      std::move(locAssem), nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem->Assem_nonzero_estimate();

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem->Fix_nonzero_err_str();
  gloAssem->Clear_subKG();

  // gloAssem->Assem_tangent_matrix(tm_RK.get(), initial_step);
  gloAssem->Assem_tangent_matrix(initial_step);
  
  // ===== Initialize the shell tangent matrix =====
  Mat K_shell;
  
  MatCreateShell( PETSC_COMM_WORLD, local_row_size, local_col_size,
    PETSC_DETERMINE, PETSC_DETERMINE, (void *)gloAssem.get(), &K_shell);

  MatShellSetOperation(K_shell, MATOP_MULT, (void(*)(void))MF_TA::MF_MatMult);

  // ===== Linear solver context =====
  auto lsolver = SYS_T::make_unique<PLinear_Solver_PETSc>();

  // ===== Linear solver context of Martrix A =====
  auto lsolver_A = SYS_T::make_unique<PLinear_Solver_PETSc>(
    1.0e-8, 1.0e-15, 1.0e30, 1000, "A_", "A_");

  // ===== Linear solver context of Schur complement matrix =====
  auto lsolver_S = SYS_T::make_unique<PLinear_Solver_PETSc>(
    1.0e-8, 1.0e-15, 1.0e30, 1000, "S_", "S_");
  
  auto solverCtx = SYS_T::make_unique<MF_TA::SolverContext>(std::move(gloAssem), std::move(lsolver_A), std::move(lsolver_S));

  // ===== Initialize the shell preconditioner =====
  PC pc_shell;

  PCCreate(PETSC_COMM_WORLD, &pc_shell);
  PCSetType( pc_shell, PCSHELL );
  PCShellSetContext(pc_shell, solverCtx.get());
  PCShellSetApply(pc_shell, MF_TA::MF_PCSchurApply);

  KSPSetPC( lsolver->ksp, pc_shell );   
  lsolver->SetOperator(K_shell); 

  // ===== Time step info ===== 
  auto timeinfo = SYS_T::make_unique<PDNTimeStep>(initial_index, initial_time, 
      initial_step);
 
  // ===== Temporal solver context =====
  auto tsolver = SYS_T::make_unique<PTime_NS_HERK_Solver_AccurateA>(
      std::move(solverCtx), std::move(lsolver), std::move(pmat), std::move(tm_RK),
      std::move(inflow_rate), std::move(dot_inflow_rate), std::move(base),
      std::move(locinfnbc), sol_bName, nlocalnode, sol_record_freq, final_time );

  tsolver->print_info();

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver->TM_NS_HERK(is_restart, std::move(sol), std::move(velo), std::move(dot_velo), 
      std::move(pres), std::move(timeinfo));

  // ===== Print complete solver info =====
  tsolver -> print_lsolver_info();

  MatDestroy(&K_shell);
  PCDestroy(&pc_shell);
  tsolver.reset();

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
