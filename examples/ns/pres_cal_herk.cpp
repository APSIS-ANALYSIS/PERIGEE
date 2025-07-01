// ==================================================================
// pressure calculate.cpp
//
// Finite element code to calculate the pressure for the 3D Navier-Stokes equations 
// program using the Variational Multiscale Formulation and Half-explicit RK time stepping.
//
// Author: Yujie Sun
// Date: June. 25 2025
// ==================================================================

#include "HDF5_Writer.hpp"
#include "ANL_Tools.hpp"
#include "FlowRateFactory.hpp"
#include "PGAssem_Block_NS_FEM_HERK.hpp"
#include "PTime_NS_HERK_Solver.hpp"
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

int main(int argc, char *argv[])
{
  // solution file name to be loaded for calculating the pressure
  std::string read_sol_bname("SOL_");
  int time_start = 0, time_step = 1, time_end = 1;

  // base name of the dot_solution file
  std::string sol_bname("dot_SOL_");

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  SYS_T::print_perigee_art();

  // Read analysis code parameter if the solver_cmd.h5 exists
  SYS_T::commPrint("===> Data from HDF5 files are read from disk.\n");
 
  hid_t solcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( solcmd_file );

  double initial_step = cmd_h5r -> read_doubleScalar("/","init_step");
  int nqp_vol         = cmd_h5r -> read_intScalar("/", "nqp_vol");
  int nqp_sur         = cmd_h5r -> read_intScalar("/", "nqp_sur");
  double fluid_density = cmd_h5r -> read_doubleScalar("/", "fl_density");
  double fluid_mu = cmd_h5r -> read_doubleScalar("/", "fl_mu");

  delete cmd_h5r; H5Fclose(solcmd_file);

  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5_Reader * pcmd_h5r = new HDF5_Reader( prepcmd_file );

  std::string part_file = pcmd_h5r -> read_string("/", "part_file" );
  const std::string elemType_str = pcmd_h5r -> read_string("/","elemType");
  const FEType elemType = FE_T::to_FEType(elemType_str);
  
  delete pcmd_h5r; H5Fclose(prepcmd_file);

  // Estimate of the nonzero per row for the sparse matrix
  int nz_estimate = 300;

  // Stabilization para for Darcy problem
  double L0 = 0.1;
  double cu = 2.0;
  double cp = 2.0;

  // dot_inflow_file
  std::string dot_inflow_file("dot_inflow_fourier_series.txt");

  // Yaml options
  bool is_loadYaml = true;
  std::string yaml_file("./runscript.yml");

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis. \n", size);

  SYS_T::print_fatal_if( size!= ANL_T::get_cpu_size(part_file, rank),
      "Error: Assigned CPU number does not match the partition. \n");

  // ===== Yaml Arguments =====
  SYS_T::GetOptionBool("-is_loadYaml", is_loadYaml);
  SYS_T::GetOptionString("-yaml_file", yaml_file);

  if (is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionInt("-nqp_vol", nqp_vol);
  SYS_T::GetOptionInt("-nqp_sur", nqp_sur);
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);
  SYS_T::GetOptionReal("-fl_density", fluid_density);
  SYS_T::GetOptionReal("-fl_mu", fluid_mu);
  SYS_T::GetOptionReal("-L0", L0);
  SYS_T::GetOptionReal("-cu", cu);
  SYS_T::GetOptionReal("-cp", cp);
  SYS_T::GetOptionString("-dot_inflow_file", dot_inflow_file);
  //   SYS_T::GetOptionString("-lpn_file", lpn_file);
  SYS_T::GetOptionReal("-dt", initial_step);
  SYS_T::GetOptionString("-read_sol_name", read_sol_bname);
  SYS_T::GetOptionString("-sol_name", sol_bname);
  SYS_T::GetOptionString("-part_file", part_file);

  // time stepping parameters
  // Assuming SOL_900000000 corresponds to a time of 0.0s
  double initial_time = time_start * initial_step;
  int initial_index = time_start;
  double final_time = initial_time + time_end * initial_step;

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-nqp_vol:", nqp_vol);
  SYS_T::cmdPrint("-nqp_sur:", nqp_sur);
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-fl_density:", fluid_density);
  SYS_T::cmdPrint("-fl_mu:", fluid_mu);
  SYS_T::cmdPrint("-L0:", L0);
  SYS_T::cmdPrint("-cu:", cu);
  SYS_T::cmdPrint("-cp:", cp);
  SYS_T::cmdPrint("-dot_inflow_file:", dot_inflow_file);
//   SYS_T::cmdPrint("-lpn_file:", lpn_file);
  SYS_T::cmdPrint("-read_sol_name:", read_sol_bname);
  SYS_T::cmdPrint("-sol_name:", sol_bname);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-init_time:", initial_time);
  SYS_T::cmdPrint("-dt:", initial_step);
  SYS_T::cmdPrint("-init_index:", initial_index);
  SYS_T::cmdPrint("-fina_time:", final_time);

  // ===== Record important solver options =====
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("solver_pres_cmd.h5",
        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

    cmdh5w->write_doubleScalar("fl_density", fluid_density);
    cmdh5w->write_doubleScalar("fl_mu", fluid_mu);
    cmdh5w->write_doubleScalar("init_step", initial_step);
    cmdh5w->write_intScalar("nqp_vol", nqp_vol);
    cmdh5w->write_intScalar("nqp_sur", nqp_sur);
    // cmdh5w->write_string("lpn_file", lpn_file);
    cmdh5w->write_string("dot_inflow_file", dot_inflow_file);
    cmdh5w->write_string("sol_bName", sol_bname);
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

  // dot_inflow rate
  auto dot_inflow_rate = FlowRateFactory::createFlowRate(dot_inflow_file);
  
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

  // ===== HERK Local Assembly routine =====
  auto locAssem = SYS_T::make_unique<PLocAssem_Block_VMS_NS_HERK>(
        ANL_T::get_elemType(part_file, rank), nqp_vol, nqp_sur, nullptr,
        fluid_density, fluid_mu, L0, cu, cp );

  // ===== Initial condition =====
  // dot_sol stores boundary base of dot_velocity
   std::unique_ptr<PDNSolution> base =
    SYS_T::make_unique<PDNSolution_NS>( pNode.get(), fNode.get(), locinfnbc.get(), 1 ); 
 
  // sol stores velocity and pressure (pressure equals to 0)
  std::unique_ptr<PDNSolution> sol =
    SYS_T::make_unique<PDNSolution_NS>( pNode.get(), 0 );

  // dot_sol stores dot_velocity and pressure
  std::unique_ptr<PDNSolution> dot_sol =
    SYS_T::make_unique<PDNSolution_NS>( pNode.get(), 0 );
  
  // pres sol stores pressure
  std::unique_ptr<PDNSolution> pres =
    SYS_T::make_unique<PDNSolution_P>( pNode.get(), 0, true, "pres" );
  
  // dot_velo stores dot velocity
  std::unique_ptr<PDNSolution> dot_velo =
    SYS_T::make_unique<PDNSolution_V>( pNode.get(), 0, true, "dot_velo" );

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

  gloAssem->Assem_tangent_matrix(1.0);

  // ===== Initialize the shell tangent matrix =====
  Mat K_shell;
  
  MatCreateShell( PETSC_COMM_WORLD, local_row_size, local_col_size,
    PETSC_DETERMINE, PETSC_DETERMINE, (void *)gloAssem.get(), &K_shell);

  MatShellSetOperation(K_shell, MATOP_MULT, (void(*)(void))MF_T::MF_MatMult);

  // ===== Linear solver context =====
  auto lsolver = SYS_T::make_unique<PLinear_Solver_PETSc>();

  // ===== Linear solver context of Martrix A =====
  auto lsolver_A = SYS_T::make_unique<PLinear_Solver_PETSc>(
    1.0e-8, 1.0e-15, 1.0e30, 1000, "A_", "A_");

  lsolver_A->SetOperator(gloAssem->subK[3]); 

  // ===== Linear solver context of Schur complement matrix =====
  auto lsolver_S = SYS_T::make_unique<PLinear_Solver_PETSc>(
    1.0e-8, 1.0e-15, 1.0e30, 1000, "S_", "S_");
  
  Mat S_approx;
  
  MF_T::SetupApproxSchur(gloAssem.get(), S_approx);

  lsolver_S->SetOperator(S_approx);
  
  auto solverCtx = SYS_T::make_unique<MF_T::SolverContext>(gloAssem.get(), std::move(lsolver_A), std::move(lsolver_S));

  // ===== Initialize the shell preconditioner =====
  PC pc_shell;

  PCCreate(PETSC_COMM_WORLD, &pc_shell);
  PCSetType( pc_shell, PCSHELL );
  PCShellSetContext(pc_shell, solverCtx.get());
  PCShellSetApply(pc_shell, MF_T::MF_PCSchurApply);

  KSPSetPC( lsolver->ksp, pc_shell ); 
  
  lsolver->SetOperator(K_shell); 
 
  // ===== Temporal solver context =====
  auto tsolver = SYS_T::make_unique<PTime_NS_HERK_Solver>(
      std::move(gloAssem), std::move(lsolver), std::move(pmat), nullptr,
      nullptr, std::move(dot_inflow_rate), std::move(base),
      std::move(locinfnbc), sol_bname, nlocalnode, 1, final_time );

  tsolver->print_info();

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  
  // Read the sol
  std::ostringstream time_index;
  for (int time = time_start; time<=time_end; time += time_step)
  {
    std::string name_to_read(read_sol_bname);
    std::string name_to_write(sol_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());

    SYS_T::commPrint("Time %f: Read %s, Write %s.\n", 
      time*initial_step, name_to_read.c_str(), name_to_write.c_str() );
    
    SYS_T::file_check(name_to_read);
      sol->ReadBinary(name_to_read);

    tsolver->Cal_NS_pres(sol.get(), dot_velo.get(), pres.get(), dot_sol.get(), time, initial_step);
  }

  // ===== Print complete solver info =====
  tsolver -> print_lsolver_info();
  solverCtx -> lsolver_A -> print_info();
  solverCtx -> lsolver_S -> print_info();

  MatDestroy(&S_approx);
  MatDestroy(&K_shell);
  PCDestroy(&pc_shell);
  tsolver.reset();
  solverCtx.reset();
  sol.reset(); dot_velo.reset(); 
  pres.get(); dot_sol.get();

  PetscFinalize();
  return EXIT_SUCCESS;  
}

// EOF
