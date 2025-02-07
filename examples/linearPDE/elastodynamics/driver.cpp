// ============================================================================
// driver_elastodynamics.cpp
// 
// Finite element code for 3D elastodynamics equations using non-overshoot
// generalized alpha time stepping.
//
// Date: Oct. 24 2023
// ============================================================================
#include "HDF5_Writer.hpp"
#include "ANL_Tools.hpp"
#include "PLocAssem_Elastodynamics_GenAlpha.hpp"
#include "PGAssem_LinearPDE_GenAlpha.hpp"
#include "PTime_LinearPDE_Solver.hpp"

int main(int argc, char *argv[])
{
  // material parameters
  double rho = 1.0, module_E = 1.0e+7, nu = 0.3;

  // number of quadrature points
  int nqp_vol = 5, nqp_sur = 4;

  // generalized-alpha rho_inf
  double genA_rho_inf = 0.5;
  bool is_backward_Euler = false;

  // Estimate of the nonzero per row for the sparse matrix
  int nz_estimate = 300;

  // part file location
  std::string part_file("./apart/part");
  
  // nonlinear solver parameters
  double nl_rtol = 1.0e-3; // convergence criterion relative tolerance
  double nl_atol = 1.0e-6; // convergence criterion absolute tolerance
  double nl_dtol = 10.0;   // divergence criterion
  int nl_maxits = 20;      // maximum number if nonlinear iterations
  int nl_refreq = 4;       // frequency of tangent matrix renewal
  int nl_threshold = 4;    // threshold of tangent matrix renewal

  // time stepping parameters
  double initial_time = 0.0; // time of the initial condition
  double initial_step = 0.1; // time step size
  int initial_index = 0;     // indiex of the initial condition
  double final_time = 1.0;   // final time
  std::string sol_bName("SOL_"); // base name of the solution file
  int ttan_renew_freq = 1;   // frequency of tangent matrix renewal
  int sol_record_freq = 1;   // frequency of recording the solution

  // Restart options
  bool is_restart = false;
  int restart_index = 0;     // restart solution time index
  double restart_time = 0.0; // restart time
  double restart_step = 1.0e-3; // restart simulation time step size
  std::string restart_u_name = "SOL_U_";
  std::string restart_v_name = "SOL_V_";

  // Yaml options
  bool is_loadYaml = true;
  std::string yaml_file("./linearPDE_input.yml");

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

 if(is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionReal("-rho", rho);
  SYS_T::GetOptionReal("-youngs_module", module_E);
  SYS_T::GetOptionReal("-poissons_ratio", nu);
  SYS_T::GetOptionInt("-nqp_vol", nqp_vol);
  SYS_T::GetOptionInt("-nqp_sur", nqp_sur);
  SYS_T::GetOptionReal("-rho_inf", genA_rho_inf);
  SYS_T::GetOptionBool("-is_backward_Euler", is_backward_Euler);
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionReal("-nl_rtol", nl_rtol);
  SYS_T::GetOptionReal("-nl_atol", nl_atol);
  SYS_T::GetOptionReal("-nl_dtol", nl_dtol);
  SYS_T::GetOptionInt("-nl_maxits", nl_maxits);
  SYS_T::GetOptionInt("-nl_refreq", nl_refreq);
  SYS_T::GetOptionInt("-nl_threshold", nl_threshold);
  SYS_T::GetOptionReal("-init_time", initial_time);
  SYS_T::GetOptionReal("-fina_time", final_time);
  SYS_T::GetOptionReal("-init_step", initial_step);
  SYS_T::GetOptionInt("-init_index", initial_index);
  SYS_T::GetOptionInt("-ttan_freq", ttan_renew_freq);
  SYS_T::GetOptionInt("-sol_rec_freq", sol_record_freq);
  SYS_T::GetOptionString("-sol_name", sol_bName);
  SYS_T::GetOptionBool("-is_restart", is_restart);
  SYS_T::GetOptionInt("-restart_index", restart_index);
  SYS_T::GetOptionReal("-restart_time", restart_time);
  SYS_T::GetOptionReal("-restart_step", restart_step);
  SYS_T::GetOptionString("-restart_u_name", restart_u_name);
  SYS_T::GetOptionString("-restart_v_name", restart_v_name);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint("-nqp_vol:", nqp_vol);
  SYS_T::cmdPrint("-nqp_sur:", nqp_sur);
  if( is_backward_Euler )
    SYS_T::commPrint(   "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(    "-rho_inf:",         genA_rho_inf);
  
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-nl_rtol:", nl_rtol);
  SYS_T::cmdPrint("-nl_atol:", nl_atol);
  SYS_T::cmdPrint("-nl_dtol:", nl_dtol);
  SYS_T::cmdPrint("-nl_maxits:", nl_maxits);
  SYS_T::cmdPrint("-nl_refreq:", nl_refreq);
  SYS_T::cmdPrint("-nl_threshold:", nl_threshold);
  SYS_T::cmdPrint("-init_time:", initial_time);
  SYS_T::cmdPrint("-init_step:", initial_step);
  SYS_T::cmdPrint("-init_index:", initial_index);
  SYS_T::cmdPrint("-fina_time:", final_time);
  SYS_T::cmdPrint("-ttan_freq:", ttan_renew_freq);
  SYS_T::cmdPrint("-sol_rec_freq:", sol_record_freq);
  SYS_T::cmdPrint("-sol_name:", sol_bName);
  if(is_restart)
  {
    SYS_T::commPrint("-is_restart: true \n");
    SYS_T::cmdPrint("-restart_index:", restart_index);
    SYS_T::cmdPrint("-restart_time:", restart_time);
    SYS_T::cmdPrint("-restart_step:", restart_step);
    SYS_T::cmdPrint("-restart_u_name:", restart_u_name);
    SYS_T::cmdPrint("-restart_v_name:", restart_v_name);
  }
  else SYS_T::commPrint("-is_restart: false \n");

  // ===== Record important solver options =====
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("solver_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

    cmdh5w->write_doubleScalar("init_step", initial_step);
    cmdh5w->write_intScalar("sol_record_freq", sol_record_freq);
    cmdh5w->write_intScalar("nqp_vol", nqp_vol);
    cmdh5w->write_intScalar("nqp_sur", nqp_sur);
    cmdh5w->write_doubleScalar("youngs_module", module_E);
    cmdh5w->write_doubleScalar("poissons_ratio", nu);

    delete cmdh5w; H5Fclose(cmd_file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  SYS_T::print_fatal_if( size != ANL_T::get_cpu_size(part_file, rank),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis.\n", size);

  auto locIEN = SYS_T::make_unique<ALocal_IEN>(part_file, rank);
  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_file, rank);
  auto fNode = SYS_T::make_unique<FEANode>(part_file, rank);
  auto pNode = SYS_T::make_unique<APart_Node>(part_file, rank);
  auto locnbc = SYS_T::make_unique<ALocal_NBC>(part_file, rank);
  auto locebc = SYS_T::make_unique<ALocal_EBC>(part_file, rank);

  // ===== Generate a sparse matrix for the enforcement of essential BCs
  auto pmat = SYS_T::make_unique<Matrix_PETSc>(pNode.get(), locnbc.get());

  pmat->gen_perm_bc(pNode.get(), locnbc.get());

  // ===== Generalized-alpha =====
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  auto tm_galpha = is_backward_Euler
    ? SYS_T::make_unique<TimeMethod_GenAlpha>(1.0, 1.0, 1.0)
    : SYS_T::make_unique<TimeMethod_GenAlpha>(genA_rho_inf, false);

  tm_galpha->print_info();

  // ===== Local Assembly Routine =====
  std::unique_ptr<IPLocAssem> locAssem_ptr = 
    SYS_T::make_unique<PLocAssem_Elastodynamics_GenAlpha>(
      ANL_T::get_elemType(part_file, rank), nqp_vol, nqp_sur,
      rho, module_E, nu, tm_galpha.get(), locebc->get_num_ebc());

  // ===== Initial condition =====
  std::unique_ptr<PDNSolution> disp =
    SYS_T::make_unique<PDNSolution_Elastodynamics>( pNode.get(), 0 );
  
  std::unique_ptr<PDNSolution> velo =
    SYS_T::make_unique<PDNSolution_Elastodynamics>( pNode.get(), 0 );

  std::unique_ptr<PDNSolution> dot_disp =
    SYS_T::make_unique<PDNSolution_Elastodynamics>( pNode.get(), 0 );

  std::unique_ptr<PDNSolution> dot_velo =
    SYS_T::make_unique<PDNSolution_Elastodynamics>( pNode.get(), 0 );
  
  if( is_restart )
  {
    initial_index = restart_index;
    initial_time  = restart_time;
    initial_step  = restart_step;

    // Read sol file
    SYS_T::file_check(restart_u_name);
    SYS_T::file_check(restart_v_name);

    disp->ReadBinary(restart_u_name);//str
    velo->ReadBinary(restart_v_name);

    // generate the corresponding dot_sol file name
    std::string restart_dot_u_name = "dot_";
    std::string restart_dot_v_name = "dot_";
    restart_dot_u_name.append(restart_u_name);
    restart_dot_v_name.append(restart_v_name);

    // Read dot_sol file
    SYS_T::file_check(restart_dot_u_name);
    SYS_T::file_check(restart_dot_v_name);

    dot_disp->ReadBinary(restart_dot_u_name);
    dot_velo->ReadBinary(restart_dot_v_name);

    SYS_T::commPrint("===> Read sol from disk as a restart run... \n");
    SYS_T::commPrint("     restart_u_name: %s \n", restart_u_name.c_str());
    SYS_T::commPrint("     restart_v_name: %s \n", restart_v_name.c_str());
    SYS_T::commPrint("     restart_dot_u_name: %s \n", restart_dot_u_name.c_str());
    SYS_T::commPrint("     restart_dot_v_name: %s \n", restart_dot_v_name.c_str());
    SYS_T::commPrint("     restart_time: %e \n", restart_time);
    SYS_T::commPrint("     restart_index: %d \n", restart_index);
    SYS_T::commPrint("     restart_step: %e \n", restart_step);
  }

  // ===== Time step info =====
  auto timeinfo = SYS_T::make_unique<PDNTimeStep>(initial_index, initial_time, initial_step);

  // ===== Global assembly =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  std::unique_ptr<IPGAssem> gloAssem =
    SYS_T::make_unique<PGAssem_LinearPDE_GenAlpha>(
      std::move(locIEN), std::move(locElem), std::move(fNode),
      std::move(pNode), std::move(locnbc), std::move(locebc),
      std::move(locAssem_ptr),
      ANL_T::get_nLocBas(part_file, rank), nz_estimate );  

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem->Assem_nonzero_estimate();

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem->Fix_nonzero_err_str();
  gloAssem->Clear_KG();
  
  // ===== Initialize the dot_sol vector by solving mass matrix =====
  if( is_restart == false )
  {
    auto lsolver_acce = SYS_T::make_unique<PLinear_Solver_PETSc>(
        1.0e-14, 1.0e-85, 1.0e30, 1000, "mass_", "mass_" );

    KSPSetType(lsolver_acce->ksp, KSPGMRES);
    KSPGMRESSetOrthogonalization(lsolver_acce->ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPGMRESSetRestart(lsolver_acce->ksp, 500);

    PC preproc; lsolver_acce->GetPC(&preproc);
    PCSetType( preproc, PCHYPRE );
    PCHYPRESetType( preproc, "boomeramg" );
   
    SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n"); 
    gloAssem->Assem_mass_residual( disp.get() );

    lsolver_acce->Solve( gloAssem->K, gloAssem->G, dot_velo.get() );
  
    dot_velo -> ScaleValue( -1.0 );

    // dot_u = v
    dot_disp -> Copy( velo.get() );

    SYS_T::commPrint("\n===> Consistent initial acceleration is obtained. \n");
    lsolver_acce -> print_info();
    SYS_T::commPrint(" The mass matrix lsolver is destroyed.\n");
  }

  // ===== Linear solver context =====
  auto lsolver = SYS_T::make_unique<PLinear_Solver_PETSc>();
  
  // ===== Nonlinear solver context =====
  auto nsolver = SYS_T::make_unique<PNonlinear_LinearPDE_Solver>(
      std::move(gloAssem), std::move(lsolver), std::move(pmat),
      std::move(tm_galpha),
      nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq, nl_threshold );

  nsolver->print_info();

  // ===== Temporal solver context =====
  auto tsolver = SYS_T::make_unique<PTime_LinearPDE_Solver>(
      std::move(nsolver), sol_bName,
      sol_record_freq, ttan_renew_freq, final_time );

  tsolver->print_info();
  
  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");

  tsolver->TM_GenAlpha_Elastodynamics(is_restart,
      std::move(dot_disp), std::move(dot_velo),
      std::move(disp), std::move(velo), std::move(timeinfo) );

  // ===== Print complete solver info =====
  tsolver -> print_lsolver_info();

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
