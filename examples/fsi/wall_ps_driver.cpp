// ============================================================================
// wall_ps_tet4_driver.cpp
//
// Wall mechanics solver for generating the prestress field.
//
// Date: Jan 28 2022
// ============================================================================
#include "HDF5_Tools.hpp"
#include "ANL_Tools.hpp"
#include "MaterialModel_vol_Incompressible.hpp"
#include "MaterialModel_vol_ST91.hpp"
#include "MaterialModel_vol_M94.hpp"
#include "MaterialModel_ich_NeoHookean.hpp"
#include "MaterialModel_ich_GOH14.hpp"
#include "PLocAssem_2x2Block_VMS_Incompressible.hpp"
#include "PLocAssem_2x2Block_VMS_Hyperelasticity.hpp"
#include "PGAssem_Wall_Prestress.hpp"
#include "PTime_FSI_Solver.hpp"

int main( int argc, char *argv[] )
{
  // solution file name to be loaded for prestressing
  std::string restart_velo_name = "SOL_velo_re";
  std::string restart_pres_name = "SOL_pres_re";

  // (Pseudo-) time integration parameters
  double genA_rho_inf = 0.0;
  bool is_backward_Euler = true;
  const bool is_load_ps = false;

  // Estimate of num nonzeros per row for the sparse tangent matrix
  int nz_estimate = 300;

  // Prestress tolerance
  double prestress_disp_tol = 1.0e-6;

  // Nonlinear solver parameters
  double nl_rtol    = 1.0e-3;        // convergence criterion relative tolerance
  double nl_atol    = 1.0e-6;        // convergence criterion absolute tolerance
  double nl_dtol    = 1.0e3;         // divergence criterion
  int    nl_maxits  = 20;            // maximum number if nonlinear iterations
  int    nl_refreq  = 4;             // frequency of tangent matrix renewal
  int    nl_threshold = 4;             // threshold of tangent matrix renewal

  // Time stepping parameters
  double initial_time = 0.0;         // time of initial condition
  double initial_step = 0.1;         // time step size
  int    initial_index = 0;          // index of initial condition
  double final_time = 1.0;           // end time of simulation
  bool   is_record_sol = false;      // bool flag to decide if one wants to record the solution
  std::string sol_bName("PS_");      // base name of the solution file
  int    ttan_renew_freq = 1;        // frequency of tangent matrix renewal
  int    sol_record_freq = 1;        // frequency for recording the solution

  // Solid properties
  double solid_density = 1.0;
  double solid_E = 2.0e6;
  double solid_nu = 0.5;

  // We assume that a 3D solver has been called (to generate the wall traction)
  // and a suite of command line arguments has been saved to disk
  hid_t solver_cmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( solver_cmd_file );

  const int nqp_vol     = cmd_h5r -> read_intScalar(    "/", "nqp_vol");
  const int nqp_sur     = cmd_h5r -> read_intScalar(    "/", "nqp_sur");

  delete cmd_h5r; H5Fclose(solver_cmd_file);

  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5_Reader * pcmd_h5r = new HDF5_Reader( prepcmd_file );

  const std::string part_v_file = pcmd_h5r -> read_string(    "/", "part_file_v" );
  const std::string part_p_file = pcmd_h5r -> read_string(    "/", "part_file_p" );
  const int fsiBC_type          = pcmd_h5r -> read_intScalar( "/", "fsiBC_type" );
  const std::string elemType_str = pcmd_h5r -> read_string("/","elemType");
  const FEType elemType = FE_T::to_FEType(elemType_str);

  delete pcmd_h5r; H5Fclose(prepcmd_file);

  // Initialize PETSc
#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  // Assert that the fsiBC type is 2, which clamped the lumen nodes
  SYS_T::print_fatal_if( fsiBC_type != 2, "Error: fsiBC_type should be 2.\n" );

  // Clean potentially pre-existing hdf5 files of prestress saved in the folder
  // named as ps_data
  if(rank == 0 )
  {
    if( SYS_T::directory_exist("ps_data") )
    {
      std::cout<<"Clean the folder ps_data.\n";
      SYS_T::execute("rm -rf ps_data");
    }

    SYS_T::execute("mkdir ps_data");
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  const std::string ps_file_name("./ps_data/prestress");

  SYS_T::GetOptionString("-restart_velo_name",   restart_velo_name);
  SYS_T::GetOptionString("-restart_pres_name",   restart_pres_name);
  SYS_T::GetOptionReal(  "-rho_inf",             genA_rho_inf);
  SYS_T::GetOptionBool(  "-is_backward_Euler",   is_backward_Euler);
  SYS_T::GetOptionInt(   "-nz_estimate",         nz_estimate);
  SYS_T::GetOptionReal(  "-prestress_disp_tol",  prestress_disp_tol);
  SYS_T::GetOptionReal(  "-nl_rtol",             nl_rtol);
  SYS_T::GetOptionReal(  "-nl_atol",             nl_atol);
  SYS_T::GetOptionReal(  "-nl_dtol",             nl_dtol);
  SYS_T::GetOptionInt(   "-nl_maxits",           nl_maxits);
  SYS_T::GetOptionInt(   "-nl_refreq",           nl_refreq);
  SYS_T::GetOptionInt(   "-nl_rethred",          nl_threshold);
  SYS_T::GetOptionBool(  "-is_backward_Euler",   is_backward_Euler);
  SYS_T::GetOptionReal(  "-init_time",           initial_time);
  SYS_T::GetOptionReal(  "-fina_time",           final_time);
  SYS_T::GetOptionReal(  "-init_step",           initial_step);
  SYS_T::GetOptionInt(   "-init_index",          initial_index);
  SYS_T::GetOptionInt(   "-ttan_freq",           ttan_renew_freq);
  SYS_T::GetOptionBool(  "-is_record_sol",       is_record_sol);
  SYS_T::GetOptionInt(   "-sol_rec_freq",        sol_record_freq);
  SYS_T::GetOptionReal(  "-sl_density",          solid_density);
  SYS_T::GetOptionReal(  "-sl_E",                solid_E);
  SYS_T::GetOptionReal(  "-sl_nu",               solid_nu);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint(      "part_v_file:",          part_v_file);
  SYS_T::cmdPrint(      "part_p_file:",          part_p_file);
  SYS_T::cmdPrint(       "-prestress_disp_tol:", prestress_disp_tol);
  SYS_T::cmdPrint(       "-nl_rtol:",            nl_rtol);
  SYS_T::cmdPrint(       "-nl_atol:",            nl_atol);
  SYS_T::cmdPrint(       "-nl_dtol:",            nl_dtol);
  SYS_T::cmdPrint(       "-nl_maxits:",          nl_maxits);
  SYS_T::cmdPrint(       "-nl_refreq:",          nl_refreq);
  SYS_T::cmdPrint(       "-nl_rethred",          nl_threshold);

  if( is_backward_Euler )
    SYS_T::commPrint(    "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(     "-rho_inf:",            genA_rho_inf);

  SYS_T::cmdPrint(       "-init_time:",          initial_time);
  SYS_T::cmdPrint(       "-init_step:",          initial_step);
  SYS_T::cmdPrint(       "-init_index:",         initial_index);
  SYS_T::cmdPrint(       "-fina_time:",          final_time);
  SYS_T::cmdPrint(       "-ttan_freq:",          ttan_renew_freq);

  if( is_record_sol )
    SYS_T::cmdPrint(     "-sol_rec_freq:",       sol_record_freq);
  else
    SYS_T::commPrint(    "-is_record_sol: false \n");

  SYS_T::cmdPrint(     "-sl_density:",         solid_density);
  SYS_T::cmdPrint(     "-sl_E:",               solid_E);
  SYS_T::cmdPrint(     "-sl_nu:",              solid_nu);

  // ====== Record important parameters ======
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("wall_ps_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);
    
    cmdh5w -> write_string(      "ps_file_name",       ps_file_name);
    cmdh5w -> write_doubleScalar("prestress_disp_tol", prestress_disp_tol );

    delete cmdh5w; H5Fclose(cmd_file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ====== Data for Analysis ======
  auto fNode = SYS_T::make_unique<FEANode>(part_v_file, rank);

  auto locIEN_v = SYS_T::make_unique<ALocal_IEN>(part_v_file, rank);

  auto locIEN_p = SYS_T::make_unique<ALocal_IEN>(part_p_file, rank);

  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_v_file, rank);

  std::unique_ptr<APart_Node> pNode_v = SYS_T::make_unique<APart_Node_FSI>(part_v_file, rank);

  std::unique_ptr<APart_Node> pNode_p = SYS_T::make_unique<APart_Node_FSI>(part_p_file, rank);

  std::unique_ptr<APart_Node> pNode_v_time = SYS_T::make_unique<APart_Node_FSI>(part_v_file, rank);

  std::unique_ptr<APart_Node> pNode_p_time = SYS_T::make_unique<APart_Node_FSI>(part_p_file, rank);

  std::unique_ptr<APart_Node> pNode_v_nlinear = SYS_T::make_unique<APart_Node_FSI>(part_v_file, rank);

  auto locebc_v = SYS_T::make_unique<ALocal_EBC>(part_v_file, rank);

  auto locebc_p = SYS_T::make_unique<ALocal_EBC>(part_p_file, rank);

  auto locnbc_v = SYS_T::make_unique<ALocal_NBC>(part_v_file, rank, "/nbc/MF");

  auto locnbc_p = SYS_T::make_unique<ALocal_NBC>(part_p_file, rank, "/nbc/MF");

  auto ps_data = SYS_T::make_unique<Tissue_prestress>(locElem.get(), nqp_vol, rank, is_load_ps, ps_file_name);
 
  SYS_T::commPrint("===> Mesh HDF5 files are read from disk.\n");

  // Group APart_Node and ALocal_NBC into a vector
  std::vector<APart_Node *> pNode_list { pNode_v.get(), pNode_p.get() };

  std::vector<ALocal_NBC *> locnbc_list { locnbc_v.get(), locnbc_p.get() };

  std::vector<APart_Node *> pNode_m_list { pNode_v.get() };

  // ===== Basic Checking =====
  SYS_T::print_fatal_if( size!= ANL_T::get_cpu_size(part_v_file, rank),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis. \n", size);

  // ===== Generate the IS for pres and velo =====
  const int idx_v_start = HDF5_T::read_intScalar( SYS_T::gen_partfile_name(part_v_file, rank).c_str(), "/DOF_mapper", "start_idx" );
  const int idx_p_start = HDF5_T::read_intScalar( SYS_T::gen_partfile_name(part_p_file, rank).c_str(), "/DOF_mapper", "start_idx" );

  const int idx_v_len = pNode_v->get_dof() * pNode_v -> get_nlocalnode();
  const int idx_p_len = pNode_p->get_dof() * pNode_p -> get_nlocalnode();

  auto is_array_velo = SYS_T::make_unique<PetscInt>( static_cast<size_t>(idx_v_len) );
  for (int ii = 0; ii < idx_v_len; ++ii) is_array_velo[ii] = idx_v_start + ii;

  auto is_array_pres = SYS_T::make_unique<PetscInt>( static_cast<size_t>(idx_p_len) );
  for (int ii = 0; ii < idx_p_len; ++ii) is_array_pres[ii] = idx_p_start + ii;

  IS is_velo, is_pres;
  ISCreateGeneral(PETSC_COMM_WORLD, idx_v_len, is_array_velo.get(), PETSC_COPY_VALUES, &is_velo);
  ISCreateGeneral(PETSC_COMM_WORLD, idx_p_len, is_array_pres.get(), PETSC_COPY_VALUES, &is_pres);

  // ===== Generate a sparse matrix for strong enforcement of essential BC
  std::vector<int> start_idx{ idx_v_start, idx_p_start };

  auto pmat = SYS_T::make_unique<Matrix_PETSc>( idx_v_len + idx_p_len );
  pmat->gen_perm_bc( pNode_list, locnbc_list, start_idx );

  // ===== Generate the generalized-alpha method
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  auto tm_galpha = is_backward_Euler
    ? SYS_T::make_unique<TimeMethod_GenAlpha>(1.0, 1.0, 1.0)
    : SYS_T::make_unique<TimeMethod_GenAlpha>(genA_rho_inf, false);

  tm_galpha->print_info();

  // ===== Local assembly =====
  std::unique_ptr<IPLocAssem_2x2Block> locAssem_solid = nullptr;

  const double solid_mu = solid_E/(2.0+2.0*solid_nu);
  std::unique_ptr<IMaterialModel_ich> imodel = SYS_T::make_unique<MaterialModel_ich_NeoHookean>(solid_mu);

  if( solid_nu == 0.5 )
  {
    std::unique_ptr<IMaterialModel_vol> vmodel = SYS_T::make_unique<MaterialModel_vol_Incompressible>(solid_density);
    std::unique_ptr<MaterialModel_Mixed_Elasticity> matmodel = SYS_T::make_unique<MaterialModel_Mixed_Elasticity>(std::move(vmodel), std::move(imodel));
    locAssem_solid = SYS_T::make_unique<PLocAssem_2x2Block_VMS_Incompressible>(
        elemType, nqp_vol, nqp_sur, tm_galpha.get(), std::move(matmodel) );
  }
  else
  {
    const double solid_lambda = solid_nu * solid_E / ((1+solid_nu) * (1-2.0*solid_nu));
    const double solid_kappa  = solid_lambda + 2.0 * solid_mu / 3.0;
    std::unique_ptr<IMaterialModel_vol> vmodel = SYS_T::make_unique<MaterialModel_vol_M94>(solid_density, solid_kappa);
    std::unique_ptr<MaterialModel_Mixed_Elasticity> matmodel = SYS_T::make_unique<MaterialModel_Mixed_Elasticity>(std::move(vmodel), std::move(imodel));
    locAssem_solid = SYS_T::make_unique<PLocAssem_2x2Block_VMS_Hyperelasticity>(
        elemType, nqp_vol, nqp_sur, tm_galpha.get(), std::move(matmodel) );
  }

  // ===== Initial conditions =====
  std::unique_ptr<PDNSolution> velo =
    SYS_T::make_unique<PDNSolution_V>( pNode_v.get(), 0, true, "velo" );

  std::unique_ptr<PDNSolution> disp =
    SYS_T::make_unique<PDNSolution_V>( pNode_v.get(), 0, true, "disp" );

  std::unique_ptr<PDNSolution> pres =
    SYS_T::make_unique<PDNSolution_P>( pNode_p.get(), 0, true, "pres" );  

  std::unique_ptr<PDNSolution> dot_velo =
    SYS_T::make_unique<PDNSolution_V>( pNode_v.get(), 0, true, "dot_velo" );

  std::unique_ptr<PDNSolution> dot_disp =
    SYS_T::make_unique<PDNSolution_V>( pNode_v.get(), 0, true, "dot_disp" );

  std::unique_ptr<PDNSolution> dot_pres =
    SYS_T::make_unique<PDNSolution_P>( pNode_p.get(), 0, true, "dot_pres" ); 

  // Read sol file
  SYS_T::file_check(restart_velo_name);
  velo->ReadBinary(restart_velo_name);

  SYS_T::file_check(restart_pres_name);
  pres->ReadBinary(restart_pres_name);

  // Read dot_sol file
  std::string restart_dot_velo_name = "dot_";
  restart_dot_velo_name.append(restart_velo_name);
  SYS_T::file_check(restart_dot_velo_name);
  dot_velo->ReadBinary(restart_dot_velo_name);
 
  std::string restart_dot_pres_name = "dot_";
  restart_dot_pres_name.append(restart_pres_name);
  SYS_T::file_check(restart_dot_pres_name);
  dot_pres->ReadBinary(restart_dot_pres_name);

  SYS_T::commPrint("===> Read sol from disk as a restart run: \n");
  SYS_T::commPrint("     restart_velo_name:     %s \n", restart_velo_name.c_str());
  SYS_T::commPrint("     restart_dot_velo_name: %s \n", restart_dot_velo_name.c_str());
  SYS_T::commPrint("     restart_pres_name:     %s \n", restart_pres_name.c_str());
  SYS_T::commPrint("     restart_dot_pres_name: %s \n", restart_dot_pres_name.c_str());

  // ===== Time step info =====
  auto timeinfo = SYS_T::make_unique<PDNTimeStep>(initial_index, initial_time, 
      initial_step);

  // ===== Global assembly routine =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  std::unique_ptr<IPGAssem> gloAssem = SYS_T::make_unique<PGAssem_Wall_Prestress>(
      std::move(locIEN_v), std::move(locIEN_p), std::move(locElem), std::move(fNode), 
      std::move(pNode_v), std::move(pNode_p), std::move(locnbc_v), std::move(locnbc_p), 
      std::move(locebc_v), std::move(locebc_p), std::move(locAssem_solid), 
      std::move(ps_data), nz_estimate );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem->Assem_nonzero_estimate();

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem->Fix_nonzero_err_str();
  gloAssem->Clear_KG();

  // ===== Linear and nonlinear solver context =====
  auto lsolver = SYS_T::make_unique<PLinear_Solver_PETSc>();

  PC upc; lsolver->GetPC(&upc);
  PCFieldSplitSetIS(upc, "u", is_velo);
  PCFieldSplitSetIS(upc, "p", is_pres);

  // ===== Nonlinear solver context =====
  auto nsolver = SYS_T::make_unique<PNonlinear_FSI_Solver>(
      std::move(lsolver), std::move(pmat), std::move(tm_galpha), 
      std::move(pNode_v_nlinear), nl_rtol, nl_atol, nl_dtol, 
      nl_maxits, nl_refreq, nl_threshold );
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->print_info();

  // ===== Temporal solver context =====
  auto tsolver = SYS_T::make_unique<PTime_FSI_Solver>(
      std::move(nsolver), std::move(pNode_v_time), std::move(pNode_p_time), 
      sol_bName, sol_record_freq, ttan_renew_freq, final_time );
  SYS_T::commPrint("===> Time marching solver setted up:\n");
  tsolver->print_info();

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver -> TM_FSI_Prestress( is_record_sol, prestress_disp_tol, is_velo, is_pres,
      std::move(dot_disp), std::move(dot_velo), std::move(dot_pres), std::move(disp), std::move(velo), std::move(pres),
      std::move(timeinfo), gloAssem.get() );

  // ===== Record the wall prestress to h5 file =====
  ps_data -> write_prestress_hdf5();

  // Print complete solver info
  tsolver -> print_lsolver_info();

  // ==========================================================================
  // Clean the memory
  ISDestroy(&is_velo); ISDestroy(&is_pres);
  tsolver.reset(); gloAssem.reset();
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF 
