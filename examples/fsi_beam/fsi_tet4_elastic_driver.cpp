// ==================================================================
// fsi_tet4_elastic_driver.cpp
//
// 4-node tetrahedral finite element anslysis code for 3D FSI problem
// with the VMS formulation:
// 
// Elastostatic mesh update algorithm is invoked.
//
// Date: Oct. 12 2017
// Author: Ju Liu
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node_FSI.hpp"
#include "ALocal_Elem_wTag.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "ALocal_EBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "MaterialModel_NeoHookean_Incompressible_Mixed.hpp"
#include "MaterialModel_NeoHookean_M94_Mixed.hpp"
#include "PDNSolution_Mixed_UPV_3D.hpp"
#include "PLocAssem_Tet4_VMS_Seg_Incompressible.hpp"
#include "PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha.hpp"
#include "PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha.hpp"
#include "PLocAssem_Tet4_FSI_Mesh_Elastostatic.hpp"
#include "PTime_Seg_Solver.hpp"

int main( int argc, char * argv[] )
{
  int nqp_tet = 5, nqp_tri = 4;
  
  double genA_rho_inf = 0.5;

  std::string part_file("part");

  double nl_rtol = 1.0e-3;
  double nl_atol = 1.0e-6;
  double nl_dtol = 10.0;
  int nl_maxits = 20;
  int nl_refreq = 4;

  double initial_time = 0.0;
  double initial_step = 0.1;
  int initial_index = 0;
  double final_time = 1.0;
  std::string sol_bName("SOL_");
  int ttan_renew_freq = 1;
  int sol_record_freq = 1;

  bool is_restart = false;
  int restart_index = 0;
  double restart_time = 0.0;
  double restart_step = 1.0e-3;
  std::string restart_name = "SOL_";

  PetscMPIInt rank, size; 

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
 
  // ===== Command Line Argument ===== 
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-nqp_tet", nqp_tet);
  SYS_T::GetOptionInt("-nqp_tri", nqp_tri);
  SYS_T::GetOptionReal("-rho_inf", genA_rho_inf );
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionReal("-nl_rtol", nl_rtol);
  SYS_T::GetOptionReal("-nl_atol", nl_atol);
  SYS_T::GetOptionReal("-nl_dtol", nl_dtol);
  SYS_T::GetOptionInt("-nl_maxits", nl_maxits);
  SYS_T::GetOptionInt("-nl_refreq", nl_refreq);
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
  SYS_T::GetOptionString("-restart_name", restart_name);

  // ===== Record the command line arguments =====
  if( rank == 0 && !is_restart )
  {
    hid_t acmd_fid = H5Fcreate("analysis_cmd.h5",
        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(acmd_fid);

    cmdh5w->write_intScalar("nqp_tet", nqp_tet);
    cmdh5w->write_intScalar("nqp_tri", nqp_tri);
    cmdh5w->write_string("part_file", part_file);
    cmdh5w->write_doubleScalar("nl_rtol", nl_rtol);
    cmdh5w->write_doubleScalar("nl_atol", nl_atol);
    cmdh5w->write_doubleScalar("nl_dtol", nl_dtol);
    cmdh5w->write_intScalar("nl_maxits", nl_maxits);
    cmdh5w->write_intScalar("nl_refreq", nl_refreq);
    cmdh5w->write_doubleScalar("init_time", initial_time);
    cmdh5w->write_doubleScalar("init_step", initial_step);
    cmdh5w->write_doubleScalar("fina_time", final_time);

    delete cmdh5w; H5Fclose(acmd_fid);
  }

  // ===== Print the command line arguments on screen =====
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-nqp_tet:", nqp_tet);
  SYS_T::cmdPrint("-nqp_tri:", nqp_tri);
  SYS_T::cmdPrint("-rho_inf:", genA_rho_inf);
  SYS_T::cmdPrint("-nl_rtol:", nl_rtol);
  SYS_T::cmdPrint("-nl_atol:", nl_atol);
  SYS_T::cmdPrint("-nl_dtol:", nl_dtol);
  SYS_T::cmdPrint("-nl_maxits:", nl_maxits);
  SYS_T::cmdPrint("-nl_refreq:", nl_refreq);
  SYS_T::cmdPrint("-init_time:", initial_time);
  SYS_T::cmdPrint("-init_step:", initial_step);
  SYS_T::cmdPrint("-init_index:", initial_index);
  SYS_T::cmdPrint("-fina_time:", final_time);
  SYS_T::cmdPrint("-ttan_freq:", ttan_renew_freq);
  SYS_T::cmdPrint("-sol_rec_freq:", sol_record_freq);
  SYS_T::cmdPrint("-sol_name:", sol_bName);
  if(is_restart)
  {
    PetscPrintf(PETSC_COMM_WORLD, "-is_restart: true \n");
    SYS_T::cmdPrint("-restart_index:", restart_index);
    SYS_T::cmdPrint("-restart_time:", restart_time);
    SYS_T::cmdPrint("-restart_step:", restart_step);
    SYS_T::cmdPrint("-restart_name:", restart_name);
  }
  else PetscPrintf(PETSC_COMM_WORLD, "-is_restart: false \n");

  // ===== Main Data Structure =====
  SYS_T::commPrint("===> Reading mesh files ... ");
  FEANode * fNode = new FEANode(part_file, rank);
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  ALocal_Elem * locElem = new ALocal_Elem_wTag(part_file, rank);
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);
  ALocal_NodalBC * mesh_locnbc = new ALocal_NodalBC(part_file, rank, "mesh_nbc");
  ALocal_EBC * locebc = new ALocal_EBC(part_file, rank);
  ALocal_EBC * mesh_locebc = new ALocal_EBC(part_file, rank, "mesh_ebc");
  APart_Node * pNode = new APart_Node_FSI(part_file, rank, locElem, locIEN);
  SYS_T::commPrint("Done! \n");

  ALocal_Inflow_NodalBC * locinfnbc = NULL;
  ICVFlowRate * inflow_rate_ptr= NULL;

  // ===== Check the parallel setup =====
  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  // ===== Generate quadrature rules and FEM container =====
  PetscPrintf(PETSC_COMM_WORLD,
      "===> %d processor(s) are assigned for FEM analysis. \n", size);

  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv = new FEAElement_Tet4( quadv-> get_num_quadPts() );
  FEAElement * elements = new FEAElement_Triangle3_3D_der0(
      quads-> get_num_quadPts() );

  // ===== Generate a sparse matrix for strong enforcement of essential BC
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);
  pmat->gen_perm_bc(pNode, locnbc);

  Matrix_PETSc * mmat = new Matrix_PETSc(pNode, mesh_locnbc);
  mmat->gen_perm_bc(pNode, mesh_locnbc);

  // ===== Generate the Generalized-alpha method
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");
  const bool genA_is2ndSystem = false;
  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(
      genA_rho_inf, genA_is2ndSystem);
  tm_galpha_ptr->print_info();

  // ===== Solid material model =====
  SYS_T::commPrint("===> Setup the Material model.\n");
  const double mat_in_rho0 = 84.7458; 
  const double mat_in_E = 2.11865e9;

  IMaterialModel * matmodel = new MaterialModel_NeoHookean_M94_Mixed(
      mat_in_rho0, mat_in_E, 0.35 );

  IPLocAssem * locAssem_solid_ptr
    = new PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha(
        matmodel, tm_galpha_ptr, GMIptr->get_nLocBas(),
        quadv->get_num_quadPts(), elements->get_nLocBas() );

  // ===== Fluid material model =====
  IPLocAssem * locAssem_fluid_ptr = new PLocAssem_Tet4_ALE_VMS_NS_3D_GenAlpha(
      tm_galpha_ptr, GMIptr->get_nLocBas(),
      quadv->get_num_quadPts(), elements->get_nLocBas() );

  // ===== Mesh motion model =====
  // Mesh Young's modulus is chosen to be 1, nu to be 0.3
  IPLocAssem * locAssem_mesh_ptr = new PLocAssem_Tet4_FSI_Mesh_Elastostatic(
      1.0, 0.3 );

  // ===== Initial condition =====
  // case 1: steady inflow
  PDNSolution * sol = new PDNSolution_Mixed_UPV_3D( pNode, fNode, 1 );

  PDNSolution * dot_sol = new PDNSolution_Mixed_UPV_3D( pNode, fNode, 0 );

  if( is_restart )
  {
    initial_index = restart_index;
    initial_time  = restart_time;
    initial_step  = restart_step;

    SYS_T::file_exist_check(restart_name.c_str());

    sol->ReadBinary(restart_name.c_str());
    
    PetscPrintf(PETSC_COMM_WORLD, "===> Read sol from disk as a restart run... \n");
    PetscPrintf(PETSC_COMM_WORLD, "     restart_name: %s \n", restart_name.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "     restart_time: %e \n", restart_time);
    PetscPrintf(PETSC_COMM_WORLD, "     restart_index: %d \n", restart_index);
    PetscPrintf(PETSC_COMM_WORLD, "     restart_step: %e \n", restart_step);
  }

  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // ===== Global assembly routine =====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_FSI_FEM( locAssem_fluid_ptr,
      locAssem_solid_ptr, GMIptr, locElem, locIEN, pNode, locnbc, locebc );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_fluid_ptr, 
      locAssem_solid_ptr, locIEN, pNode, locnbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // ===== Global assembly for mesh equation =====
  SYS_T::commPrint("===> Initializing Mat K_mesh and Vec G_mesh ... \n");
  IPGAssem * gloAssem_mesh_ptr = new PGAssem_Seg_FEM( locAssem_mesh_ptr,
      GMIptr, locElem, locIEN, pNode, mesh_locnbc, mesh_locebc );

  SYS_T::commPrint("===> Assembly nonzero estimate for K_mesh ... \n");
  gloAssem_mesh_ptr->Assem_nonzero_estimate( locElem, locAssem_mesh_ptr, locIEN,
      pNode, mesh_locnbc );

  SYS_T::commPrint("===> Matrix K_mesh nonzero structure fixed. \n");
  gloAssem_mesh_ptr->Fix_nonzero_err_str();
  gloAssem_mesh_ptr->Clear_KG();

  // ===== Initialize the Gen-alpha dot-sol by mass matrix =====
  SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n");
  PLinear_Solver_PETSc * lsolver_acce = new PLinear_Solver_PETSc(
      1.0e-14, 1.0e-85, 1.0e30, 1000, "ls_mass_", "pc_mass_" );

  KSPSetType(lsolver_acce->ksp, KSPGMRES);
  KSPGMRESSetOrthogonalization(lsolver_acce->ksp,
      KSPGMRESModifiedGramSchmidtOrthogonalization);
  KSPGMRESSetRestart(lsolver_acce->ksp, 500);

  lsolver_acce->Monitor();

  PC preproc; lsolver_acce->GetPC(&preproc);
  PCSetType( preproc, PCASM ); PCASMSetOverlap(preproc, 1);

  lsolver_acce->Info();
  gloAssem_ptr->Assem_mass_residual( sol, locElem, locAssem_fluid_ptr, 
      locAssem_solid_ptr, elementv,
      elements, quadv, quads, locIEN, pNode, fNode, locnbc, locebc );

  // Solution for the mass-residual system
  PDNSolution * dot_pres_velo = new PDNSolution_P_V_Mixed_3D( pNode, fNode, 0 );

  lsolver_acce->Solve( gloAssem_ptr->K, gloAssem_ptr->G, dot_pres_velo);

  SEG_SOL_T::PlusAiPV(0.0, -1.0, -1.0, dot_pres_velo, dot_sol);

  SEG_SOL_T::PlusAiVPV(1.0, 0.0, 0.0, sol, dot_sol);

  delete lsolver_acce; delete dot_pres_velo;
  SYS_T::commPrint("\n===> Consistent initial acceleration is obtained.");
  SYS_T::commPrint("\n===> The mass matrix lsolver is destroyed. \n\n");

  // ===== Linear solver context =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  // Linear solver preconditioner field split
  PC upc; lsolver->GetPC(&upc);
  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields); // A_00 for velo
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);  // A_11 for pres

  // Mesh motion equation linear solver
  PLinear_Solver_PETSc * mesh_lsolver = new PLinear_Solver_PETSc(
      1.0e-5, 1.0e-55, 1.0e30, 500, "mesh_", "mesh_" );

  // Fix the discrete laplacian for the mesh_lsolver
  gloAssem_mesh_ptr->Assem_tangent_residual( sol, sol, 0.0,
      timeinfo->get_step(), locElem, locAssem_mesh_ptr, elementv,
      elements, quadv, quads, locIEN, pNode, fNode, mesh_locnbc,
      mesh_locebc );
  mesh_lsolver -> SetOperator( gloAssem_mesh_ptr->K ); 
  
  PC mesh_pc; mesh_lsolver->GetPC(&mesh_pc);
  PCFieldSplitSetBlockSize( mesh_pc,3 );

  mesh_lsolver->Info();
  SYS_T::commPrint("===> mesh solver LHS setted up.\n");

  // ===== Nonlinear solver context =====
  PNonlinear_Seg_Solver * nsolver = new PNonlinear_Seg_Solver(
      pNode, fNode, nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq);
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->print_info();

  // ===== Time solver context =====
  PTime_Seg_Solver * tsolver = new PTime_Seg_Solver( sol_bName,
      sol_record_freq, ttan_renew_freq, final_time );
  SYS_T::commPrint("===> Time marching solver setted up:\n");
  tsolver->print_info();

  // ===== FEM analysis =====
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver->TM_FSI_GenAlpha_ElasticStiffen(
      is_restart, dot_sol, sol, tm_galpha_ptr,
      timeinfo, inflow_rate_ptr,locElem, locIEN, pNode, fNode, locnbc, 
      locinfnbc, mesh_locnbc, 
      locebc, mesh_locebc, pmat, mmat, elementv, elements, quadv, quads,
      locAssem_fluid_ptr, locAssem_solid_ptr, locAssem_mesh_ptr, 
      gloAssem_ptr, gloAssem_mesh_ptr, lsolver, mesh_lsolver, nsolver);

  // Linear solver finally print all information
  lsolver -> Info();
  mesh_lsolver -> Info();

  // ===== Clean up allocations =====
  delete tsolver; delete nsolver; delete lsolver; delete mesh_lsolver;
  delete gloAssem_ptr; delete gloAssem_mesh_ptr;
  delete timeinfo; delete sol; delete dot_sol;
  delete locAssem_solid_ptr; 
  delete locAssem_fluid_ptr; 
  delete locAssem_mesh_ptr; 
  delete matmodel; 
  delete elementv; delete elements; delete quadv; delete quads;
  delete pmat; delete mmat; delete tm_galpha_ptr;
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  delete locElem; delete locnbc; delete mesh_locnbc;
  delete locebc; delete pNode; delete mesh_locebc;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
