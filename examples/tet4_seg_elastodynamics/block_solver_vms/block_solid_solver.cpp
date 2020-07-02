// ==================================================================
// block_solid_solver.cpp
//
// This is the code that I investigate the block solver method for 
// VMS-hyperelastodynamics. The element will be 4-node tet and 10-node
// tet.
//
// Author: Ju Liu
// Date Created: Feb 20 2018
// ==================================================================
#include "FEANode.hpp"
#include "ALocal_IEN.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node.hpp"
#include "ALocal_Elem_wTag.hpp"
#include "ALocal_EBC_wIntPts.hpp"
#include "ALocal_NodalBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "Matrix_PETSc.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "MaterialModel_NeoHookean_ST91_Mixed.hpp"
#include "PLocAssem_2x2Block_VMS_HE.hpp"
#include "PDNSolution_Disp_3D.hpp"
#include "PDNSolution_Pres_3D.hpp"
#include "PGAssem_2x2Block_HED.hpp"
#include "PDNTimeStep.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "PLinear_Solver_2x2Block_BIPN.hpp"
#include "PLinear_Solver_2x2Block_BIPN_SIMPLE.hpp"
#include "PTime_Solver_2x2Block_HED.hpp"

int main( int argc, char * argv[] )
{
  int nqp_tet = 5, nqp_tri = 4;
  std::string part_file("part");

  // BIPN parameters
  double bipn_rtol = 1.0e-2;
  int bipn_maxits = 10;
  int ainvcase = 0;

  // Nonlinear solver options 
  double nl_rtol = 1.0e-3;
  double nl_atol = 1.0e-6;
  double nl_dtol = 10.0;
  int nl_maxits = 20;
  int nl_refreq = 4;

  // Temporal solver options
  double initial_time = 0.0;
  double initial_step = 0.1;
  int initial_index = 0;
  double final_time = 1.0;
  std::string sol_d_bName("SOL_d_");
  std::string sol_p_bName("SOL_p_");
  std::string sol_v_bName("SOL_v_");
  int ttan_renew_freq = 1;
  int sol_record_freq = 1;

  // Restart options
  bool is_restart = false;
  int restart_index = 0;
  double restart_time = 0.0;
  double restart_step = 1.0e-3;
  std::string restart_d_name = "SOL_d_";
  std::string restart_p_name = "SOL_p_";
  std::string restart_v_name = "SOL_v_";

  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-nqp_tet", nqp_tet);
  SYS_T::GetOptionInt("-nqp_tri", nqp_tri);
  SYS_T::GetOptionReal("-bipn_rtol", bipn_rtol);
  SYS_T::GetOptionInt("-bipn_maxits", bipn_maxits);
  SYS_T::GetOptionInt("-bipn_ainv", ainvcase);
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
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionBool("-is_restart", is_restart);
  SYS_T::GetOptionInt("-restart_index", restart_index);
  SYS_T::GetOptionReal("-restart_time", restart_time);
  SYS_T::GetOptionReal("-restart_step", restart_step);

  SYS_T::cmdPrint("-nqp_tet:", nqp_tet);
  SYS_T::cmdPrint("-nqp_tri:", nqp_tri);
  SYS_T::cmdPrint("-bipn_rtol:", bipn_rtol);
  SYS_T::cmdPrint("-bipn_maxits:", bipn_maxits);
  SYS_T::cmdPrint("-bipn_ainv:", ainvcase);
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
  SYS_T::cmdPrint("-part_file", part_file);
  SYS_T::cmdPrint("sol_d_bName", sol_d_bName);
  SYS_T::cmdPrint("sol_p_bName", sol_p_bName);
  SYS_T::cmdPrint("sol_v_bName", sol_v_bName);
  if(is_restart)
  {
    PetscPrintf(PETSC_COMM_WORLD, "-is_restart: true \n");
    SYS_T::cmdPrint("-restart_index:", restart_index);
    SYS_T::cmdPrint("-restart_time:", restart_time);
    SYS_T::cmdPrint("-restart_step:", restart_step);
    SYS_T::cmdPrint("restart_d_name", restart_d_name);
    SYS_T::cmdPrint("restart_p_name", restart_p_name);
    SYS_T::cmdPrint("restart_v_name", restart_v_name);
  }
  else PetscPrintf(PETSC_COMM_WORLD, "-is_restart: false \n");

  // ===== Main Data Structure =====
  SYS_T::commPrint("===> Reading mesh files ... ");

  FEANode * fNode = new FEANode(part_file, rank);

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file, rank);

  ALocal_EBC_wIntPts * locebc = new ALocal_EBC_wIntPts(part_file, rank);

  //ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  APart_Node * pNode = new APart_Node(part_file, rank);

  ALocal_IEN * locIEN = new ALocal_IEN( part_file, rank );

  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  SYS_T::commPrint("Done! \n");

  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(), "Error: Assigned CPU number does not match the partition. \n");

  delete PartBasic;

  PetscPrintf(PETSC_COMM_WORLD, "===> %d processor(s) are assigned for FEM analysis. \n", size);

  // FEA objects
  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv = new FEAElement_Tet4( nqp_tet );
  FEAElement * elements = new FEAElement_Triangle3_3D_der0( nqp_tri );

  // Matrix for Essential BC
  //Matrix_PETSc * pmat = new Matrix_PETSc( pNode, locnbc );
  //pmat -> gen_perm_bc( pNode, locnbc );

  // Gen-alpha objects
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");
  const double genA_spectrium = 0.5;
  const bool genA_is2ndSystem = false;
  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(
      genA_spectrium, genA_is2ndSystem);
  tm_galpha_ptr->print_info();

  // Material Model
  SYS_T::commPrint("===> Setup the Material model.\n");
  const double mat_in_rho = 1.0;
  const double mat_in_E = 2.405e8;
  const double mat_in_nu = 0.4998;
  IMaterialModel * matmodel = new MaterialModel_NeoHookean_ST91_Mixed(
      mat_in_rho, mat_in_E, mat_in_nu );

  IPLocAssem_2x2Block * locAssem_ptr = new PLocAssem_2x2Block_VMS_HE(
      matmodel, tm_galpha_ptr, GMIptr->get_nLocBas(),
      quadv->get_num_quadPts(), elements->get_nLocBas() );

  PDNSolution * disp = new PDNSolution_Disp_3D(pNode, fNode, 0);
  PDNSolution * pres = new PDNSolution_Pres_3D(pNode, fNode, 0);
  PDNSolution * velo = new PDNSolution_Disp_3D(pNode, fNode, 0);

  PDNSolution * dot_disp = new PDNSolution_Disp_3D(pNode, fNode, 0);
  PDNSolution * dot_pres = new PDNSolution_Pres_3D(pNode, fNode, 0);
  PDNSolution * dot_velo = new PDNSolution_Disp_3D(pNode, fNode, 0);

  if( is_restart )
  {
    initial_index = restart_index;
    initial_time  = restart_time;
    initial_step  = restart_step;

    SYS_T::file_check(restart_d_name.c_str());
    SYS_T::file_check(restart_p_name.c_str());
    SYS_T::file_check(restart_v_name.c_str());

    disp->ReadBinary(restart_d_name.c_str());
    pres->ReadBinary(restart_p_name.c_str());
    velo->ReadBinary(restart_v_name.c_str());
    
    PetscPrintf(PETSC_COMM_WORLD, "===> Read sol from disk as a restart run... \n");
    PetscPrintf(PETSC_COMM_WORLD, "     restart_d_name: %s \n", restart_d_name.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "     restart_p_name: %s \n", restart_p_name.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "     restart_v_name: %s \n", restart_v_name.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "     restart_time: %e \n", restart_time);
    PetscPrintf(PETSC_COMM_WORLD, "     restart_index: %d \n", restart_index);
    PetscPrintf(PETSC_COMM_WORLD, "     restart_step: %e \n", restart_step);
  }

  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // Global Assembly Routine
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem_2x2Block * gloAssem_ptr = new PGAssem_2x2Block_HED(
      locAssem_ptr, locIEN, pNode, locnbc, locebc );

  gloAssem_ptr->Assem_nonzero_estimate( locAssem_ptr, locIEN, locnbc );
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();
  
  // Mass matrix solver
  SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n");
  PLinear_Solver_PETSc * ls_mass0 = new PLinear_Solver_PETSc(
      1.0e-14, 1.0e-85, 1.0e30, 1000, "ls_mass0_", "pc_mass0_" );
  PLinear_Solver_PETSc * ls_mass1 = new PLinear_Solver_PETSc(
      1.0e-14, 1.0e-85, 1.0e30, 1000, "ls_mass1_", "pc_mass1_" );

  KSPSetType(ls_mass0->ksp, KSPGMRES);
  KSPGMRESSetOrthogonalization(ls_mass0->ksp,
      KSPGMRESModifiedGramSchmidtOrthogonalization);
  KSPGMRESSetRestart(ls_mass0->ksp, 500);

  KSPSetType(ls_mass1->ksp, KSPGMRES);
  KSPGMRESSetOrthogonalization(ls_mass1->ksp,
      KSPGMRESModifiedGramSchmidtOrthogonalization);
  KSPGMRESSetRestart(ls_mass1->ksp, 500);

  ls_mass0->Monitor();
  ls_mass1->Monitor();

  PC lsm0_pc; ls_mass0->GetPC( &lsm0_pc );
  PCSetType( lsm0_pc, PCML );

  PC lsm1_pc; ls_mass1->GetPC( &lsm1_pc );
  PCSetType( lsm1_pc, PCML );

  gloAssem_ptr -> Assem_mass_residual( disp, pres, velo, locAssem_ptr, 
      elementv, elements, quadv, quads, locIEN, fNode, locnbc, locebc );

  // Initialize dot_pres
  ls_mass0->Solve( gloAssem_ptr->K_00, gloAssem_ptr->G_0, dot_pres );
  SYS_T::commPrint("\n\n");
  // Initialize dot_velo
  ls_mass1->Solve( gloAssem_ptr->K_11, gloAssem_ptr->G_1, dot_velo );
 
  // dot_disp = velo
  dot_disp -> ScaleValue( 0.0 );
  dot_disp -> PlusAX( velo, 1.0 );

  delete ls_mass0; delete ls_mass1;
  SYS_T::commPrint("\n===> Consistent initial acceleration is obtained.");
  SYS_T::commPrint("\n===> The mass matrix lsolver is destroyed. \n\n");

  // FEM.6 Block Linear Solver
  const double ll_rtol = 1.0e-5;
  const double ll_atol = 1.0e-50;
  const double ll_dtol = 1.0e20;
  const int ll_maxits = 200;
  
  IPLinear_Solver_2x2Block * lsolver = new PLinear_Solver_2x2Block_BIPN_SIMPLE(
      ll_rtol, ll_atol, ll_dtol, ll_maxits,
      ll_rtol, ll_atol, ll_dtol, ll_maxits,
      bipn_rtol, bipn_maxits, gloAssem_ptr, ainvcase ); 

  // FEM.7 Nonlinear solver
  PNonlinear_Solver_2x2Block_HED * nsolver = 
    new PNonlinear_Solver_2x2Block_HED( pNode, fNode, nl_rtol, nl_atol,
      nl_dtol, nl_maxits, nl_refreq );
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->print_info();

  // FEM.8 Time solver
  PTime_Solver_2x2Block_HED * tsolver =
    new PTime_Solver_2x2Block_HED( pNode, fNode, 
        sol_d_bName, sol_p_bName, sol_v_bName, 
        sol_record_freq, ttan_renew_freq, final_time );
  SYS_T::commPrint("===> Time marching solver setted up:\n");
  tsolver->print_info();

  // FEM.9 FEM analysis
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver->TM_GenAlpha_Solve(is_restart, dot_disp, dot_pres, dot_velo, 
      disp, pres, velo, tm_galpha_ptr, timeinfo, locIEN, fNode, 
      locnbc, locebc, elementv, elements, quadv, quads,
      locAssem_ptr, gloAssem_ptr, lsolver, nsolver);

  lsolver -> print_info();
  
  // ===== PETSc Finalize =====
  delete tsolver;  delete nsolver; delete lsolver;
  delete gloAssem_ptr; delete timeinfo;
  delete disp; delete pres; delete velo;
  delete dot_disp; delete dot_pres; delete dot_velo;
  delete locAssem_ptr; delete tm_galpha_ptr; delete matmodel; // delete pmat;
  delete elementv; delete elements; delete quadv; delete quads; //delete locElem;
  delete fNode; delete GMIptr; delete locebc; delete pNode; delete locIEN;
  delete locnbc;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
