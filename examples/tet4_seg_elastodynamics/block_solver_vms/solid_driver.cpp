// solid_driver.cpp
// local projection stabilization FEM for solid dynamics
// lowest order LPS.
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
#include "MaterialModel_NeoHookean_ST91_Mixed.hpp"
#include "PLocAssem_LPS0_Seg_GenAlpha.hpp"
#include "PDNSolution_Mixed_U_Hyperelastic_3D.hpp"
#include "PDNSolution_P_V_Mixed_Hyperelastic_3D.hpp"
#include "PGAssem_Seg_FEM.hpp"
#include "PTime_Seg_Solver.hpp"

int main(int argc, char *argv[])
{
  int nqp_tet = 29, nqp_tri = 13;
  const std::string part_file("part");

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

  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-nqp_tet", nqp_tet);
  SYS_T::GetOptionInt("-nqp_tri", nqp_tri);
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

  SYS_T::cmdPrint("-nqp_tet:", nqp_tet);
  SYS_T::cmdPrint("-nqp_tri:", nqp_tri);
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
  
  // ===== Main Data Structure =====
  SYS_T::commPrint("===> Reading mesh files ... ");
  
  FEANode * fNode = new FEANode(part_file, rank);
  
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file, rank);
  
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  
  ALocal_EBC_wIntPts * locebc = new ALocal_EBC_wIntPts(part_file, rank);
  
  APart_Node * pNode = new APart_Node(part_file, rank);
  
  ALocal_IEN * locIEN = new ALocal_IEN( part_file, rank );
  
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);
  
  SYS_T::commPrint("Done! \n");
  
  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  PetscPrintf(PETSC_COMM_WORLD,
      "===> %d processor(s) are assigned for FEM analysis. \n", size);

  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv = new FEAElement_Tet4( nqp_tet );
  FEAElement * elements = new FEAElement_Triangle3_3D_der0( nqp_tri );
  
  Matrix_PETSc * pmat = new Matrix_PETSc( pNode, locnbc );
  pmat -> gen_perm_bc( pNode, locnbc );
  
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");
  const double genA_spectrium = 0.5;
  const bool genA_is2ndSystem = false;
  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(
      genA_spectrium, genA_is2ndSystem);
  tm_galpha_ptr->print_info();

  SYS_T::commPrint("===> Setup the Material model.\n");
  const double mat_in_rho = 1000.0;
  const double mat_in_E = 2.405e7;
  const double mat_in_nu = 0.4998;
  IMaterialModel * matmodel = new MaterialModel_NeoHookean_ST91_Mixed(
      mat_in_rho, mat_in_E, mat_in_nu );

  // Local assembly method
  IPLocAssem * locAssem_ptr = new PLocAssem_LPS0_Seg_GenAlpha(
      matmodel, tm_galpha_ptr, GMIptr->get_nLocBas(),
      quadv->get_num_quadPts(), elements->get_nLocBas() );

  // Initial condition
  PDNSolution * sol = new PDNSolution_Mixed_U_Hyperelastic_3D( pNode,
      locAssem_ptr, fNode, 0 );

  PDNSolution * dot_sol = new PDNSolution_Mixed_U_Hyperelastic_3D( pNode,
      locAssem_ptr, fNode, 0 );

  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // Global asssembly routine
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_Seg_FEM( locAssem_ptr,
      GMIptr, locElem, locIEN, pNode, locnbc, locebc );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr, locIEN,
      pNode, locnbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // Initial condition 
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
  gloAssem_ptr->Assem_mass_residual( sol, locElem, locAssem_ptr, elementv,
      elements, quadv, quads, locIEN, pNode, fNode, locnbc, locebc );

  PDNSolution * dot_pres_velo = new PDNSolution_P_V_Mixed_Hyperelastic_3D( pNode,
      locAssem_ptr, fNode, 0 );
  
  lsolver_acce->Solve( gloAssem_ptr->K, gloAssem_ptr->G, dot_pres_velo);
  
  SEG_SOL_T::PlusAiPV(0.0, -1.0, -1.0, dot_pres_velo, dot_sol);
  SEG_SOL_T::PlusAiVPV(1.0, 0.0, 0.0, sol, dot_sol);
  delete lsolver_acce; delete dot_pres_velo;
  SYS_T::commPrint("\n===> Consistent initial acceleration is obtained.");
  SYS_T::commPrint("\n===> The mass matrix lsolver is destroyed. \n\n");

  // Linear solver context
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();
  lsolver->Info();

  PC upc; lsolver->GetPC(&upc);
  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields); // A_00 for velo
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);  // A_11 for pres

  // Nonlinear solver context
  PNonlinear_Seg_Solver * nsolver = new PNonlinear_Seg_Solver(
      pNode, locAssem_ptr, fNode,
      nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq);
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->print_info();

  // Time solver context
  PTime_Seg_Solver * tsolver = new PTime_Seg_Solver( sol_bName,
      sol_record_freq, ttan_renew_freq, final_time );
  SYS_T::commPrint("===> Time marching solver setted up:\n");
  tsolver->print_info();

  // FEM analysis
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  const bool is_restart = false;
  tsolver->TM_Seg_GenAlpha(is_restart, dot_sol, sol, tm_galpha_ptr,
      timeinfo, locElem, locIEN, pNode, fNode, locnbc, locebc, pmat,
      elementv, elements, quadv, quads,
      locAssem_ptr, gloAssem_ptr, lsolver, nsolver);

  // === Clean up memory allocation ===
  delete tsolver; delete nsolver; delete lsolver; delete gloAssem_ptr;
  delete timeinfo; delete sol; delete dot_sol; delete locAssem_ptr;
  delete pmat; delete tm_galpha_ptr; delete matmodel;
  delete elementv; delete elements; delete quadv; delete quads;
  delete pNode; delete locebc; delete locnbc; delete GMIptr;
  delete PartBasic; delete locElem; delete locIEN; delete fNode;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
