// ==================================================================
// seg_solid_driver.cpp
// This solver uses diagonal scaling for the matrices.
// and reads in the preprocessor info from a Gmsh geometry output.
//
// Date: Mar. 9 2018
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "ALocal_EBC_wIntPts.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Tet10.hpp"
#include "Matrix_PETSc.hpp"
#include "MaterialModel_NeoHookean_ST91_Mixed.hpp"
#include "MaterialModel_NeoHookean_Incompressible_Mixed.hpp"
#include "PLocAssem_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha.hpp"
#include "PLocAssem_Tet4_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha.hpp"
#include "PLocAssem_Tet4_VMS_Seg_Debug.hpp"
#include "PLocAssem_Tet4_VMS_Seg_Incompressible.hpp"
#include "PDNSolution_Mixed_U_Hyperelastic_3D.hpp"
#include "PDNSolution_P_V_Mixed_Hyperelastic_3D.hpp"

// Used for pressure mass matrix assembly
#include "PDNSolution_Disp_3D.hpp"
#include "PDNSolution_Pres_3D.hpp"
#include "PGAssem_2x2Block_HED.hpp"
#include "PLocAssem_2x2Block_VMS_HE.hpp"

#include "PGAssem_Seg_Ipt_FEM.hpp"
#include "PTime_Seg_Solver.hpp"
#include "PLinear_Solver_DiagScale_Shell.hpp"

#ifdef SLEPC  
#include "PEigen_Solver_SLEPc.hpp"
#endif

int main(int argc, char *argv[])
{
  // User defined schur pc
  bool is_usr_schur = false;

  // Number of quadrature points
  int nqp_tet = 5, nqp_tri = 4; 

  // Partition file base name
  std::string part_file("part");

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-3;
  double nl_atol = 1.0e-6;
  double nl_dtol = 10.0;
  int nl_maxits = 20;
  int nl_refreq = 4;

  // Time stepping parameters
  double initial_time = 0.0;
  double initial_step = 0.1;
  int initial_index = 0;
  double final_time = 1.0;
  std::string sol_bName("SOL_");
  int ttan_renew_freq = 1;
  int sol_record_freq = 1;

  // Restart options
  bool is_restart = false;
  int restart_index = 0;
  double restart_time = 0.0;
  double restart_step = 1.0e-3;
  std::string restart_name = "SOL_";

  PetscMPIInt rank, size;
  // ===== Initialization PETSc =====
#ifdef SLEPC  
  SlepcInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#endif

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ===== Command Line Argument =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionBool("-is_usr_schur", is_usr_schur);
  SYS_T::GetOptionInt("-nqp_tet", nqp_tet);
  SYS_T::GetOptionInt("-nqp_tri", nqp_tri);
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

  // ===== Print Command Line Argument =====
  SYS_T::cmdPrint("-part_file:", part_file);
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
  if(is_restart)
  {
    PetscPrintf(PETSC_COMM_WORLD, "-is_restart: true \n");
    SYS_T::cmdPrint("-restart_index:", restart_index);
    SYS_T::cmdPrint("-restart_time:", restart_time);
    SYS_T::cmdPrint("-restart_step:", restart_step);
    SYS_T::cmdPrint("-restart_name:", restart_name);
  }
  else PetscPrintf(PETSC_COMM_WORLD, "-is_restart: false \n");

  if(is_usr_schur) PetscPrintf(PETSC_COMM_WORLD, "-is_usr_schur: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-is_usr_schur: false \n");

  // ===== Main Data Structure =====
  SYS_T::commPrint("===> Reading mesh files ... ");
  
  FEANode * fNode = new FEANode(part_file, rank);
  
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file, rank);
  
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file ,rank); 
  
  ALocal_EBC * locebc = new ALocal_EBC_wIntPts(part_file, rank);
  
  APart_Node * pNode = new APart_Node(part_file, rank);
  SYS_T::commPrint("Done! \n");
 
  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  PetscPrintf(PETSC_COMM_WORLD, 
      "===> %d processor(s) are assigned for FEM analysis. \n", size);

  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );  
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  
  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv, * elements;
  if( GMIptr->get_elemType() == 531 )
  {
    elementv = new FEAElement_Tet4( quadv-> get_num_quadPts() );
    elements = new FEAElement_Triangle3_3D_der0( quads-> get_num_quadPts() );
  }
  else if( GMIptr->get_elemType() == 532 )
  {
    elementv = new FEAElement_Tet10( quadv-> get_num_quadPts() );
    elements = new FEAElement_Triangle6_3D_der0( quads-> get_num_quadPts() );
  }
  else
    SYS_T::print_fatal("Error: The element type is not implemented. \n");

  
  // upate matrix for strong essential bc imposition
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);
  pmat->gen_perm_bc(pNode, locnbc); 

  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");
  const double genA_spectrium = 0.5;
  const bool genA_is2ndSystem = false;
  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(
      genA_spectrium, genA_is2ndSystem);
  tm_galpha_ptr->print_info();

  SYS_T::commPrint("===> Setup the Material model.\n");
  const double mat_in_r = 1.0;
  const double mat_in_E = 2.405e8;
  const double mat_in_nu = 0.4998;
  
  IMaterialModel * matmodel = new MaterialModel_NeoHookean_ST91_Mixed(
      mat_in_r, mat_in_E, mat_in_nu );

  //IMaterialModel * matmodel = new MaterialModel_NeoHookean_Incompressible_Mixed(
  //    1.0e0 );
  //IMaterialModel * matmodel = new MaterialModel_NeoHookean_Quad_Mixed(
  //    1100, 1.7e7, 0.45 );

  IPLocAssem * locAssem_ptr 
    = new PLocAssem_VMS_Seg_Hyperelastic_3D_FEM_GenAlpha(
      matmodel, tm_galpha_ptr, GMIptr->get_nLocBas(), 
      quadv->get_num_quadPts(), elements->get_nLocBas() );  

  //IPLocAssem * locAssem_ptr = new PLocAssem_Tet4_VMS_Seg_Incompressible(
  //    matmodel, tm_galpha_ptr, GMIptr->get_nLocBas(), 
  //    quadv->get_num_quadPts(), elements->get_nLocBas() );  

  //IPLocAssem * locAssem_ptr = new PLocAssem_Tet4_VMS_Seg_Debug(
  //    matmodel, tm_galpha_ptr, GMIptr->get_nLocBas(), 
  //    quadv->get_num_quadPts(), elements->get_nLocBas() );  

  // FEA.3 Initial Condition
  // Full solution vector including U-P-V
  PDNSolution * sol = new PDNSolution_Mixed_U_Hyperelastic_3D( pNode,
      locAssem_ptr, fNode, 0 );

  PDNSolution * dot_sol = new PDNSolution_Mixed_U_Hyperelastic_3D( pNode,
      locAssem_ptr, fNode, 0 );


  if( is_restart )
  {
    initial_index = restart_index;
    initial_time  = restart_time;
    initial_step  = restart_step;

    SYS_T::file_check(restart_name.c_str());

    sol->ReadBinary(restart_name.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "===> Read sol from disk as a restart run... \n");
    PetscPrintf(PETSC_COMM_WORLD, "     restart_name: %s \n", restart_name.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "     restart_time: %e \n", restart_time);
    PetscPrintf(PETSC_COMM_WORLD, "     restart_index: %d \n", restart_index);
    PetscPrintf(PETSC_COMM_WORLD, "     restart_step: %e \n", restart_step);
  }

  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step); 
  // FEM.4 Global assembly routine
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_Seg_Ipt_FEM( locAssem_ptr,
      GMIptr, locElem, locIEN, pNode, locnbc, locebc );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr, locIEN,
      pNode, locnbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // FEM.5 Initialize the velo vector by mass matrix
  SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n");
  PLinear_Solver_PETSc * lsolver_acce = new PLinear_Solver_PETSc(
      1.0e-12, 1.0e-85, 1.0e30, 1000, "ls_mass_", "pc_mass_" );

  KSPSetType(lsolver_acce->ksp, KSPGMRES);
  KSPGMRESSetOrthogonalization(lsolver_acce->ksp,
      KSPGMRESModifiedGramSchmidtOrthogonalization);
  KSPGMRESSetRestart(lsolver_acce->ksp, 500);

  lsolver_acce->Monitor();

  PC preproc; lsolver_acce->GetPC(&preproc);
  PCSetType( preproc, PCHYPRE );

  gloAssem_ptr->Assem_mass_residual( sol, locElem, locAssem_ptr, elementv,
      elements, quadv, quads, locIEN, pNode, fNode, locnbc, locebc );

  PDNSolution * dot_pres_velo = new PDNSolution_P_V_Mixed_Hyperelastic_3D( pNode,
      locAssem_ptr, fNode, 0 );

  lsolver_acce->Solve( gloAssem_ptr->K, gloAssem_ptr->G, dot_pres_velo);

  lsolver_acce->Info();
  
  // Note: dot_sol is initialized as a zero vector. So here += equals an
  //       assignment operation. We always initialize dot_sol as a zero vector.
  // [dot_p, dot_v] = - M^-1 [Res_p, Res_v]
  SEG_SOL_T::PlusAiPV(0.0, -1.0, -1.0, dot_pres_velo, dot_sol);
  
  // dot_u = v
  SEG_SOL_T::PlusAiVPV(1.0, 0.0, 0.0, sol, dot_sol);

  delete lsolver_acce; delete dot_pres_velo;
  SYS_T::commPrint("\n===> Consistent initial acceleration is obtained.");
  SYS_T::commPrint("\n===> The mass matrix lsolver is destroyed. \n\n");


  // FEM.5.5 If we need Pressure mass matrix as a Schur preconditioner, assembly
  // it.
  Mat Mp; // Pressure mass matrix
  Mat Mu; // Velocity mass matrix
  PDNSolution * block_disp = new PDNSolution_Disp_3D(pNode, fNode, 0);
  PDNSolution * block_pres = new PDNSolution_Pres_3D(pNode, fNode, 0);
  PDNSolution * block_velo = new PDNSolution_Disp_3D(pNode, fNode, 0);

  IPLocAssem_2x2Block * block_locAssem_ptr = new PLocAssem_2x2Block_VMS_HE(
      matmodel, tm_galpha_ptr, GMIptr->get_nLocBas(),
      quadv->get_num_quadPts(), elements->get_nLocBas() );

  IPGAssem_2x2Block * block_gloAssem_ptr = new PGAssem_2x2Block_HED(
      block_locAssem_ptr, locIEN, pNode, locnbc, locebc );

  // Allocate block matrix
  block_gloAssem_ptr->Assem_nonzero_estimate( block_locAssem_ptr, locIEN, locnbc );
  block_gloAssem_ptr->Fix_nonzero_err_str();
  block_gloAssem_ptr->Clear_KG();

  // Block mass matrices -> Pressure mass matrix generated
  block_gloAssem_ptr -> Assem_mass_residual( block_disp, block_pres, block_velo, 
      block_locAssem_ptr,
      elementv, elements, quadv, quads, locIEN, fNode, locnbc, locebc );

  // Copy K_00, the pressure mass matrix to Mp
  MatConvert(block_gloAssem_ptr->K_00, MATSAME, MAT_INITIAL_MATRIX, &Mp);
  MatConvert(block_gloAssem_ptr->K_11, MATSAME, MAT_INITIAL_MATRIX, &Mu);

  // Scale with a factor 1/elastic modulus mu
  //MatScale( Mp, 1.0 / matmodel->get_elastic_mu() );

  delete block_gloAssem_ptr; delete block_locAssem_ptr;
  delete block_disp; delete block_pres; delete block_velo;


  // FEM.6 Linear solver context
  const double ll_rtol = 1.0e-5;
  const double ll_atol = 1.0e-50;
  const double ll_dtol = 1.0e20;
  const int ll_maxits = 200;
  PLinear_Solver_DiagScale * lsolver = new PLinear_Solver_DiagScale_Shell(
      ll_rtol, ll_atol, ll_dtol, ll_maxits, gloAssem_ptr, Mp, Mu );

  MatDestroy(&Mp); MatDestroy(&Mu);

  lsolver->Monitor();

  // linear solver preconditioners field split. 
  PC upc; lsolver->GetPC(&upc);

  // Enforce FIELDSPLIT type PC
  PCSetType(upc, PCFIELDSPLIT);

  // Enforce Schur PC
  PCFieldSplitSetType(upc, PC_COMPOSITE_SCHUR);

  // Schur factorization type
  PCFieldSplitSetSchurFactType(upc, PC_FIELDSPLIT_SCHUR_FACT_FULL);

  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields); // A_00 for velo
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);  // A_11 for pres

  // This defines the usr defined Schur approximation
  if(is_usr_schur)
    PCFieldSplitSetSchurPre( upc, PC_FIELDSPLIT_SCHUR_PRE_USER, Mp );


  // FEM 7. Nonlinear solver context
  PNonlinear_Seg_Solver * nsolver = new PNonlinear_Seg_Solver(
      pNode, locAssem_ptr, fNode,
      nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq);
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->print_info();

  // FEM.8 Time solver context
  PTime_Seg_Solver * tsolver = new PTime_Seg_Solver( sol_bName, 
      sol_record_freq, ttan_renew_freq, final_time );
  SYS_T::commPrint("===> Time marching solver setted up:\n");
  tsolver->print_info();

  // FEM.9 FEM analysis
  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver->TM_Seg_GenAlpha_DiagScale(is_restart, dot_sol, sol, tm_galpha_ptr,
      timeinfo, locElem, locIEN, pNode, fNode, locnbc, locebc, pmat,
      elementv, elements, quadv, quads, 
      locAssem_ptr, gloAssem_ptr, lsolver, nsolver);

  lsolver->Info();

  // ===== PETSc Finalize =====
  delete tsolver; delete nsolver; delete lsolver;
  delete gloAssem_ptr; delete timeinfo; delete locAssem_ptr;
  delete matmodel; delete tm_galpha_ptr; delete pmat;
  delete sol; delete dot_sol;
  delete elements; delete elementv; delete quads; delete quadv;
  delete pNode; delete locebc; delete locnbc; delete locElem;
  delete PartBasic; delete GMIptr; delete fNode; delete locIEN;

#ifdef SLEPC
  SlepcFinalize();
#else
  PetscFinalize();
#endif
  return 0;
}


// EOF
