// ==================================================================
// ns_tet4_driver.cpp
// 
// 4-node tetrahedral finite element analysis code for 3D Navier-Stokes
// equations with variational multiscale analysis.
//
// Ref: CMAME 197 (2007) 173-201
// Date: June 10 2017
// Author: Ju Liu
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_EBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "Matrix_PETSc.hpp"
#include "PLocAssem_Tet4_VMS_NS_3D_GenAlpha.hpp"
#include "PDNTimeStep.hpp"
#include "PTime_Tet4_NS_3D_Solver.hpp"

int main(int argc, char *argv[])
{
  int nqp_tet = 5, nqp_tri = 4;

  std::string part_file("part");

  // nonlinear solver parameters
  double nl_rtol = 1.0e-3;
  double nl_atol = 1.0e-6;
  double nl_dtol = 10.0;
  int nl_maxits = 20;
  int nl_refreq = 4;

  // time stepping parameters
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
  // ===== Initialization of PETSc =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ===== Command Line Argument =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

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

  // ===== Print Command Line Arguement =====
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
  
  // ===== Main Data Structure =====
  SYS_T::commPrint("===> Reading mesh files ... ");

  FEANode * fNode = new FEANode(part_file, rank);

  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);

  ALocal_Inflow_NodalBC * locinfnbc = new ALocal_Inflow_NodalBC(part_file, rank);

  ALocal_EBC * locebc = new ALocal_EBC(part_file, rank);

  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::commPrint("Done! \n");

  // ===== Do Basic checking =====
  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  // ===== Generate Quadrature rules and FEM container =====
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

  // ===== Generate Generalized-alpha method
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");
  const double genA_spectrium = 0.5;
  const bool genA_is2ndSystem = false;
  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(
      genA_spectrium, genA_is2ndSystem);
  tm_galpha_ptr->print_info();

  // ===== Local assembly initialization =====
  IPLocAssem * locAssem_ptr = new PLocAssem_Tet4_VMS_NS_3D_GenAlpha(
      tm_galpha_ptr, GMIptr->get_nLocBas(),
      quadv->get_num_quadPts(), elements->get_nLocBas() );

  // ===== Initial condition =====
  // The inlet velocity profile
  const double flow_vol_rate = 1.0;
  PDNSolution * base = new PDNSolution_Tet4_NS_3D( pNode,
      fNode, locinfnbc, flow_vol_rate, 2 );

  PDNSolution * sol = new PDNSolution_Tet4_NS_3D( pNode, fNode, 0 );

  PDNSolution * dot_sol = new PDNSolution_Tet4_NS_3D( pNode, fNode, 0 );

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
  // ===== Global Assembly Routine ====
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_v360_FEM( locAssem_ptr,
      GMIptr, locElem, locIEN, pNode, locnbc, locebc );

  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr, locIEN,
      pNode, locnbc );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();


  // ===== Initialize the dot_sol vector by solving mass matrix =====
  SYS_T::commPrint("===> Assembly mass matrix and residual vector.\n");
  PLinear_Solver_PETSc * lsolver_acce = new PLinear_Solver_PETSc(
      1.0e-14, 1.0e-85, 1.0e30, 1000, "ls_mass_", "pc_mass_" );

  KSPSetType(lsolver_acce->ksp, KSPGMRES);
  KSPGMRESSetOrthogonalization(lsolver_acce->ksp,
      KSPGMRESModifiedGramSchmidtOrthogonalization);
  KSPGMRESSetRestart(lsolver_acce->ksp, 500);

  PC preproc; lsolver_acce->GetPC(&preproc);
  PCSetType( preproc, PCASM );

  gloAssem_ptr->Assem_mass_residual( sol, locElem, locAssem_ptr, elementv,
      elements, quadv, quads, locIEN, pNode, fNode, locnbc, locebc );

  lsolver_acce->Solve( gloAssem_ptr->K, gloAssem_ptr->G, dot_sol);
  dot_sol->ScaleValue(-1.0);

  delete lsolver_acce;
  SYS_T::commPrint("\n===> Consistent initial acceleration is obtained.");
  SYS_T::commPrint(" The mass matrix lsolver is destroyed. \n\n");

  // ===== Linear & Nonlinear & Time solvers =====
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  PC upc; lsolver->GetPC(&upc);
  const PetscInt pfield[1] = {0}, vfields[] = {1,2,3};
  PCFieldSplitSetBlockSize(upc,4);
  PCFieldSplitSetFields(upc,"u",3,vfields,vfields);
  PCFieldSplitSetFields(upc,"p",1,pfield,pfield);

  PNonlinear_Tet4_NS_3D_Solver * nsolver = new PNonlinear_Tet4_NS_3D_Solver(
      pNode, fNode, nl_rtol, nl_atol, nl_dtol, nl_maxits, nl_refreq);
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->print_info();

  PTime_Tet4_NS_3D_Solver * tsolver = new PTime_Tet4_NS_3D_Solver( sol_bName,
      sol_record_freq, ttan_renew_freq, final_time );
  SYS_T::commPrint("===> Time marching solver setted up:\n");
  tsolver->print_info();

  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver->TM_VMS_GenAlpha(is_restart, base, dot_sol, sol, tm_galpha_ptr,
      timeinfo, locElem, locIEN, pNode, fNode, locnbc, locinfnbc, locebc, pmat,
      elementv, elements, quadv, quads,
      locAssem_ptr, gloAssem_ptr, lsolver, nsolver);

  // ===== PETSc Finalize =====
  delete tsolver; delete nsolver; delete lsolver;
  delete gloAssem_ptr; delete timeinfo;
  delete base;  delete sol; delete dot_sol;
  delete locAssem_ptr; delete tm_galpha_ptr;
  delete pmat; delete locinfnbc;
  delete elementv; delete elements; delete quadv; delete quads;
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic; 
  delete locElem; delete locnbc; delete locebc; delete pNode;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF 
