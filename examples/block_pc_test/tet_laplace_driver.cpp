// ==================================================================
// tet_laplace_driver.cpp
// 
// This is the laplace equation solver based on tets using the shell
// matrix and shell preconditioner based on the Sherman-Morrison 
// formula
//
// Author: Ju Liu
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "ALocal_Elem.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_InflowBC.hpp"
#include "ALocal_NBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "Matrix_PETSc.hpp"
#include "PLocAssem_Tet_Heat.hpp"
#include "PGAssem_Heat.hpp"
#include "PDNSolution_Heat.hpp"
#include "PLinear_Solver_PETSc.hpp"

PetscErrorCode MyMatMult(Mat shell, Vec X, Vec Y)
{
  void * ptr;
  PGAssem_Heat * user;
  MatShellGetContext(shell, &ptr);
  user = (PGAssem_Heat*) ptr;
  
  MatMult(user->K, X, Y);

  const int num_ebc = user->get_num_ebc();

  PetscScalar * dots = new PetscScalar [num_ebc];
  
  VecMDot(X, num_ebc, user->intNA, dots);  

  std::vector<double> theta{ 1.0e0, 1.0e2, 1.0e3, 2.0e3, 2.0e4, 2.0e5, 3.0e3, 3.0e4, 3.0e5, 5.0e5};

  for(int ii=0; ii<num_ebc; ++ii)
  {
    dots[ii] = dots[ii] * theta[ii];
  }
  
  VecMAXPY(Y, num_ebc, dots, user->intNA);
 
  delete [] dots; dots = nullptr;
  
  return 0;
}


PetscErrorCode MyPCSetup(PC pc)
{
  std::cout<<"SETUP \n";
  PGAssem_Heat * ptr;
  PCShellGetContext(pc,(void**)&ptr);

  KSPSetOperators(ptr->ksp_K, ptr->K, ptr->K);

  MatGetDiagonal(ptr->K, ptr->diag);

  VecReciprocal(ptr->diag);

  return 0;
}


PetscErrorCode MyPCApply(PC pc, Vec x, Vec y)
{
  std::cout<<"APPLY \n";
  PGAssem_Heat * ptr;
  PCShellGetContext(pc,(void**)&ptr);

  KSPSolve(ptr->ksp_K, x, y);
  
  const int num_ebc = 0;
  //const int num_ebc = ptr -> get_num_ebc();

  std::vector<double> theta{ 1.0e0, 1.0e2, 1.0e3, 2.0e3, 2.0e4, 2.0e5, 3.0e3, 3.0e4, 3.0e5, 5.0e5};

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    // b = diag(F)^-1 a
    //VecPointwiseMult(ptr->b, ptr->intNA[ebc_id], ptr->diag);
    KSPSolve(ptr->ksp_K, ptr->intNA[ebc_id], ptr->b);

    PetscScalar b, c, val;
    VecDot(ptr->b, x, &b);  // b = x dot b
    VecDot(ptr->b, ptr->intNA[ebc_id], &c); // c = a dot b
    
    val = -1.0 * theta[ebc_id] * b / (1.0 + theta[ebc_id] * c); // val = -1 theta b / (1 + theta c)

    VecAXPY(y, val, ptr->b); // y = y + val b
  }
  return 0;
}


int main(int argc, char *argv[])
{
  int nqp_tet = 5, nqp_tri = 4;
  int nz_estimate = 300;
  std::string part_file("part");

  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionInt("-nqp_tet", nqp_tet);
  SYS_T::GetOptionInt("-nqp_tri", nqp_tri);
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);
  SYS_T::GetOptionString("-part_file", part_file);

  SYS_T::cmdPrint("-nqp_tet:", nqp_tet);
  SYS_T::cmdPrint("-nqp_tri:", nqp_tri);
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-part_file:", part_file);

  FEANode * fNode = new FEANode(part_file, rank);
  
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);
  
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  
  ALocal_NBC * locnbc = new ALocal_NBC(part_file, rank);
  
  ALocal_InflowBC * locinfnbc = new ALocal_InflowBC(part_file, rank);
  
  ALocal_EBC * locebc = new ALocal_EBC(part_file, rank);
  
  APart_Node * pNode = new APart_Node(part_file, rank);
  
  SYS_T::commPrint("===> Data from HDF5 files are read from disk.\n");

  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis. \n", size);


  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv = nullptr;
  FEAElement * elements = nullptr;

  if( GMIptr->get_elemType() == 501 )
  {
    if( nqp_tet > 5 ) SYS_T::commPrint("Warning: the tet element is linear and you are using more than 5 quadrature points.\n");
    if( nqp_tri > 4 ) SYS_T::commPrint("Warning: the tri element is linear and you are using more than 4 quadrature points.\n");

    elementv = new FEAElement_Tet4( nqp_tet ); // elem type 501
    elements = new FEAElement_Triangle3_3D_der0( nqp_tri );
  }
  else if( GMIptr->get_elemType() == 502 )
  {
    SYS_T::print_fatal_if( nqp_tet < 29, "Error: not enough quadrature points for tets.\n" );
    SYS_T::print_fatal_if( nqp_tri < 13, "Error: not enough quadrature points for triangles.\n" );

    elementv = new FEAElement_Tet10_v2( nqp_tet ); // elem type 502
    elements = new FEAElement_Triangle6_3D_der0( nqp_tri );
  }
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);

  pmat->gen_perm_bc(pNode, locnbc);

  IPLocAssem * locAssem_ptr = new PLocAssem_Tet_Heat(
      GMIptr->get_nLocBas(), quadv->get_num_quadPts(), 
      elements->get_nLocBas(), GMIptr->get_elemType() );

  PDNSolution * sol = new PDNSolution_Heat( pNode, fNode, locinfnbc, 1 );

  PDNSolution * res = new PDNSolution( pNode, 1 );

  // Global assembly
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  PGAssem_Heat * gloAssem_ptr = new PGAssem_Heat( locAssem_ptr,
      elements, quads, GMIptr, locElem, locIEN, pNode, locnbc, locebc,
      nz_estimate );

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // Assembly matrix problem
  const double curr_time = 0.0;
  const double dt = 0.0;
  gloAssem_ptr->Assem_tangent_residual( sol, sol, sol, sol,
      curr_time, dt, locElem, locAssem_ptr, elementv, elements,
      quadv, quads, locIEN, pNode, fNode, locnbc, locebc );

  PETSc_T::MatInfo_Display_global( gloAssem_ptr->K );

  gloAssem_ptr->Assem_intNA(locAssem_ptr, elements, quads, locnbc, locebc);

  // Linear solver context
  //PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();

  //lsolver -> Solve( gloAssem_ptr->K, gloAssem_ptr->G, res );

  //sol -> PlusAX( res, -1.0 );

  //sol -> WriteBinary("SOL_900000000");

  //lsolver -> print_info();

  // Test shell matrix here
  Mat A;
  PetscInt local_row_size, local_col_size;
  MatGetLocalSize( gloAssem_ptr->K, &local_row_size, &local_col_size );

  MatCreateShell( PETSC_COMM_WORLD, local_row_size, local_col_size, PETSC_DETERMINE,
      PETSC_DETERMINE, (void *)gloAssem_ptr, &A );

  MatShellSetOperation(A, MATOP_MULT, (void(*)(void))MyMatMult);

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetTolerances(ksp, 1.0e-5, 1.0e-10, 1.0e50, 100);
  KSPSetFromOptions(ksp);

  KSPSetOperators(ksp, A, gloAssem_ptr->K);
 
  // Now set shell preconditioner
  PC pc;
  KSPGetPC(ksp,&pc);

  PCSetType(pc,PCSHELL);
  PCShellSetContext(pc, (void *)gloAssem_ptr);
  PCShellSetSetUp(pc, MyPCSetup);
  PCShellSetApply(pc, MyPCApply);
  PCShellSetName(pc, "MyPC");

  // Solve linear system 
  KSPSolve(ksp, gloAssem_ptr->G, res->solution);
  
  // Update the solution
  sol -> PlusAX( res, -1.0 );
  sol -> WriteBinary("SOL_900000000");

  KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
  KSPView(gloAssem_ptr->ksp_K, PETSC_VIEWER_STDOUT_WORLD);

  MatDestroy(&A);
  KSPDestroy(&ksp);

  // Write matrix to disk
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"mat",FILE_MODE_WRITE,&viewer);
  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  MatView(gloAssem_ptr->K,viewer);
  PetscViewerPopFormat(viewer);
  PetscViewerDestroy(&viewer);

  // Finalize
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  delete locElem; delete locnbc; delete locebc; delete pNode; 
  delete locinfnbc; delete elementv; delete elements;
  delete quads; delete quadv; delete pmat; delete locAssem_ptr;
  delete sol; delete res; delete gloAssem_ptr;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
