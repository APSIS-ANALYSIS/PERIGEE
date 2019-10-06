// ==================================================================
// le3d_driver.cpp
// ------------------------------------------------------------------
// Three-dimensional driver for the 3D linear elasticity.
//
// Date: May 9th 2017
// Author: Ju Liu
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "Matrix_PETSc.hpp"
#include "PLocAssem_LinearElastic_3D_Static.hpp"
#include "PDNSolution_LinearElastic_3D.hpp"
#include "PGAssem_v360_FEM.hpp"
#include "PTime_Solver.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature points
  int nqp_tet = 29, nqp_tri = 13;

  // Partition file base name
  std::string part_file("part");

  // Linear solver parameters
  double ll_rtol = 1.0e-5;
  double ll_atol = 1.0e-55;
  double ll_dtol = 1.0e30;
  int ll_maxits = 1000;

  std::string sol_bName("SOL_");

  PetscMPIInt rank, size;
  // ===== Initialization PETSc =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ===== Command Line Argument =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-nqp_tet", nqp_tet);
  SYS_T::GetOptionInt("-nqp_tri", nqp_tri);
  SYS_T::GetOptionString("-sol_name", sol_bName);

  // ===== Print Command Line Argument =====
  SYS_T::cmdPrint("-nqp_tet:", nqp_tet); 
  SYS_T::cmdPrint("-nqp_tri:", nqp_tri);
  SYS_T::cmdPrint("-sol_name:", sol_bName);

  // ===== Main Data Structure =====
  SYS_T::commPrint("===> Reading mesh files ... \n");

  FEANode * fNode = new FEANode(part_file, rank);
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file, rank);
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);
  ALocal_EBC * locebc = new ALocal_EBC(part_file, rank);
  APart_Node * pNode = new APart_Node(part_file, rank);
  SYS_T::commPrint("Done! \n");

  if( size != PartBasic->get_cpu_size() ) SYS_T::print_fatal("Error: Assigned CPU number does not match the partition. \n");

  PetscPrintf(PETSC_COMM_WORLD, "\n===> %d processor(s) are assigned for FEM analysis. \n", size);

  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv = new FEAElement_Tet4( quadv-> get_num_quadPts() );
  FEAElement * elements = new FEAElement_Triangle3_3D_der0(
      quads-> get_num_quadPts() );

  Matrix_PETSc * pmat = new Matrix_PETSc(pNode);
  pmat->gen_perm_bc(pNode, locnbc);

  const double mat_in_E = 1.0;
  const double mat_in_nu = 0.3;
  IPLocAssem * locAssem_ptr = new PLocAssem_LinearElastic_3D_Static(
      mat_in_E, mat_in_nu, GMIptr->get_nLocBas(),
      quadv->get_num_quadPts(), elements->get_nLocBas() );

  // This is the initial 'guess' that has the desired dirichlet bc value
  PDNSolution * disp = new PDNSolution_LinearElastic_3D(
      pNode, locAssem_ptr, fNode, -1);

  PDNSolution * update = new PDNSolution_LinearElastic_3D(
      pNode, locAssem_ptr, fNode, 0 );


  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_v360_FEM( locAssem_ptr,
      GMIptr, locElem, locIEN, pNode, locnbc, locebc );

  gloAssem_ptr->Assem_tangent_residual(
      disp, disp, 0.0, 0.0, locElem, locAssem_ptr, elementv, elements,
      quadv, quads, locIEN, pNode, fNode, locnbc, locebc );

  gloAssem_ptr->Fix_nonzero_err_str();

  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc(
      ll_rtol, ll_atol, ll_dtol, ll_maxits );

  lsolver->Solve( gloAssem_ptr->K, gloAssem_ptr->G, update );


  disp->PlusAX( *update, -1.0 );
  
  disp->WriteBinary("SOL_900000000");

  SYS_T::commPrint("\n===> Static solution obtained. \n");

  // ====== PETSc Finalize =====
  delete lsolver;
  delete gloAssem_ptr; delete update; delete disp; delete locAssem_ptr;
  delete pmat; delete quadv; delete quads; delete elementv; delete elements;
  delete fNode; delete locIEN; delete GMIptr;
  delete PartBasic; delete locElem; delete pNode;
  delete locnbc; delete locebc;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
