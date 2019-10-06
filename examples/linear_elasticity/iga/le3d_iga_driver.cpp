// ==================================================================
// le3d_iga_driver.cpp
// ------------------------------------------------------------------
// Three-dimensional driver for the 3D linear elasticity using IGA.
//
// Date: May 18 2017
// Author: Ju Liu
// ==================================================================
#include "AExtractor_3D_NURBS_xyz.hpp"
#include "ALocal_meshSize_3D_NURBS.hpp"
#include "AGlobal_Mesh_Info_1Patch_NURBS_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "QuadPts_Gauss.hpp"
#include "FEAElement_NURBS_3D_der1.hpp"
#include "PLocAssem_Linear_Elastostatics_3D_IGA.hpp"
#include "PGAssem_v360_NURBS.hpp"
#include "PDNSolution_LinearElastic_3D.hpp"
#include "PTime_Solver.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature points
  int nqpx = 3, nqpy = 3, nqpz = 3;

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
  
  SYS_T::GetOptionInt("-nqpx", nqpx);
  SYS_T::GetOptionInt("-nqpy", nqpy);
  SYS_T::GetOptionInt("-nqpz", nqpz);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionString("-sol_name", sol_bName);

  // ===== Print Command Line Argument =====
  SYS_T::cmdPrint("-nqpx:", nqpx); 
  SYS_T::cmdPrint("-nqpy:", nqpy);
  SYS_T::cmdPrint("-nqpz:", nqpz);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-sol_name:", sol_bName);

  // ===== Main Data Structure =====
  SYS_T::commPrint("===> Reading mesh files ... \n");
  HDF5_PartReader * h5reader = new HDF5_PartReader(part_file, rank);

  FEANode * fNode = new FEANode(part_file, rank);
  IALocal_meshSize * locmSize = new ALocal_meshSize_3D_NURBS(h5reader);
  IAExtractor * fExt = new AExtractor_3D_NURBS_xyz(h5reader);
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_1Patch_NURBS_3D(h5reader);
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);
  ALocal_ElemBC * locebc = new ALocal_ElemBC(h5reader);
  APart_Node * pNode = new APart_Node(part_file, rank);
  SYS_T::commPrint("Done! \n"); delete h5reader;

  if( size != PartBasic->get_cpu_size() ) SYS_T::print_fatal("Error: Assigned CPU number does not match the partition. \n");

  PetscPrintf(PETSC_COMM_WORLD, "\n===> %d processor(s) are assigned for FEM analysis. \n", size);

  SYS_T::commPrint("===> Build quadrature rules ... \n");
  IQuadPts * quad_x = new QuadPts_Gauss(nqpx);
  IQuadPts * quad_y = new QuadPts_Gauss(nqpy);
  IQuadPts * quad_z = new QuadPts_Gauss(nqpz);

  SYS_T::commPrint("===> Build quadrature weights ... \n");
  AInt_Weight * Int_w = new AInt_Weight(quad_x, quad_y, quad_z);

  SYS_T::commPrint("===> Build univariate Bezier elements ... \n");
  BernsteinBasis_Array Bena_x(GMIptr->get_xdegree(), quad_x);
  BernsteinBasis_Array Bena_y(GMIptr->get_ydegree(), quad_y);
  BernsteinBasis_Array Bena_z(GMIptr->get_zdegree(), quad_z);

  FEAElement * element = new FEAElement_NURBS_3D_der1(
      GMIptr->get_xdegree(), GMIptr->get_ydegree(),
      GMIptr->get_zdegree(), nqpx, nqpy, nqpz );

  Matrix_PETSc * pmat = new Matrix_PETSc(pNode);
  pmat->gen_perm_bc(pNode, locnbc);

  const double mat_in_E = 1.0;
  const double mat_in_nu = 0.3;
  IPLocAssem * locAssem_ptr = new PLocAssem_Linear_Elastostatics_3D_IGA(
      mat_in_E, mat_in_nu, GMIptr->get_nLocBas(), Int_w->get_num(),
      GMIptr->get_xdegree(), GMIptr->get_ydegree(), GMIptr->get_zdegree(),
      nqpx, nqpy, nqpz );

  // This is the initial 'guess' that has the desired dirichlet bc value
  PDNSolution * disp = new PDNSolution_LinearElastic_3D(
      pNode, locAssem_ptr, fNode, 0);

  PDNSolution * update = new PDNSolution_LinearElastic_3D(
      pNode, locAssem_ptr, fNode, 0 );


  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_v360_NURBS( locAssem_ptr,
      GMIptr, locElem, locIEN, pNode, locnbc );

  gloAssem_ptr->Assem_tangent_residual(
      disp, disp, 0.0, 0.0, locElem, locAssem_ptr, element,
      locIEN, pNode, fNode, Int_w, locmSize,
      &Bena_x, &Bena_y, &Bena_z, fExt, locnbc, locebc );

  gloAssem_ptr->Fix_nonzero_err_str();

  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc(
      ll_rtol, ll_atol, ll_dtol, ll_maxits );

  lsolver->Solve( gloAssem_ptr->K, gloAssem_ptr->G, update );


  disp->PlusAX( *update, -1.0 );

  disp->WriteBinary("SOL_900000000");

  SYS_T::commPrint("\n===> Static solution obtained. \n");

  // ====== PETSc Finalize =====
  delete lsolver;
  delete gloAssem_ptr; delete disp; delete update; 
  delete locAssem_ptr; delete element; delete pmat;
  delete quad_x; delete quad_y; delete quad_z; delete Int_w;
  delete fNode; delete locmSize; delete fExt; delete locIEN; delete GMIptr;
  delete PartBasic; delete locElem; delete pNode;
  delete locnbc; delete locebc;
  PetscFinalize();
  return 0;
}


// EOF
