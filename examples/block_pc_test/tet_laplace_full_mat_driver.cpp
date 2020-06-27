// ==================================================================
// tet_laplace_full_mat_driver.cpp
// 
// This is the laplace equation solver based on tets using the FULL
// matrix for comparison purposes.
//
// Author: Ju Liu
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "ALocal_Elem.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "ALocal_NodalBC.hpp"
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
  
  ALocal_NodalBC * locnbc = new ALocal_NodalBC(part_file, rank);
  
  ALocal_Inflow_NodalBC * locinfnbc = new ALocal_Inflow_NodalBC(part_file, rank);
  
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
  //gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  gloAssem_ptr->Release_nonzero_err_str();

  // Assembly matrix problem
  const double curr_time = 0.0;
  const double dt = 0.0;
  gloAssem_ptr->Assem_tangent_residual( sol, sol, sol, sol,
      curr_time, dt, locElem, locAssem_ptr, elementv, elements,
      quadv, quads, locIEN, pNode, fNode, locnbc, locebc );

  gloAssem_ptr->Assem_intNA(locAssem_ptr, elements, quads, locnbc, locebc);

  const int num_ebc = gloAssem_ptr->get_num_ebc();

  PetscInt nlocal;
  VecGetLocalSize(gloAssem_ptr->G, &nlocal);

  std::vector<double> val;
  std::vector<int> idx; 
  for(int ii=0; ii<num_ebc; ++ii)
  {
    PetscScalar * array;

    val.clear(); idx.clear();
    
    VecGetArray(gloAssem_ptr->intNA[ii], &array);

    for(int jj=0; jj<nlocal; ++jj)
    {
      if(array[jj] != 0.0)
      {
        val.push_back( array[jj] );
        idx.push_back( jj );
      }
    }

    VecRestoreArray(gloAssem_ptr->intNA[ii], &array);
  
    for(unsigned int row=0; row<val.size(); ++row)
    {
      for(unsigned int col=0; col<val.size(); ++col)
        MatSetValue(gloAssem_ptr->K, idx[row], idx[col], val[row]*val[col], ADD_VALUES);
    }
  }

  MatAssemblyBegin(gloAssem_ptr->K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(gloAssem_ptr->K, MAT_FINAL_ASSEMBLY);

  PETSc_T::MatInfo_Display_global( gloAssem_ptr->K );

  // Solve linear system 
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetTolerances(ksp, 1.0e-5, 1.0e-10, 1.0e50, 100);
  KSPSetFromOptions(ksp);

  KSPSetOperators(ksp, gloAssem_ptr->K, gloAssem_ptr->K);
  KSPSolve(ksp, gloAssem_ptr->G, res->solution);

  // Update the solution
  sol -> PlusAX( res, -1.0 );
  sol -> WriteBinary("SOL_900000000");

  KSPDestroy(&ksp);

  // write matrix on disk
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
