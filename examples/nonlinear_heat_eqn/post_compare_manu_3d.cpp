// ==================================================================
// post_compare_manu.cpp
// ------------------------------------------------------------------
// This is a postprocess driver for computing error of solution by 
// comparing with manufactured solutions.
//
// This is a parallel routine. Users have to run Pre_postprocess to 
// obtain a parallel mesh partition.
//
// Date: Dec. 12 2013
// ==================================================================
#include <cmath>
#include "Sys_Tools.hpp"
#include "QuadPts_Gauss.hpp"
#include "HDF5_PartReader.hpp"
#include "BernsteinBasis_Array.hpp"
#include "FEANode.hpp"
#include "AExtractor_3D_NURBS_xyz.hpp"
#include "FEAElement_NURBS_3D_der0_v3.hpp"
#include "FEAElement_NURBS_3D_der1_v3.hpp"
#include "APart_Node.hpp"
#include "AGlobal_Mesh_Info_1Patch_NURBS_3D.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "ALocal_meshSize_3D_NURBS.hpp"
#include "APart_Basic_Info.hpp"
#include "AInt_Weight.hpp"
#include "PostVectSolution.hpp"
#include "Post_error.hpp"
using namespace std;

int main( int argc, char * argv[] )
{
  // Number of quadrature points
  int nqpx = 3; int nqpy = 3; int nqpz = 3;

  // Solution name
  string sol_name("SOL_900000000");

  // Solution time
  double sol_time = 1.0;

  // partition file base name
  string part_file("postpart");

  // degree of freedom of the problem
  int dof = 1;
  int field = 0;

  // Mesh type
  int mesh_type = 0;

  PetscMPIInt rank, size;
  // ====== PETSc Initialize =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ====== Read command line =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  
  SYS_T::GetOptionInt("-nqpx", nqpx);
  SYS_T::GetOptionInt("-nqpy", nqpy);
  SYS_T::GetOptionInt("-nqpz", nqpz);
  SYS_T::GetOptionInt("-mesh_type", mesh_type);
  SYS_T::GetOptionInt("-field", field);
  SYS_T::GetOptionReal("-sol_time", sol_time);
  SYS_T::GetOptionString("-sol_name", sol_name);
  SYS_T::GetOptionString("-part_file", part_file);

  SYS_T::cmdPrint("-nqpx:", nqpx); 
  SYS_T::cmdPrint("-nqpy:",nqpy);
  SYS_T::cmdPrint("-nqpz:", nqpz);
  SYS_T::cmdPrint("-field:", field);
  SYS_T::cmdPrint("-mesh_type:", mesh_type);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-sol_name:", sol_name);
  SYS_T::cmdPrint("-sol_time:", sol_time);

  // ===== Read Partition file =====
  HDF5_PartReader * h5reader = new HDF5_PartReader(part_file, rank);
  SYS_T::commPrint("===> Reading mesh files ... \n");
  // get global mesh info
  IAGlobal_Mesh_Info * gInfo_ptr = new AGlobal_Mesh_Info_1Patch_NURBS_3D(h5reader);

  // get node partition info
  APart_Node * pNode = new APart_Node(part_file, rank);  

  // get degree-of-freedom from partition files
  //h5reader->get_GMI_dofNum( dof );
  //assert(field < dof);

  // get control points
  FEANode * fNode = new FEANode(part_file, rank);

  // get mesh size
  IALocal_meshSize * locmSize = new ALocal_meshSize_3D_NURBS(h5reader);

  // get extraction operator
  IAExtractor * fExt = new AExtractor_3D_NURBS_xyz(h5reader);

  // get LIEN
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // get local element info
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // get partition info
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  SYS_T::synPrint(" done\t", rank);
  delete h5reader;

  if(size != PartBasic->get_cpu_size())
    MPI_Abort(PETSC_COMM_WORLD, 1);

  PetscPrintf(PETSC_COMM_WORLD, "\n===> %d processor(s) are assigned for:", size);
  PetscPrintf(PETSC_COMM_WORLD, "Postprocessing - compute error heat manufactured solutions.\n");

  // build finite elements -----------------------------------
  SYS_T::commPrint("\n===> Build quadrature rules ... \n");
  IQuadPts * quad_z = new QuadPts_Gauss(nqpz);
  IQuadPts * quad_y = new QuadPts_Gauss(nqpy);
  IQuadPts * quad_x = new QuadPts_Gauss(nqpx);

  SYS_T::commPrint("===> Build quadrature weight ... \n");
  AInt_Weight * Int_w = new AInt_Weight(quad_x, quad_y, quad_z);

  SYS_T::commPrint("===> Build univariate Bezier elements ... \n");
  BernsteinBasis_Array Bena_x(gInfo_ptr->get_xdegree(), quad_x);
  BernsteinBasis_Array Bena_y(gInfo_ptr->get_ydegree(), quad_y);
  BernsteinBasis_Array Bena_z(gInfo_ptr->get_zdegree(), quad_z);
  // -----------------------------------------------------------

  // ===== Manage Solution for postprocessing
  PostVectSolution * pSolu = new PostVectSolution( sol_name,
      "node_mapping.h5", "post_node_mapping.h5", pNode, gInfo_ptr, dof );

  // ===== Build element and calculate error
  int nLocBas = gInfo_ptr->get_nLocBas();

  int * IEN_e = new int [nLocBas];
  double * R = new double [nLocBas];
  double * Rx = new double [nLocBas];
  double * Ry = new double [nLocBas];
  double * Rz = new double [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  double * loc_sol = new double [nLocBas];

  double el2error = 0.0;
  double eh1error = 0.0;
  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    FEAElement * elem_ptr = new FEAElement_NURBS_3D_der1_v3(ee, locmSize,
        &Bena_x, &Bena_y, &Bena_z, fNode, fExt, locIEN);

    locIEN->get_LIEN_e(ee, IEN_e);
    fNode->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    pSolu->get_esol(0, nLocBas, IEN_e, loc_sol);

    el2error += POST_T::get_manu_scalar_l2_error(loc_sol, elem_ptr, ectrl_x,
        ectrl_y, ectrl_z, Int_w, R, sol_time);

    eh1error += POST_T::get_manu_scalar_h1_error(loc_sol, elem_ptr, ectrl_x,
        ectrl_y, ectrl_z, Int_w, R, Rx, Ry, Rz, sol_time);

    delete elem_ptr;
  }

  double l2error = 0.0;
  double h1error = 0.0;
  PetscBarrier(NULL);

  MPI_Reduce(&el2error, &l2error, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&eh1error, &h1error, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);

  l2error = sqrt(l2error); 
  h1error = sqrt(h1error); 

  PetscPrintf(PETSC_COMM_WORLD, "Error in L2 norm is : %e \n", l2error);
  PetscPrintf(PETSC_COMM_WORLD, "Error in H1 norm is : %e \n", h1error);
  // PetscBarrier(NULL);
  //PetscPrintf(PETSC_COMM_SELF, " --- rank %d norm is : %e \n", rank, el2error);


  delete [] IEN_e;
  delete [] R;
  delete [] ectrl_x;
  delete [] ectrl_y;
  delete [] ectrl_z;
  delete [] loc_sol;

  // ===== PETSc Finalize =====
  delete pNode;
  delete pSolu;
  delete gInfo_ptr;
  delete fNode;
  delete locmSize;
  delete fExt;
  delete locIEN;
  delete locElem;
  delete PartBasic;

  delete quad_x; delete quad_y; delete quad_z;

  PetscFinalize();
  return 0;
}

// EOF
