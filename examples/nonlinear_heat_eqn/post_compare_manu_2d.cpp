// ==================================================================
// post_compare_manu_2d.cpp
// ------------------------------------------------------------------
// Thie postprocess driver computes the error of numerical solutions
// with exact solutions.
//
// It is parallel, which requires the user to run pre_postprocess2d
// first to get the mesh partition for postprocessor.
//
// Date: April 18 2014
// ==================================================================
#include <cmath>
#include <vector>

#include "Sys_Tools.hpp"
#include "IQuadPts.hpp"
#include "QuadPts_Gauss.hpp"
#include "HDF5_PartReader.hpp"

#include "BernsteinBasis_Array.hpp"
#include "FEANode.hpp"
#include "IAExtractor.hpp"
#include "AExtractor_2D_NURBS_xy.hpp"
#include "FEAElement.hpp"
#include "FEAElement_NURBS_2D_der1.hpp"
#include "FEAElement_NURBS_2D_der2.hpp"
#include "APart_Node.hpp"
#include "IAGlobal_Mesh_Info.hpp"
#include "AGlobal_Mesh_Info_1Patch_NURBS_2D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "IALocal_meshSize.hpp"
#include "ALocal_meshSize_2D_NURBS.hpp"
#include "AInt_Weight.hpp"
#include "PostVectSolution.hpp"
#include "Post_error.hpp"

using namespace std;

int main( int argc, char * argv[] )
{
  // quadrature pts
  int nqpx = 3; int nqpy = 3;

  // solution name
  string sol_name("SOL_900000000");

  // solution time
  double sol_time = 1.0;

  // partition file base name
  string part_file("postpart");

  // degree of freedom of the problem
  int dof = 1;
  int field = 0;

  PetscMPIInt rank, size;

  // ======= PETSc Initialize =======
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ======= Read command line arguments =======
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  
  SYS_T::GetOptionInt("-nqpx", nqpx);
  SYS_T::GetOptionInt("-nqpy", nqpy);
  SYS_T::GetOptionInt("-field", field);
  SYS_T::GetOptionReal("-sol_time", sol_time);
  SYS_T::GetOptionString("-sol_name", sol_name);
  SYS_T::GetOptionString("-part_file", part_file);

  SYS_T::cmdPrint("-nqpx:", nqpx);
  SYS_T::cmdPrint("-nqpy:",nqpy);
  SYS_T::cmdPrint("-field:", field);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-sol_name:", sol_name);
  SYS_T::cmdPrint("-sol_time:", sol_time);

  // ======= Read Partition file =======
  HDF5_PartReader * h5reader = new HDF5_PartReader(part_file, rank);
  SYS_T::commPrint("===> Reading mesh files ... \n");

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_1Patch_NURBS_2D(h5reader);
  APart_Node * pNode = new APart_Node(part_file, rank);
  FEANode * fNode = new FEANode(part_file, rank);
  IALocal_meshSize * locmSize = new ALocal_meshSize_2D_NURBS(h5reader);
  IAExtractor * fExt = new AExtractor_2D_NURBS_xy(h5reader);
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  SYS_T::synPrint(" done\t", rank);

  if(size != PartBasic->get_cpu_size())
    MPI_Abort(PETSC_COMM_WORLD, 1);

  int dof_check;
  h5reader->get_GMI_dofNum(dof_check);
  if(dof_check != dof)
  {
    PetscPrintf(PETSC_COMM_WORLD, 
        "Error: dof in part file does not match with the given # of dof. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  if(field >= dof)
  {
    PetscPrintf(PETSC_COMM_WORLD, "Error: given field id is out of the range {0, ..., dof-1}.\n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }


  PetscPrintf(PETSC_COMM_WORLD, "\n===> %d processor(s) are assigned for:", size);
  PetscPrintf(PETSC_COMM_WORLD, "Postprocessing - compute error 2d heat manufactured solutions.\n");

  delete h5reader;

  // ======= build finite elements =======
  SYS_T::commPrint("===> Build quadrature rules ... \n");
  IQuadPts * quad_x = new QuadPts_Gauss(nqpx);
  IQuadPts * quad_y = new QuadPts_Gauss(nqpy);

  SYS_T::commPrint("===> Build quadrature weights ... \n");
  AInt_Weight * Int_w = new AInt_Weight(quad_x, quad_y);

  SYS_T::commPrint("===> Build univariate Bezier elements ... \n");
  BernsteinBasis_Array Bena_x(GMIptr->get_xdegree(), quad_x);
  BernsteinBasis_Array Bena_y(GMIptr->get_ydegree(), quad_y);

  // ======= Manage solution for postprocessing =======
  PostVectSolution * pSolu = new PostVectSolution( sol_name,
      "node_mapping.h5", "post_node_mapping.h5", pNode, GMIptr, dof );


  // ======= Build element and calculate error ========
  int nLocBas = GMIptr->get_nLocBas();

  int * IEN_e = new int [nLocBas];
  double * R = new double [nLocBas];
  double * Rx = new double [nLocBas];
  double * Ry = new double [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  double * loc_sol = new double [nLocBas];

  double el2error = 0.0;
  double eh1error = 0.0;
  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    FEAElement * elem_ptr = new FEAElement_NURBS_2D_der1(ee, locmSize,
        &Bena_x, &Bena_y, fNode, fExt, locIEN);

    locIEN->get_LIEN_e(ee, IEN_e);
    fNode->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    pSolu->get_esol(0, nLocBas, IEN_e, loc_sol);

    el2error += POST_T::get_manu_scalar_l2_error(loc_sol, elem_ptr, ectrl_x,
        ectrl_y, Int_w, R, sol_time);

    eh1error += POST_T::get_manu_scalar_h1_error(loc_sol, elem_ptr, ectrl_x,
        ectrl_y, Int_w, R, Rx, Ry, sol_time);

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

  delete [] IEN_e;
  delete [] R; delete [] Rx; delete [] Ry;
  delete [] ectrl_x;
  delete [] ectrl_y;
  delete [] ectrl_z;
  delete [] loc_sol;

  delete pNode;
  delete pSolu;
  delete GMIptr;
  delete fNode;
  delete locmSize;
  delete fExt;
  delete locIEN;
  delete locElem;
  delete PartBasic;

  delete quad_x; delete quad_y;

  PetscFinalize();
  return 0;
}


// EOF
