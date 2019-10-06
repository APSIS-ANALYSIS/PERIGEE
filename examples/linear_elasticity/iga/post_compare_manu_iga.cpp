// ==================================================================
// post_compare_manu_iga.cpp
// ------------------------------------------------------------------
// This is the postprocessing driver for computing the solution error
// by comparing with given manufactured solutions.
//
// This is a parallel routine that requires postpart_xxxxx.h5 files.
//
// Date: May 18 2017
// ==================================================================
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
#include "Post_error_tet4_le.hpp"

int main( int argc, char * argv[] )
{
  int nqpx = 3; int nqpy = 3; int nqpz = 3;
  std::string sol_name("SOL_900000000");
  std::string part_file("postpart");
  const int dof = 3;

  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
 
  SYS_T::GetOptionInt("-nqpx", nqpx);
  SYS_T::GetOptionInt("-nqpy", nqpy);
  SYS_T::GetOptionInt("-nqpz", nqpz); 
  SYS_T::GetOptionString("-sol_name", sol_name);
  SYS_T::GetOptionString("-part_file", part_file);

  SYS_T::cmdPrint("-nqpx:", nqpx);
  SYS_T::cmdPrint("-nqpy:",nqpy);
  SYS_T::cmdPrint("-nqpz:", nqpz);
  SYS_T::cmdPrint("-sol_name:", sol_name);
  SYS_T::cmdPrint("-part_file:", part_file);

  HDF5_PartReader * h5reader = new HDF5_PartReader( part_file, rank );

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_1Patch_NURBS_3D(h5reader);

  APart_Node * pNode = new APart_Node(part_file, rank);

  FEANode * fNode = new FEANode(part_file, rank);

  IALocal_meshSize * locmSize = new ALocal_meshSize_3D_NURBS(h5reader);

  IAExtractor * fExt = new AExtractor_3D_NURBS_xyz(h5reader);

  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  delete h5reader;
  PetscPrintf(PETSC_COMM_WORLD, "\n===> %d processor(s) are assigned for:", size);
  SYS_T::commPrint("Postprocessing - compute error from manufactured solutions.\n");

  SYS_T::commPrint("\n===> Build quadrature rules ... \n");
  IQuadPts * quad_z = new QuadPts_Gauss(nqpz);
  IQuadPts * quad_y = new QuadPts_Gauss(nqpy);
  IQuadPts * quad_x = new QuadPts_Gauss(nqpx);

  SYS_T::commPrint("===> Build quadrature weight ... \n");
  AInt_Weight * Int_w = new AInt_Weight(quad_x, quad_y, quad_z);

  SYS_T::commPrint("===> Build univariate Bezier elements ... \n");
  BernsteinBasis_Array Bena_x(GMIptr->get_xdegree(), quad_x);
  BernsteinBasis_Array Bena_y(GMIptr->get_ydegree(), quad_y);
  BernsteinBasis_Array Bena_z(GMIptr->get_zdegree(), quad_z);

  PostVectSolution * pSolu = new PostVectSolution( sol_name,
      "node_mapping.h5", "post_node_mapping.h5", pNode, GMIptr, dof );
  
  int nLocBas = GMIptr->get_nLocBas();

  int * IEN_e = new int [nLocBas];
  double * R = new double [nLocBas];
  double * Rx = new double [nLocBas];
  double * Ry = new double [nLocBas];
  double * Rz = new double [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  double * loc_u = new double [nLocBas];
  double * loc_v = new double [nLocBas];
  double * loc_w = new double [nLocBas];

  double el2disp = 0.0, el2sigma = 0.0;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    FEAElement * elem_ptr = new FEAElement_NURBS_3D_der1_v3(ee, locmSize,
        &Bena_x, &Bena_y, &Bena_z, fNode, fExt, locIEN);

    if(elem_ptr->is_sizeNonzero())
    {
      locIEN -> get_LIEN_e(ee, IEN_e);
      fNode -> get_ctrlPts_xyz( nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z );

      pSolu -> get_esol(0, nLocBas, IEN_e, loc_u);
      pSolu -> get_esol(1, nLocBas, IEN_e, loc_v);
      pSolu -> get_esol(2, nLocBas, IEN_e, loc_w);

      el2disp += POST_T_TET4_LE::get_manu_disp_error(loc_u, loc_v, loc_w,
          elem_ptr, ectrl_x, ectrl_y, ectrl_z, Int_w, R);

      el2sigma += POST_T_TET4_LE::get_manu_stress_error( loc_u, loc_v, loc_w,
          elem_ptr, ectrl_x, ectrl_y, ectrl_z, Int_w, R, Rx, Ry, Rz );
    }

    delete elem_ptr;
  }

  double l2disp = 0.0;
  MPI_Reduce(&el2disp, &l2disp, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2disp = sqrt(l2disp);
  l2disp /= 0.6124;
  PetscPrintf(PETSC_COMM_WORLD, "Error in L2 norm of disp is : %e \n", l2disp);

  double l2sigma = 0.0;
  MPI_Reduce(&el2sigma, &l2sigma, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2sigma = sqrt(l2sigma);
  l2sigma /= 5.1575;
  PetscPrintf(PETSC_COMM_WORLD, "Error in L2 norm of sigma is: %e \n", l2sigma);

  delete Int_w; 
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  delete locElem; delete pNode;
  PetscFinalize();
  return 0;
}


// EOF
