// ==================================================================
// post_compare_manu.cpp
// ------------------------------------------------------------------
// This is the postprocessing driver for computing the solution error
// by comparing with given manufactured solutions.
//
// This is a parallel routine that requires postpart_xxxxx.h5 files.
//
// Date: May 11 2017
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node.hpp"
#include "ALocal_Elem.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "AInt_Weight.hpp"
#include "FEAElement_Tet4.hpp"
#include "PostVectSolution.hpp"
#include "Post_error_tet4_le.hpp"

int main( int argc, char * argv[] )
{
  int nqp_tet = 29; // 4 (2), 5 (3), 17 (5), or 29 (6)

  std::string sol_name("SOL_900000000");
  std::string part_file("postpart");
  const int dof = 4;

  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  SYS_T::GetOptionInt("-nqp", nqp_tet);
  SYS_T::GetOptionString("-sol_name", sol_name);

  SYS_T::cmdPrint("-nqp:", nqp_tet);
  SYS_T::cmdPrint("-sol_name:", sol_name);

  FEANode * fNode = new FEANode(part_file, rank);
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file, rank);
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  APart_Node * pNode = new APart_Node(part_file, rank);

  PetscPrintf(PETSC_COMM_WORLD, "\n===> %d processor(s) are assigned for:", size);
  SYS_T::commPrint("Postprocessing - compute error from manufactured solutions.\n");

  IQuadPts * quad = new QuadPts_Gauss_Tet( nqp_tet );
  AInt_Weight * Int_w = new AInt_Weight( quad );
  FEAElement * elem = new FEAElement_Tet4( quad -> get_num_quadPts() );

  PostVectSolution * pSolu = new PostVectSolution( sol_name,
      "node_mapping.h5", "post_node_mapping.h5", pNode, GMIptr, dof );

  // There are four local basis functions in tet element
  double R[4]; double Rx[4]; double Ry[4]; double Rz[4];

  int IEN_e[4];
  double ectrl_x[4], ectrl_y[4], ectrl_z[4];
  double loc_u[4], loc_v[4], loc_w[4], loc_p[4];
  double el2disp = 0.0;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    locIEN -> get_LIEN_e(ee, IEN_e);
    fNode -> get_ctrlPts_xyz( 4, IEN_e, ectrl_x, ectrl_y, ectrl_z );

    elem -> buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );
  
    pSolu -> get_esol(0, 4, IEN_e, loc_p);
    pSolu -> get_esol(1, 4, IEN_e, loc_u);
    pSolu -> get_esol(2, 4, IEN_e, loc_v);
    pSolu -> get_esol(3, 4, IEN_e, loc_w);
    
    el2disp += POST_T_TET4_LE::get_manu_disp_error(loc_u, loc_v, loc_w,
        elem, ectrl_x, ectrl_y, ectrl_z, Int_w, R);
  }

  double l2disp = 0.0;
  MPI_Reduce(&el2disp, &l2disp, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2disp = sqrt(l2disp);
  PetscPrintf(PETSC_COMM_WORLD, "Error in L2 norm of disp is : %e \n", l2disp);
  
  delete elem; delete Int_w; delete quad; delete locElem; delete pNode;
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  PetscFinalize();
  return 0;
}

// EOF
