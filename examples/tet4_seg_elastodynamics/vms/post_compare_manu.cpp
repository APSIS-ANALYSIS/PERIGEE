// ==================================================================
// post_compare_manu.cpp
// ------------------------------------------------------------------
// This is the postprocessing driver for computing the solution error
// by comparing with given manufactured solutions.
//
// This is a parallel routine that requires postpart_xxxxx.h5 files.
//
// Date: Feb. 9 2017
// ==================================================================
#include "HDF5_PartReader.hpp"
#include "FEANode.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Tet4.hpp"
#include "PostVectSolution.hpp"
#include "Post_error_tet4_eladyn.hpp"

int main( int argc, char * argv[] )
{
  int nqp_tet = 29; // 4 (2), 5 (3), 17 (5), or 29 (6)

  std::string sol_name("SOL_900000000");
  double sol_time = 1.0;
  std::string part_file("postpart");
  const int dof = 7;

  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  SYS_T::GetOptionInt("-nqp", nqp_tet);
  SYS_T::GetOptionReal("-sol_time", sol_time);
  SYS_T::GetOptionString("-sol_name", sol_name);
  SYS_T::GetOptionString("-part_file", part_file);

  SYS_T::cmdPrint("-nqp:", nqp_tet);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-sol_name:", sol_name);
  SYS_T::cmdPrint("-sol_time:", sol_time);

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
  double loc_p[4], loc_u[4], loc_v[4], loc_w[4];
  double loc_vu[4], loc_vv[4], loc_vw[4];
  double el2pres = 0.0, el2F = 0.0, el2disp = 0.0, el2velo = 0.0; 
  double el2sigma = 0.0;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    locIEN -> get_LIEN_e(ee, IEN_e);
    fNode -> get_ctrlPts_xyz( 4, IEN_e, ectrl_x, ectrl_y, ectrl_z );

    elem -> buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );
  
    pSolu -> get_esol(0, 4, IEN_e, loc_u);
    pSolu -> get_esol(1, 4, IEN_e, loc_v);
    pSolu -> get_esol(2, 4, IEN_e, loc_w);
    
    pSolu -> get_esol(3, 4, IEN_e, loc_p);
    
    pSolu -> get_esol(4, 4, IEN_e, loc_vu);
    pSolu -> get_esol(5, 4, IEN_e, loc_vv);
    pSolu -> get_esol(6, 4, IEN_e, loc_vw);
    
    el2pres += POST_T_TET4_ELADYN::get_manu_scalar_l2_error(loc_p, elem, ectrl_x,
        ectrl_y, ectrl_z, Int_w, R, sol_time);

    el2F += POST_T_TET4_ELADYN::get_manu_matrix_l2_error(loc_u, loc_v, loc_w,
        elem, ectrl_x, ectrl_y, ectrl_z, Int_w, R, Rx, Ry, Rz, sol_time);
    
    el2disp += POST_T_TET4_ELADYN::get_manu_disp_error(loc_u, loc_v, loc_w,
        elem, ectrl_x, ectrl_y, ectrl_z, Int_w, R, sol_time);
  
    el2velo += POST_T_TET4_ELADYN::get_manu_velo_error(loc_vu, loc_vv, loc_vw,
        elem, ectrl_x, ectrl_y, ectrl_z, Int_w, R, sol_time);
    
    el2sigma += POST_T_TET4_ELADYN::get_manu_sigma_l2_error( loc_u, loc_v, loc_w,
        elem, ectrl_x, ectrl_y, ectrl_z, Int_w, R, Rx, Ry, Rz, sol_time );
  }


  double l2disp = 0.0;
  MPI_Reduce(&el2disp, &l2disp, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2disp = sqrt(l2disp);
  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in L2 norm of disp is : %e \n", l2disp);
  
  double l2velo = 0.0;
  MPI_Reduce(&el2velo, &l2velo, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2velo = sqrt(l2velo);
  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in L2 norm of velo is : %e \n", l2velo);

  double l2p = 0.0;
  MPI_Reduce(&el2pres, &l2p, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2p = sqrt(l2p);
  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in L2 norm of pressure is : %e \n", l2p);
  
  double l2F = 0.0;
  MPI_Reduce(&el2F, &l2F, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2F = sqrt(l2F);
  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in L2 norm of F is : %e \n", l2F);

  double l2sigma = 0.0;
  MPI_Reduce(&el2sigma, &l2sigma, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2sigma = sqrt(l2sigma);
  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in L2 norm of sigma is: %e \n", l2sigma);

  delete elem; delete Int_w; delete quad;
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  delete locElem; delete pNode;
  PetscFinalize();
  return 0;
}

// EOF
