// ==================================================================
// Tet4 NS postprocessing compare with manufactured solution routine.
// Date: June 28 2017 
// ==================================================================
#include "HDF5_PartReader.hpp"
#include "FEANode.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Tet4.hpp"
#include "PostVectSolution.hpp"
#include "Post_error_tet4_ns.hpp"

int main( int argc, char * argv[] )
{
  int nqp_tet = 29; // 4 (2), 5 (3), 17 (5), or 29 (6)

  std::string sol_name("SOL_900000000");
  double sol_time = 1.0;
  std::string part_file("postpart");
  const int dof = 4;

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

  double R[4], Rx[4], Ry[4], Rz[4]; 
  int IEN_e[4];
  double ectrl_x[4], ectrl_y[4], ectrl_z[4];
  double loc_p[4], loc_u[4], loc_v[4], loc_w[4];
  double el2pres = 0.0, el2velo = 0.0, eh1pres = 0.0;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    locIEN -> get_LIEN_e(ee, IEN_e);
    fNode -> get_ctrlPts_xyz( 4, IEN_e, ectrl_x, ectrl_y, ectrl_z );

    elem -> buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

    pSolu -> get_esol(0, 4, IEN_e, loc_p);
    pSolu -> get_esol(1, 4, IEN_e, loc_u);
    pSolu -> get_esol(2, 4, IEN_e, loc_v);
    pSolu -> get_esol(3, 4, IEN_e, loc_w);

    el2pres += POST_T_TET4_NS::get_pres_l2_error(loc_p, elem, ectrl_x,
        ectrl_y, ectrl_z, Int_w, R, sol_time);
    
    eh1pres += POST_T_TET4_NS::get_pres_h1_error(loc_p, elem, ectrl_x,
        ectrl_y, ectrl_z, Int_w, R, Rx, Ry, Rz, sol_time);

    el2velo += POST_T_TET4_NS::get_velo_l2_error(loc_u, loc_v, loc_w,
        elem, ectrl_x, ectrl_y, ectrl_z, Int_w, R, sol_time);
  }

  double l2p = 0.0;
  MPI_Reduce(&el2pres, &l2p, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2p = sqrt(l2p);
  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in L2 norm of pressure is : %e \n", l2p);

  double h1p = 0.0;
  MPI_Reduce(&eh1pres, &h1p, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  h1p = sqrt(h1p);
  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in H1 semi-norm of pressure is : %e \n", h1p);

  double l2velo = 0.0;
  MPI_Reduce(&el2velo, &l2velo, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2velo = sqrt(l2velo);
  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in L2 norm of velo is : %e \n", l2velo);

  delete pSolu;  delete elem; delete Int_w; delete quad;
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  delete locElem; delete pNode;
  PetscFinalize();
  return 0;
}


// EOF
