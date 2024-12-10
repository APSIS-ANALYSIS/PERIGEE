// ==================================================================
// post_compare_manu.cpp
//
// Error analysis code for 3D elastodynamics.
//
// Data: Oct. 25 2023
// ==================================================================
#include "AGlobal_Mesh_Info.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "APart_Node.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "PostVectSolution.hpp"
#include "Post_error_elastodynamics.hpp"

int main( int argc, char * argv[] )
{
  int nqp_tet = 29; // 4 (2), 5 (3), 17 (5), or 29 (6)
  int nqp_hex_1D = 4;

  double sol_time = 1.0;
  std::string sol_name("SOL_900000000");
  std::string part_file("./ppart/part");
  const int dof = 3;

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  SYS_T::GetOptionString("-sol_name", sol_name);
  SYS_T::GetOptionReal("-sol_time",   sol_time);

  SYS_T::cmdPrint("-sol_name:", sol_name);
  SYS_T::cmdPrint("-sol_time:", sol_time);

  SYS_T::file_check(sol_name.c_str());

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  FEANode * fNode = new FEANode(part_file, rank);

  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  AGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info(part_file,rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::print_fatal_if( size != PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for:", size);
  SYS_T::commPrint("Postprocessing - compute error from manufactured solutions.\n");

  IQuadPts * quadv = nullptr;
  FEAElement * elementv = nullptr;
  const FEType elemType = GMIptr -> get_elemType();

  if( elemType == FEType::Tet4 )
  {
    quadv = new QuadPts_Gauss_Tet( nqp_tet );
    elementv = new FEAElement_Tet4( nqp_tet );
  }
  else if( elemType == FEType::Tet10 )
  {
    quadv = new QuadPts_Gauss_Tet( nqp_tet );
    elementv = new FEAElement_Tet10( nqp_tet );
  }
  else if( elemType == FEType::Hex8 )
  {
    quadv = new QuadPts_Gauss_Hex( nqp_hex_1D );
    elementv = new FEAElement_Hex8( nqp_hex_1D * nqp_hex_1D * nqp_hex_1D );
  }
  else if( elemType == FEType::Hex27 )
  {
    quadv = new QuadPts_Gauss_Hex( nqp_hex_1D );
    elementv = new FEAElement_Hex27( nqp_hex_1D * nqp_hex_1D * nqp_hex_1D );
  }
  else SYS_T::print_fatal("Error: Invalid element type");

  PostVectSolution * pSolu = new PostVectSolution( sol_name,
      "node_mapping.h5", "post_node_mapping.h5", pNode, GMIptr->get_nFunc(), dof );

  int * IEN_e = new int[elementv->get_nLocBas()];
  double * ectrl_x = new double[elementv->get_nLocBas()];
  double * ectrl_y = new double[elementv->get_nLocBas()];
  double * ectrl_z = new double[elementv->get_nLocBas()];
  double * loc_disp_x = new double[elementv->get_nLocBas()];
  double * loc_disp_y = new double[elementv->get_nLocBas()];
  double * loc_disp_z = new double[elementv->get_nLocBas()];

  double subdomain_l2 = 0.0;
  double subdomain_H1 = 0.0;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    locIEN -> get_LIEN(ee, IEN_e);
    fNode -> get_ctrlPts_xyz( elementv->get_nLocBas(), IEN_e, ectrl_x, ectrl_y, ectrl_z );

    elementv -> buildBasis( quadv, ectrl_x, ectrl_y, ectrl_z );
    
    pSolu -> get_esol( 0, elementv->get_nLocBas(), IEN_e, loc_disp_x );
    pSolu -> get_esol( 1, elementv->get_nLocBas(), IEN_e, loc_disp_y );
    pSolu -> get_esol( 2, elementv->get_nLocBas(), IEN_e, loc_disp_z );

    // Calculate the error
    subdomain_l2 += POST_ERROR_E::get_manu_sol_errorL2(
      sol_time, loc_disp_x, loc_disp_y, loc_disp_z, elementv, ectrl_x, ectrl_y, ectrl_z, quadv );

    subdomain_H1 += POST_ERROR_E::get_manu_sol_errorH1(
      sol_time, loc_disp_x, loc_disp_y, loc_disp_z, elementv, ectrl_x, ectrl_y, ectrl_z, quadv );
  }
  
  double l2_error = 0.0; double H1_error = 0.0;
  MPI_Reduce(&subdomain_l2, &l2_error, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2_error = std::sqrt( l2_error );
 
  MPI_Reduce(&subdomain_H1, &H1_error, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  H1_error = std::sqrt( H1_error );

  SYS_T::commPrint("Error in L2 norm is : %e \n", l2_error);
  SYS_T::commPrint("Error in H1 norm is : %e \n", H1_error);

  delete elementv; delete quadv; delete locElem; delete pNode;
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] loc_disp_x; loc_disp_x = nullptr;
  delete [] loc_disp_y; loc_disp_y = nullptr;
  delete [] loc_disp_z; loc_disp_z = nullptr;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
