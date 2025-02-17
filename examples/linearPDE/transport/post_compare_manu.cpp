// ==================================================================
// post_compare_manu.cpp
//
// Error analysis code for 3D transport.
//
// Data: Oct. 25 2023
// ==================================================================
#include "AGlobal_Mesh_Info.hpp"
#include "ANL_Tools.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "APart_Node.hpp"
#include "FEANode.hpp"
#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"
#include "PostVectSolution.hpp"
#include "Post_error_transport.hpp"

int main( int argc, char * argv[] )
{
  int nqp_vol = 29;
  double sol_time = 1.0;
  std::string sol_name("SOL_900000000");
  std::string part_file("./ppart/part");
  const int dof = 1;

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  SYS_T::GetOptionString("-sol_name", sol_name);
  SYS_T::GetOptionReal("-sol_time",   sol_time);
  SYS_T::GetOptionInt("-nqp_vol",     nqp_vol);

  SYS_T::cmdPrint("-nqp_vol:", nqp_vol);
  SYS_T::cmdPrint("-sol_name:", sol_name);
  SYS_T::cmdPrint("-sol_time:", sol_time);

  SYS_T::file_check(sol_name.c_str());

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  auto locIEN  = SYS_T::make_unique<ALocal_IEN>(part_file, rank);
  auto locElem = SYS_T::make_unique<ALocal_Elem>(part_file, rank);
  auto fNode   = SYS_T::make_unique<FEANode>(part_file, rank);
  auto pNode   = SYS_T::make_unique<APart_Node>(part_file, rank);
  auto GMIptr  = SYS_T::make_unique<AGlobal_Mesh_Info>(part_file,rank);

  SYS_T::print_fatal_if( size != ANL_T::get_cpu_size(part_file, rank),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for:", size);
  SYS_T::commPrint("Postprocessing - compute error from manufactured solutions.\n");

  const FEType elemType = GMIptr -> get_elemType();
  std::unique_ptr<IQuadPts> quadv = QuadPtsFactory::createVolQuadrature(elemType, nqp_vol);
  std::unique_ptr<FEAElement> elementv = ElementFactory::createVolElement(elemType, nqp_vol);

  PostVectSolution * pSolu = new PostVectSolution( sol_name,
      "node_mapping.h5", "post_node_mapping.h5", pNode.get(), GMIptr->get_nFunc(), dof );

  // There are four local basis functions in tet element
  int * IEN_e = new int[elementv->get_nLocBas()];
  double * ectrl_x = new double[elementv->get_nLocBas()];
  double * ectrl_y = new double[elementv->get_nLocBas()];
  double * ectrl_z = new double[elementv->get_nLocBas()];
  double * loc_sol = new double[elementv->get_nLocBas()];

  double subdomain_l2 = 0.0;
  double subdomain_H1 = 0.0;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    locIEN -> get_LIEN(ee, IEN_e);
    fNode -> get_ctrlPts_xyz( elementv->get_nLocBas(), IEN_e, ectrl_x, ectrl_y, ectrl_z );

    elementv -> buildBasis( quadv.get(), ectrl_x, ectrl_y, ectrl_z );
    
    pSolu -> get_esol( 0, elementv->get_nLocBas(), IEN_e, loc_sol );

    // Calculate the error
    subdomain_l2 += POST_ERROR_T::get_manu_sol_errorL2(
      sol_time, loc_sol, elementv.get(), ectrl_x, ectrl_y, ectrl_z, quadv.get() );

    subdomain_H1 += POST_ERROR_T::get_manu_sol_errorH1(
      sol_time, loc_sol, elementv.get(), ectrl_x, ectrl_y, ectrl_z, quadv.get() );
  }
  
  double l2_error = 0.0; double H1_error = 0.0;
  MPI_Reduce(&subdomain_l2, &l2_error, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2_error = std::sqrt( l2_error );
 
  MPI_Reduce(&subdomain_H1, &H1_error, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  H1_error = std::sqrt( H1_error );

  SYS_T::commPrint("Error in L2 norm is : %e \n", l2_error);
  SYS_T::commPrint("Error in H1 norm is : %e \n", H1_error);

  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] loc_sol; loc_sol = nullptr;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
