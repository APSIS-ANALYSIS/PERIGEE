#include "AGlobal_Mesh_Info.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "APart_Node.hpp"
#include "ANL_Tools.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEANode.hpp"
#include "PostVectSolution.hpp"
#include "Post_error_ES.hpp"
#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"

int main( int argc, char * argv[] )
{
  int nqp_v = 5;

  std::string sol_name("SOL_900000000");
  std::string part_file("./part");
  
  const int dof = 4;
  double tt = 0.0;

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  // Get options from command line
  SYS_T::GetOptionString("-sol_name", sol_name);
  SYS_T::GetOptionInt("-nqp_v", nqp_v);
  SYS_T::GetOptionReal("-time", tt);
  
  // print options 
  SYS_T::cmdPrint("-sol_name:", sol_name);
  SYS_T::cmdPrint("-nqp_v", nqp_v);
  SYS_T::cmdPrint("-time", tt);
  
  SYS_T::file_check(sol_name.c_str());

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  FEANode * fNode = new FEANode(part_file, rank);

  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  AGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info(part_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::print_fatal_if( size != ANL_T::get_cpu_size(part_file, rank),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for:", size);
  SYS_T::commPrint("Postprocessing - compute energy from solutions.\n");

  int n_loc = 4;

  if(  GMIptr->get_elemType() == FE_T::to_FEType("Tet4") )
    n_loc = 4;
  else if( GMIptr->get_elemType() == FE_T::to_FEType("Tet10") )
    n_loc = 10;
  else if( GMIptr->get_elemType() == FE_T::to_FEType("Hex8") )
    n_loc = 8;
  else if( GMIptr->get_elemType() == FE_T::to_FEType("Hex27") )
    n_loc = 27;
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  const FEType elemType(GMIptr->get_elemType());
  const int nqpv = nqp_v;

  const std::unique_ptr<FEAElement> elementv = ElementFactory::createVolElement(elemType, nqpv);
  const std::unique_ptr<IQuadPts> quadv = QuadPtsFactory::createVolQuadrature(elemType, nqpv);
  
  PostVectSolution * pSolu = new PostVectSolution( sol_name,
      "node_mapping.h5", "node_mapping.h5", pNode, GMIptr->get_nFunc(), dof );

  const int nLoc_max = 27;
  int IEN_e[nLoc_max];
  double ectrl_x[nLoc_max], ectrl_y[nLoc_max], ectrl_z[nLoc_max];
  double loc_sol_ux[nLoc_max], loc_sol_uy[nLoc_max], loc_sol_uz[nLoc_max];
  double loc_sol_p[nLoc_max];
  
  double subdomain_l2_err = 0.0;
  double subdomain_l2_sol = 0.0;
  double subdomain_p_err = 0.0;
  double subdomain_p_sol = 0.0;
  double meshsize = 0.0;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    locIEN -> get_LIEN(ee, IEN_e);
    fNode -> get_ctrlPts_xyz( n_loc, IEN_e, ectrl_x, ectrl_y, ectrl_z );

    elementv -> buildBasis( quadv.get(), ectrl_x, ectrl_y, ectrl_z );
    const double hh = elementv -> get_h( ectrl_x, ectrl_y, ectrl_z );
    if(hh > meshsize) meshsize = hh;
    
    pSolu -> get_esol( 0, n_loc, IEN_e, loc_sol_p );
    pSolu -> get_esol( 1, n_loc, IEN_e, loc_sol_ux );
    pSolu -> get_esol( 2, n_loc, IEN_e, loc_sol_uy );
    pSolu -> get_esol( 3, n_loc, IEN_e, loc_sol_uz );

    // Compute the norms of error and solution
    subdomain_l2_err += POST_ERROR_C::get_manu_sol_u_errorL2(
      tt, loc_sol_ux, loc_sol_uy, loc_sol_uz, elementv.get(), ectrl_x, ectrl_y, ectrl_z, quadv.get() );

    subdomain_l2_sol += POST_ERROR_C::get_manu_sol_u_L2(
      tt, elementv.get(), ectrl_x, ectrl_y, ectrl_z, quadv.get() ); 

    subdomain_p_err += POST_ERROR_C::get_manu_sol_p_errorL2(
      tt, loc_sol_p, elementv.get(), ectrl_x, ectrl_y, ectrl_z, quadv.get() );

    subdomain_p_sol += POST_ERROR_C::get_manu_sol_p_L2(
      tt, elementv.get(), ectrl_x, ectrl_y, ectrl_z, quadv.get() );
  }
  
  double l2_err = 0.0;  double l2_sol = 0.0; 
  double p_err = 0.0;  double p_sol = 0.0;

  MPI_Reduce(&subdomain_l2_err, &l2_err, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&subdomain_l2_sol, &l2_sol, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2_err = std::sqrt( l2_err );
  l2_sol = std::sqrt( l2_sol );

  MPI_Reduce(&subdomain_p_err, &p_err, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&subdomain_p_sol, &p_sol, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  p_err = std::sqrt( p_err );
  p_sol = std::sqrt( p_sol );
  
  SYS_T::commPrint("Relative error for vleo in L2 norm is : %e \n", l2_err/l2_sol);
  SYS_T::commPrint("Relative error for pres in L2 norm is : %e \n", p_err/p_sol);
  SYS_T::commPrint("Meshsize is : %e \n", meshsize);

  delete locElem; delete pNode; delete pSolu;
  delete fNode; delete locIEN; delete GMIptr;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
