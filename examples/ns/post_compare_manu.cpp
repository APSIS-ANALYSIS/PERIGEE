#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "APart_Node.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Tet4.hpp"
#include "PostVectSolution.hpp"
#include "Post_error_Poiseuille.hpp"

int main( int argc, char * argv[] )
{
  int nqp_tet = 29; // 4 (2), 5 (3), 17 (5), or 29 (6)

  double sol_time = 1.0;
  std::string sol_name("SOL_900000000");
  std::string part_file("./part");
  const int dof = 4;

  double FlowRate = 0.125663706144; // 0.04 * PI
  double Radius = 0.1;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  SYS_T::GetOptionString("-sol_name", sol_name);
  SYS_T::GetOptionReal("-sol_time",   sol_time);
  SYS_T::GetOptionReal("-flowrate",   FlowRate);
  SYS_T::GetOptionReal("-radius",     Radius);

  SYS_T::cmdPrint("-sol_name:", sol_name);
  SYS_T::cmdPrint("-sol_time:", sol_time);
  SYS_T::cmdPrint("-flowrate:", FlowRate);
  SYS_T::cmdPrint("-radius:",   Radius);

  SYS_T::file_check(sol_name.c_str());

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  FEANode * fNode = new FEANode(part_file, rank);

  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::print_fatal_if( size != PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for:", size);
  SYS_T::commPrint("Postprocessing - compute error from manufactured solutions.\n");

  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );

  FEAElement * elementv = new FEAElement_Tet4( nqp_tet );

  PostVectSolution * pSolu = new PostVectSolution( sol_name,
      "node_mapping.h5", "post_node_mapping.h5", pNode, GMIptr->get_nFunc(), dof );

  // There are four local basis functions in tet element
  int IEN_e[4];
  double ectrl_x[4], ectrl_y[4], ectrl_z[4];
  double loc_sol_ux[4], loc_sol_uy[4], loc_sol_uz[4];
  double subdomain_l2 = 0.0;
  double subdomain_H1 = 0.0;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    locIEN -> get_LIEN(ee, IEN_e);
    fNode -> get_ctrlPts_xyz( 4, IEN_e, ectrl_x, ectrl_y, ectrl_z );

    elementv -> buildBasis( quadv, ectrl_x, ectrl_y, ectrl_z );
    
    pSolu -> get_esol( 1, 4, IEN_e, loc_sol_ux );
    pSolu -> get_esol( 2, 4, IEN_e, loc_sol_uy );
    pSolu -> get_esol( 3, 4, IEN_e, loc_sol_uz );

    // Calculate the error
    subdomain_l2 += POST_ERROR_P::get_manu_sol_u_error(
      loc_sol_ux, loc_sol_uy, loc_sol_uz, elementv, ectrl_x, ectrl_y, ectrl_z, quadv, FlowRate, Radius );

    subdomain_H1 += POST_ERROR_P::get_manu_sol_u_errorH1(
      loc_sol_ux, loc_sol_uy, loc_sol_uz, elementv, ectrl_x, ectrl_y, ectrl_z, quadv, FlowRate, Radius );
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
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF