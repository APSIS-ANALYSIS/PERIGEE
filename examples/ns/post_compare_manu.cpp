#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "APart_Node.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "PostVectSolution.hpp"
#include "Post_error_Cubic.hpp"

int main( int argc, char * argv[] )
{
  int nqp_tet = 5; // 4 (2), 5 (3), 17 (5), or 29 (6)
  int nqp_vol_1D = 2;
  int nqp_vol = 8;

  std::string sol_name("SOL_900000000");
  std::string part_file("./part");
  
  const int dof = 4;
  double Q_thd = 1.0;
  double T_thd = 1.0;
  double tt = 0.0;
  double RR = 1.0;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);

  // Get options from command line
  SYS_T::GetOptionString("-sol_name", sol_name);
  SYS_T::GetOptionReal("-Q_thd",  Q_thd);
  SYS_T::GetOptionReal("-t_thd",  T_thd);
  SYS_T::GetOptionReal("-time", tt);
  SYS_T::GetOptionReal("-Radius", RR);
  
  // print options 
  SYS_T::cmdPrint("-sol_name:", sol_name);
  SYS_T::cmdPrint("-Q_thd",  Q_thd);
  SYS_T::cmdPrint("-T_thd",  T_thd);
  SYS_T::cmdPrint("-time", tt);
  SYS_T::cmdPrint("-Radius", RR);
  
  SYS_T::file_check(sol_name.c_str());
  double Qt = 0.0;
  if (tt < T_thd) Qt = Q_thd * 0.5 * (1-std::cos(MATH_T::PI*tt/T_thd));
  else Qt  = Q_thd;

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

  FEAElement * elementv = nullptr;
  IQuadPts * quadv = nullptr;
  int n_loc = 4;

  if( GMIptr->get_elemType() == 501 )
  {
    nqp_tet = 5;
    elementv = new FEAElement_Tet4( nqp_tet ); // elem type 501
    quadv = new QuadPts_Gauss_Tet( nqp_tet );
    n_loc = 4;
  }
  else if( GMIptr->get_elemType() == 502 )
  {
    nqp_tet = 17;
    elementv = new FEAElement_Tet10_v2( nqp_tet ); // elem type 502
    quadv = new QuadPts_Gauss_Tet( nqp_tet );
    n_loc = 10;
  }
  else if( GMIptr->get_elemType() == 601 )
  {
    nqp_vol_1D = 2;
    nqp_vol = 8;
    elementv = new FEAElement_Hex8( nqp_vol ); // elem type 601
    SYS_T::cmdPrint("-ElemType:", GMIptr->get_elemType() );
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
    n_loc = 8;
  }
  else if( GMIptr->get_elemType() == 602 )
  {
    nqp_vol_1D = 4;
    nqp_vol = 64;
    elementv = new FEAElement_Hex27( nqp_vol ); // elem type 602
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
    n_loc = 27;
  }
  else SYS_T::print_fatal("Error: Element type not supported.\n");


  PostVectSolution * pSolu = new PostVectSolution( sol_name,
      "node_mapping.h5", "post_node_mapping.h5", pNode, GMIptr->get_nFunc(), dof );

  // 
  const int nLoc_max = 27;
  int IEN_e[nLoc_max];
  double ectrl_x[nLoc_max], ectrl_y[nLoc_max], ectrl_z[nLoc_max];
  double loc_sol_ux[nLoc_max], loc_sol_uy[nLoc_max], loc_sol_uz[nLoc_max];

  double subdomain_l2_err = 0.0;
  double subdomain_l2_sol = 0.0;
  double subdomain_H1_err = 0.0;
  double subdomain_H1_sol = 0.0;
  double meshsize = 0.0;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
  {
    locIEN -> get_LIEN(ee, IEN_e);
    fNode -> get_ctrlPts_xyz( n_loc, IEN_e, ectrl_x, ectrl_y, ectrl_z );

    elementv -> buildBasis( quadv, ectrl_x, ectrl_y, ectrl_z );
    const double hh = elementv -> get_h( ectrl_x, ectrl_y, ectrl_z );
    if(hh > meshsize) meshsize = hh;
    
    pSolu -> get_esol( 1, n_loc, IEN_e, loc_sol_ux );
    pSolu -> get_esol( 2, n_loc, IEN_e, loc_sol_uy );
    pSolu -> get_esol( 3, n_loc, IEN_e, loc_sol_uz );

    // Compute the norms of error and solution
    subdomain_l2_err += POST_ERROR_C::get_manu_sol_u_errorL2(
      loc_sol_ux, loc_sol_uy, loc_sol_uz, elementv, ectrl_x, ectrl_y, ectrl_z, quadv,
      Qt, RR );

    subdomain_l2_sol += POST_ERROR_C::get_manu_sol_u_L2(
      elementv, ectrl_x, ectrl_y, ectrl_z, quadv,Qt, RR ); 

    subdomain_H1_err += POST_ERROR_C::get_manu_sol_u_errorH1(
      loc_sol_ux, loc_sol_uy, loc_sol_uz, elementv, ectrl_x, ectrl_y, ectrl_z, quadv,
      Qt, RR );

    subdomain_H1_sol += POST_ERROR_C::get_manu_sol_u_H1(
      elementv, ectrl_x, ectrl_y, ectrl_z, quadv, Qt, RR );
  }
  
  double l2_err = 0.0;  double l2_sol = 0.0; 
  double H1_err = 0.0;  double H1_sol = 0.0;

  MPI_Reduce(&subdomain_l2_err, &l2_err, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&subdomain_l2_sol, &l2_sol, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  l2_err = std::sqrt( l2_err );
  l2_sol = std::sqrt( l2_sol );
 
  MPI_Reduce(&subdomain_H1_err, &H1_err, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&subdomain_H1_sol, &H1_sol, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  H1_err = std::sqrt( H1_err );
  H1_sol = std::sqrt( H1_sol );
  
  SYS_T::commPrint("Relative error in L2 norm is : %e \n", l2_err/l2_sol);
  SYS_T::commPrint("Relative error in H1 norm is : %e \n", H1_err/H1_sol);
  SYS_T::commPrint("Meshsize is : %e \n", meshsize);

  delete elementv; delete quadv; delete locElem; delete pNode;
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
