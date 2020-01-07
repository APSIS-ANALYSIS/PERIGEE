#include "PDNSolution_Mixed_UPV_3D.hpp"

PDNSolution_Mixed_UPV_3D::PDNSolution_Mixed_UPV_3D(
    const class APart_Node * const &pNode,
    const class FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution( pNode )
{
  if( pNode->get_dof() != 7 ) SYS_T::print_fatal("Error: PDNSolution_Mixed_UPV_3D : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode, fNode_ptr );
      break;
    case 2:
      Init_pressure( pNode, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Mixed_UPV_3D: No such type of initial condition. \n");
      break;
  }
}


PDNSolution_Mixed_UPV_3D::PDNSolution_Mixed_UPV_3D(
    const class APart_Node * const &pNode,
    const class FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc,
    const int &type ) : PDNSolution( pNode )
{
  if( pNode->get_dof() != 7 ) SYS_T::print_fatal("Error: PDNSolution_Mixed_UPV_3D : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 1:
      Init_flow_parabolic( pNode, fNode_ptr, infbc );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Mixed_UPV_3D: No such type of initial condition. \n");
      break;
  }
}


PDNSolution_Mixed_UPV_3D::PDNSolution_Mixed_UPV_3D(
    const APart_Node * const &pNode,
    const int &type ) : PDNSolution( pNode )
{
  if( pNode->get_dof() != 7 ) SYS_T::print_fatal("Error: PDNSolution_Mixed_UPV_3D : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Mixed_UPV_3D: No such type of initial condition. \n");
      break;
  }
}


PDNSolution_Mixed_UPV_3D::~PDNSolution_Mixed_UPV_3D()
{}


void PDNSolution_Mixed_UPV_3D::Init_zero(
    const APart_Node * const &pNode_ptr )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
  SYS_T::commPrint("                       disp_y = 0.0 \n");
  SYS_T::commPrint("                       disp_z = 0.0 \n");
  SYS_T::commPrint("                       pres   = 0.0 \n");
  SYS_T::commPrint("                       velo_x = 0.0 \n");
  SYS_T::commPrint("                       velo_y = 0.0 \n");
  SYS_T::commPrint("                       velo_z = 0.0 \n");
}


void PDNSolution_Mixed_UPV_3D::Init_flow_parabolic( 
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();
  double x, y, z, r, vel;

  const double vmax = 2.0 / infbc->get_fularea();

  const double out_nx = infbc->get_outvec(0);
  const double out_ny = infbc->get_outvec(1);
  const double out_nz = infbc->get_outvec(2);

  if( infbc->get_Num_LD() > 0)
  {
    for(int ii=0; ii<nlocalnode; ++ii)
    {
      if( infbc->is_inLDN(pNode_ptr->get_node_loc(ii)) )
      {
        location[0] = pNode_ptr->get_node_loc(ii) * 7;
        location[1] = location[0] + 1;
        location[2] = location[0] + 2;
        location[3] = location[0] + 3;
        location[4] = location[0] + 4;
        location[5] = location[0] + 5;
        location[6] = location[0] + 6;

        x = fNode_ptr->get_ctrlPts_x(ii);
        y = fNode_ptr->get_ctrlPts_y(ii);
        z = fNode_ptr->get_ctrlPts_z(ii);

        r = infbc->get_radius(x,y,z);
        vel = vmax * (1.0 - r*r);

        value[4] = vel * (-1.0) * out_nx;
        value[5] = vel * (-1.0) * out_ny;
        value[6] = vel * (-1.0) * out_nz;

        VecSetValues(solution, 7, location, value, INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
  SYS_T::commPrint("                       velo_x = parabolic \n");
  SYS_T::commPrint("                       velo_y = parabolic \n");
  SYS_T::commPrint("                       velo_z = parabolic \n");
  PetscPrintf(PETSC_COMM_WORLD,
      "                       flow rate 1.0 .\n");
  PetscPrintf(PETSC_COMM_WORLD,
      "                       max speed %e.\n", vmax);
  PetscPrintf(PETSC_COMM_WORLD,
      "                       active area is %e.\n", infbc->get_actarea() );
  PetscPrintf(PETSC_COMM_WORLD,
      "                       full area is %e.\n", infbc->get_fularea() );
  PetscPrintf(PETSC_COMM_WORLD,
      "                       direction [%e %e %e].\n", out_nx, out_ny, out_nz);
}


void PDNSolution_Mixed_UPV_3D::Init_pressure(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode_fluid();

  value[3] = 5.0 * 1333.2239; // prescribe the pressure value 
  // 5 mmHg to dyns/cm^2
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int loc_node = pNode_ptr->get_node_loc_fluid(ii);
    location[0] = pNode_ptr->get_node_loc( loc_node ) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
  SYS_T::commPrint("                       disp_y = 0.0 \n");
  SYS_T::commPrint("                       disp_z = 0.0 \n");
  SYS_T::commPrint("                       pres   = prescribed value \n");
  SYS_T::commPrint("                       velo_x = 0.0 \n");
  SYS_T::commPrint("                       velo_y = 0.0 \n");
  SYS_T::commPrint("                       velo_z = 0.0 \n");
}

// EOF
