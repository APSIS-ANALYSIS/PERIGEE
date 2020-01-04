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
    case 1:
      Init_inflow( pNode, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Mixed_UPV_3D: No such type of initional condition. \n");
      break;
  }
}


PDNSolution_Mixed_UPV_3D::~PDNSolution_Mixed_UPV_3D()
{}


void PDNSolution_Mixed_UPV_3D::Init_zero(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
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


void PDNSolution_Mixed_UPV_3D::Init_inflow(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double val = 51.3; // inflow 51.3 cm / sec

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    double x = fNode_ptr->get_ctrlPts_x(ii);

    if( MATH_T::equals(x, -5.0, 1.0e-5) ) value[4] = val;
    else value[4] = 0.0;

    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
  SYS_T::commPrint("                       disp_y = 0.0 \n");
  SYS_T::commPrint("                       disp_z = 0.0 \n");
  SYS_T::commPrint("                       pres   = 0.0 \n");
  SYS_T::commPrint("                       velo_x = val \n");
  SYS_T::commPrint("                       velo_y = 0.0 \n");
  SYS_T::commPrint("                       velo_z = 0.0 \n");
}

// EOF
