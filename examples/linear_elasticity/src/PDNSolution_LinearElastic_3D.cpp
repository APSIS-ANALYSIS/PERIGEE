#include "PDNSolution_LinearElastic_3D.hpp"

PDNSolution_LinearElastic_3D::PDNSolution_LinearElastic_3D(
    const class APart_Node * const &pNode,
    const class IPLocAssem * const &locassem_ptr,
    const class FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution( pNode )
{
  switch(type)
  {
    case -1:
      Init_test_1( pNode, locassem_ptr, fNode_ptr );
      break;
    case 0:
      Init_zero_solu( pNode, locassem_ptr, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_LinearElastic_3D: No such type of initional condition. \n");
  }
}



PDNSolution_LinearElastic_3D::~PDNSolution_LinearElastic_3D()
{}


void PDNSolution_LinearElastic_3D::Init_test_1(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[3];
  double value[3] = {0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();
  
  double x, y, z;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);

    if(z == 1.0) value[2] = x*y;
    else value[2] = 0.0;

    location[0] = pNode_ptr->get_node_loc(ii) * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: Test case 1. \n");
}


void PDNSolution_LinearElastic_3D::Init_zero_solu(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[3];
  double value[3] = {0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
  SYS_T::commPrint("                       disp_y = 0.0 \n");
  SYS_T::commPrint("                       disp_z = 0.0 \n");
}



// EOF
