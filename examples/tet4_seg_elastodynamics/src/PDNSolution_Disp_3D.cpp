#include "PDNSolution_Disp_3D.hpp"

PDNSolution_Disp_3D::PDNSolution_Disp_3D(
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution( pNode, 3 )
{
  switch( type )
  {
    case 0:
      Init_zero_solu( pNode, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Disp_3D: No such type of initial condition. \n");
      break;
  }
}


PDNSolution_Disp_3D::~PDNSolution_Disp_3D()
{}


void PDNSolution_Disp_3D::Init_zero_solu( 
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[3];
  double value[3] = { 0.0, 0.0, 0.0 };
  int nlocalnode = pNode_ptr -> get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
  
    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);  VecAssemblyEnd(solution);
  GhostUpdate();
  //SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
  //SYS_T::commPrint("                       disp_y = 0.0 \n");
  //SYS_T::commPrint("                       disp_z = 0.0 \n");
}

// EOF
