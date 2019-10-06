#include "PDNSolution_Mixed_LinearElastic_3D.hpp"

PDNSolution_Mixed_LinearElastic_3D::PDNSolution_Mixed_LinearElastic_3D(
    const APart_Node * const &pNode,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution( pNode )
{
  switch( type )
  {
    case 0:
      Init_zero_solu( pNode, locassem_ptr, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Mixed_LinearElastic_3D: No such type of initional condition. \n");
  }
}


PDNSolution_Mixed_LinearElastic_3D::~PDNSolution_Mixed_LinearElastic_3D()
{}


void PDNSolution_Mixed_LinearElastic_3D::Init_zero_solu( 
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4] = {0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr -> get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }
  VecAssemblyBegin(solution);  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
  SYS_T::commPrint("                       disp_x = 0.0 \n");
  SYS_T::commPrint("                       disp_y = 0.0 \n");
  SYS_T::commPrint("                       disp_z = 0.0 \n");
}

// EOF
