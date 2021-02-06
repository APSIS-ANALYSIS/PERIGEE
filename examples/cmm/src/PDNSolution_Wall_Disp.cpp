#include "PDNSolution_Wall_Disp.hpp"

PDNSolution_Wall_Disp::PDNSolution_Wall_Disp( 
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const int &type, const bool &isprint = true )
: PDNSolution( pNode, 3 ), is_print( isprint )
{
  switch( type )
  {
    case 0:
      Init_zero( pNode, fNode_ptr );
    default:
      SYS_T::print_fatal("Error: in PDNSolution_Wall_Disp, No such type of initial condition. \n");
      break;
  }
}


PDNSolution_Wall_Disp::~PDNSolution_Wall_Disp()
{}


void PDNSolution_Wall_Disp::Init_zero( const APart_Node * const &pNode_ptr )
{
  int location[3];
  const double value[3] = {0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);

  GhostUpdate();

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
    SYS_T::commPrint("                       disp_y = 0.0 \n");
    SYS_T::commPrint("                       disp_z = 0.0 \n");
  }
}

// EOF
