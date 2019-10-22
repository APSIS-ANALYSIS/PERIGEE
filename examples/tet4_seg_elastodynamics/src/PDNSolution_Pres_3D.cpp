#include "PDNSolution_Pres_3D.hpp"

PDNSolution_Pres_3D::PDNSolution_Pres_3D(
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution( pNode, 1 )
{
  switch( type )
  {
    case 0:
      Init_zero_solu( pNode, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Pres_3D: No such type of initial condition. \n");
      break;
  }
}


PDNSolution_Pres_3D::~PDNSolution_Pres_3D()
{}


void PDNSolution_Pres_3D::Init_zero_solu( const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
{
  int location;
  const double value = 0.0;
  int nlocalnode = pNode_ptr -> get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location = pNode_ptr -> get_node_loc(ii);
    VecSetValue(solution, location, value, INSERT_VALUES);
  }
  
  VecAssemblyBegin(solution);  VecAssemblyEnd(solution);
  GhostUpdate();
  //SYS_T::commPrint("===> Initial solution: pres = 0.0 \n");
}

// EOF
