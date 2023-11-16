#include "PDNSolution_Elastodynamics.hpp"

PDNSolution_Elastodynamics:: PDNSolution_Elastodynamics( 
    const APart_Node * const &pNode,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 3, "Error: PDNSolution_Elastodynamics : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Elastodynamics: No such type of initial condition. \n");
      break;
  }
}

void PDNSolution_Elastodynamics::Init_zero( const APart_Node * const &pNode )
{
  const double value[3] = { 0.0, 0.0, 0.0 };

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode -> get_node_loc(ii);
    const int location[3] = { pos, pos+1, pos+2 };

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  Assembly_GhostUpdate();

  if( is_print )
  {
    std::ostringstream ss;
    ss<<"===> Initial "<<sol_name<<" solution vector: \n";
    SYS_T::commPrint(ss.str().c_str());
    SYS_T::commPrint("     ux = 0.0 \n");
    SYS_T::commPrint("     uy = 0.0 \n");
    SYS_T::commPrint("     uz = 0.0 \n");
  }
}

// EOF