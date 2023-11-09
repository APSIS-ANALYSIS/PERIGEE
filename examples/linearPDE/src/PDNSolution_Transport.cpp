#include "PDNSolution_Transport.hpp"

PDNSolution_Transport:: PDNSolution_Transport( 
    const APart_Node * const &pNode,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 1, "Error: PDNSolution_Transport : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Transport: No such type of initial condition. \n");
      break;
  }
}

void PDNSolution_Transport::Init_zero( const APart_Node * const &pNode )
{
  const double value[1] = { 0.0 };

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode -> get_node_loc(ii);
    const int location[1] = { pos };

    VecSetValues(solution, 1, location, value, INSERT_VALUES);
  }

  Assembly_GhostUpdate();

  if( is_print )
  {
    std::ostringstream ss;
    ss<<"===> Initial "<<sol_name<<" solution vector: \n";
    SYS_T::commPrint(ss.str().c_str());
    SYS_T::commPrint("     val = 0.0 \n");
  }
}

// EOF
