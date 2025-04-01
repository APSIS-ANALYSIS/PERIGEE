#include "PDNSolution_P.hpp"

PDNSolution_P::PDNSolution_P( const APart_Node * const &pNode,
    const int &type, const int &input_dof_num,
    const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode, input_dof_num ), sol_name( in_name ), is_print( isprint )
{
  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_P: No such type of initial condition. \n");
      break;
  }
}

void PDNSolution_P::Init_zero( const APart_Node * const &pNode )
{
  const double val = 0.0;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode -> get_node_loc(ii);
    VecSetValue(solution, pos, val, INSERT_VALUES);
  }
  
  Assembly_GhostUpdate();

  if( is_print )
  {
    std::ostringstream ss;
    ss<<"===> Initial "<<sol_name<<" solution vector: \n";
    SYS_T::commPrint(ss.str().c_str());
    SYS_T::commPrint("     val_x = 0.0 \n");
    SYS_T::commPrint("     val_y = 0.0 \n");
    SYS_T::commPrint("     val_z = 0.0 \n");
  }
}

// EOF
