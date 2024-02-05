#include "PDNSolution_Smooth_Stress.hpp"

PDNSolution_Smooth_Stress::PDNSolution_Smooth_Stress(
    const APart_Node * const &pNode,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode, 6 ), sol_name( in_name ), is_print( isprint )
{
  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      break;
  }
}

void PDNSolution_Smooth_Stress::Init_zero( const APart_Node * const &pNode )
{
  const double value[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode ->  get_node_loc(ii) * 6;
    const int location[6] = { pos, pos+1, pos+2, pos+3, pos+4, pos+5 };

    VecSetValues(solution, 6, location, value, INSERT_VALUES);
  }

  Assembly_GhostUpdate();

  if( is_print )
  {
    std::ostringstream ss;
    ss<<"===> Initial "<<sol_name<<" solution vector: \n";
    SYS_T::commPrint(ss.str().c_str());
    SYS_T::commPrint("All components of the solution are zero \n");
  }
}

//EOF