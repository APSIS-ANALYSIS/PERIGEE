#include "PDNSolution_V.hpp"

PDNSolution_V::PDNSolution_V( const APart_Node * const &pNode )
: PDNSolution( pNode, 3 )
{
  Init_zero( pNode );
}

void PDNSolution_V::Init_zero( const APart_Node * const &pNode )
{
  double value[3] = {0.0, 0.0, 0.0};

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode -> get_node_loc(ii) * 3;
    const int location[3] = { pos, pos + 1, pos + 2 };

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  Assembly_GhostUpdate();
}


// EOF
