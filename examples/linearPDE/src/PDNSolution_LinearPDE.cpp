#include "PDNSolution_LinearPDE.hpp"

PDNSolution_LinearPDE:: PDNSolution_LinearPDE( 
    const APart_Node * const &pNode,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 1 || pNode->get_dof() != 3, 
    "Error: PDNSolution_LinearPDE : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_LinearPDE: No such type of initial condition. \n");
      break;
  }
}

PDNSolution_LinearPDE::~PDNSolution_LinearPDE()
{}

void PDNSolution_LinearPDE::Init_zero( const APart_Node * const &pNode )
{
  const int dof = pNode->get_dof;
  const double value[ dof ] = {0.0};

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode -> get_node_loc(ii) * dof;
    int location[dof] = {0};
    for(int ii=0; ii<dof; ++ii)
    {
      location[dof+ii] = pos + ii;
    }

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
