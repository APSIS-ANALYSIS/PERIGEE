#include "PDNSolution_Solid.hpp"

PDNSolution_Solid::PDNSolution_Solid( const APart_Node * const &pNode,
    const int &in_dof_num, const int &type,
    const bool &isprint, const std::string &in_name )
: PDNSolution( pNode, in_dof_num ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( in_dof_num < 1, "Error: PDNSolution_Solid : invalid dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode, in_dof_num );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Solid: No such type of initial condition. \n");
      break;
  }
}

void PDNSolution_Solid::Init_zero( const APart_Node * const &pNode,
    const int &in_dof_num )
{
  std::vector<int> location(in_dof_num);
  std::vector<double> value(static_cast<size_t>(in_dof_num), 0.0);

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode->get_node_loc(ii) * in_dof_num;
    for(int kk=0; kk<in_dof_num; ++kk) location[kk] = pos + kk;

    VecSetValues(solution, in_dof_num, location.data(), value.data(), INSERT_VALUES);
  }

  Assembly_GhostUpdate();

  if( is_print )
  {
    std::ostringstream ss;
    ss<<"===> Initial "<<sol_name<<" solution vector: \n";
    SYS_T::commPrint(ss.str().c_str());
    for(int kk=0; kk<in_dof_num; ++kk)
      SYS_T::commPrint("     val[%d] = 0.0 \n", kk);
  }
}

// EOF
