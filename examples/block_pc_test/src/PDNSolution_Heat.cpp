#include "PDNSolution_Heat.hpp"

PDNSolution_Heat::PDNSolution_Heat( 
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const ALocal_InflowBC * const &infbc,
    const int &type, const bool &isprint ) 
: PDNSolution( pNode ), is_print(isprint)
{
  if( pNode->get_dof() != 1 ) SYS_T::print_fatal("Error: PDNSolution_Heat : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    case 1:
      Init_flow_parabolic( pNode, fNode_ptr, infbc );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Heat: No such type of initional condition.\n");
      break;
  }
}


PDNSolution_Heat::PDNSolution_Heat( 
    const APart_Node * const &pNode,
    const int &type, const bool &isprint ) 
: PDNSolution( pNode ), is_print( isprint )
{
  if( pNode->get_dof() != 1 ) SYS_T::print_fatal("Error: PDNSolution_Heat : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Heat : No such type of initional condition.\n");
      break;
  }
}


PDNSolution_Heat::~PDNSolution_Heat()
{}


void PDNSolution_Heat::Init_zero(const APart_Node * const &pNode_ptr)
{
  int location[1];
  const double value[1] = {0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 1;

    VecSetValues(solution, 1, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

  if( is_print )
    SYS_T::commPrint("===> Initial solution: val = 0.0 \n");
}


void PDNSolution_Heat::Init_flow_parabolic(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_InflowBC * const &infbc )
{
  int location[1];
  double value[1] = {0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 1;

    VecSetValues(solution, 1, location, value, INSERT_VALUES);
  }

  // If there are inflow nodes, set their value to be parabolic 
  if( infbc->get_Num_LD() > 0)
  {
    for(int ii=0; ii<nlocalnode; ++ii)
    {
      if( infbc->is_inLDN(pNode_ptr->get_node_loc(ii)) )
      {
        location[0] = pNode_ptr->get_node_loc(ii) * 1;

        const Vector_3 pt = fNode_ptr -> get_ctrlPts_xyz(ii);
        const double r =  infbc -> get_radius( pt );

        const double vmax = 1.0;
        const double vel = vmax * (1.0 - r*r);

        // -1.0 is multiplied to make the flow direction inward
        value[0] = vel;

        VecSetValues(solution, 1, location, value, INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

  if(is_print)
    SYS_T::commPrint("===> Initial solution: pres   = parabolic \n");
}

// EOF
