#include "PDNSolution_V.hpp"

PDNSolution_V::PDNSolution_V( const APart_Node * const &pNode,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 3, "Error: PDNSolution_V : the APart_Node gives wrong dof number. \n");
  
  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_V: No such type of initial condition. \n");
      break;
  }
}

PDNSolution_V::PDNSolution_V( const APart_Node * const &pNode,
    const int &type, const int &input_dof_num,
    const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode, input_dof_num ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 3, "Error: PDNSolution_V : the APart_Node gives wrong dof number. \n");
  
  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_V: No such type of initial condition. \n");
      break;
  }
}

PDNSolution_V::PDNSolution_V( const APart_Node * const &pNode,
    const FEANode * const &fNode,
    const ALocal_InflowBC * const &infbc,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 3, "Error: PDNSolution_V : the APart_Node gives wrong dof number. \n");
  
  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    case 1:
      Init_flow_parabolic( pNode, fNode, infbc );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_V: No such type of initial condition. \n");
      break;
  }
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

void PDNSolution_V::Init_flow_parabolic( const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_InflowBC * const &infbc )
{
  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 3;
    const int location[3] = { pos, pos + 1, pos + 2 };
    const double value[3] = {0.0, 0.0, 0.0};

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  const int num_nbc = infbc -> get_num_nbc();

  for(int nbc_id = 0; nbc_id < num_nbc; ++nbc_id)
  {
    const double vmax = 2.0 / infbc->get_fularea( nbc_id );
    const double out_nx = infbc->get_outvec( nbc_id ).x();
    const double out_ny = infbc->get_outvec( nbc_id ).y();
    const double out_nz = infbc->get_outvec( nbc_id ).z();

    // If this sub-domain contains inflow nodes, set their values based on the
    // parabolic flow profile
    if( infbc->get_Num_LD( nbc_id ) > 0)
    {
      for(int ii=0; ii<nlocalnode; ++ii)
      {
        if( infbc->is_inLDN( nbc_id, pNode_ptr->get_node_loc(ii) ) )
        {
          const int pos = pNode_ptr->get_node_loc(ii) * 3;
          const int location[3] = { pos, pos + 1, pos + 2 };

          const Vector_3 pt = fNode_ptr -> get_ctrlPts_xyz(ii);
          const double r =  infbc -> get_radius( nbc_id, pt );
          const double vel = vmax * (1.0 - r*r);

          const double value[3] = { vel * out_nx, vel * out_ny, vel * out_nz };

          VecSetValues(solution, 3, location, value, INSERT_VALUES);
        }
      }
    }
  }

  Assembly_GhostUpdate();

  if( is_print )
  {
    std::ostringstream ss;
    ss<<"===> Initial "<<sol_name<<" solution vector: \n";
    SYS_T::commPrint(ss.str().c_str());
    for(int nbc_id=0; nbc_id < num_nbc; ++nbc_id)
    {
      SYS_T::commPrint("     -- nbc_id = %d \n", nbc_id);
      SYS_T::commPrint("        max speed %e.\n", 2.0 / infbc->get_fularea( nbc_id ) );
      SYS_T::commPrint("        active area is %e.\n", infbc->get_actarea(nbc_id) );
      SYS_T::commPrint("        full area is %e.\n", infbc->get_fularea(nbc_id) );
      SYS_T::commPrint("        outward normal direction [%e %e %e].\n",
          infbc->get_outvec( nbc_id ).x(), infbc->get_outvec( nbc_id ).y(), infbc->get_outvec( nbc_id ).z() );
    }
  }
}

// EOF
