#include "PDNSolution_NS.hpp"

PDNSolution_NS::PDNSolution_NS( 
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const ALocal_InflowBC * const &infbc,
    const int &type, const bool &isprint ) 
: PDNSolution( pNode ), is_print(isprint)
{
  if( pNode->get_dof() != 4 ) SYS_T::print_fatal("Error: PDNSolution_NS : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    case 1:
      Init_flow_parabolic( pNode, fNode_ptr, infbc );
      break;
    case 2:
      Init_pipe_parabolic( pNode, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_NS: No such type of initial condition.\n");
      break;
  }
}

PDNSolution_NS::PDNSolution_NS( 
    const APart_Node * const &pNode,
    const int &type, const bool &isprint ) 
: PDNSolution( pNode ), is_print( isprint )
{
  if( pNode->get_dof() != 4 ) SYS_T::print_fatal("Error: PDNSolution_NS : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_NS : No such type of initial condition.\n");
      break;
  }
}

PDNSolution_NS::~PDNSolution_NS()
{}

void PDNSolution_NS::Init_zero(const APart_Node * const &pNode_ptr)
{
  const double value[4] = {0.0, 0.0, 0.0, 0.0};

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 4;
    const int location[4] = { pos, pos + 1, pos +2, pos + 3 };
    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  Assembly_GhostUpdate();

  if( is_print )
  {
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = 0.0 \n");
    SYS_T::commPrint("                       velo_y = 0.0 \n");
    SYS_T::commPrint("                       velo_z = 0.0 \n");
  }
}

void PDNSolution_NS::Init_flow_parabolic(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_InflowBC * const &infbc )
{
  double value[4] = {0.0, 0.0, 0.0, 0.0};

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 4;
    const int location[4] = { pos, pos + 1, pos +2, pos + 3 };
    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = parabolic \n");
    SYS_T::commPrint("                       velo_y = parabolic \n");
    SYS_T::commPrint("                       velo_z = parabolic \n");
    SYS_T::commPrint("                       flow rate 1.0 .\n");
  }
  
  const int num_nbc = infbc->get_num_nbc();

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    // Maximum speed formula is 
    //             2.0 x flow rate (1.0) / surface area
    // Here I use the unit flow rate, and the actual flow rate is adjusted
    // based on the CVFlowRate class.
    const double vmax = 2.0 / infbc->get_fularea(nbc_id);

    const double out_nx = infbc->get_outvec(nbc_id).x();
    const double out_ny = infbc->get_outvec(nbc_id).y();
    const double out_nz = infbc->get_outvec(nbc_id).z();

    // If there are inflow nodes, set their value to be parabolic flow
    if( infbc->get_Num_LD(nbc_id) > 0)
    {
      for(int ii=0; ii<nlocalnode; ++ii)
      {
        if( infbc->is_inLDN(nbc_id, pNode_ptr->get_node_loc(ii)) )
        {
          const int pos = pNode_ptr->get_node_loc(ii) * 4;
          const int location[4] = { pos, pos + 1, pos +2, pos + 3 };

          const Vector_3 pt = fNode_ptr -> get_ctrlPts_xyz(ii);
          const double r = infbc -> get_radius( nbc_id, pt );
          const double vel = vmax * (1.0 - r*r);

          value[1] = vel * out_nx;
          value[2] = vel * out_ny;
          value[3] = vel * out_nz;

          VecSetValues(solution, 4, location, value, INSERT_VALUES);
        }
      }
    }

    if(is_print)
    {
      SYS_T::commPrint("                       -- nbc_id = %d \n", nbc_id);
      SYS_T::commPrint("                          max speed %e.\n", vmax);
      SYS_T::commPrint("                          active area is %e.\n", infbc->get_actarea(nbc_id) );
      SYS_T::commPrint("                          full area is %e.\n", infbc->get_fularea(nbc_id) );
      SYS_T::commPrint("                          outward normal direction [%e %e %e].\n", out_nx, out_ny, out_nz);
    }
  }

  Assembly_GhostUpdate();
}

void PDNSolution_NS::Init_pipe_parabolic(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
{
  double value[4] = {0.0, 0.0, 0.0, 0.0};

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 4;
    const int location[4] = { pos, pos + 1, pos +2, pos + 3 };

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  // Maximum speed formula is 
  //             2.0 x flow rate (1.0) / surface area
  // Here I use the unit flow rate, and the actual flow rate is adjusted
  // based on the CVFlowRate class.
  const double vmax = 2.0 * 100.0 / (3.14 * 4);

  const double out_nx = 0.0; 
  const double out_ny = 0.0;
  const double out_nz = 1.0;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 4;
    const int location[4] = { pos, pos + 1, pos +2, pos + 3 };

    const double x = fNode_ptr->get_ctrlPts_x(ii);
    const double y = fNode_ptr->get_ctrlPts_y(ii);
    //const double z = fNode_ptr->get_ctrlPts_z(ii);

    const double r = std::sqrt(x*x + y*y) / 2.0;
    const double vel = vmax * (1.0 - r*r);

    // -1.0 is multiplied to make the flow direction inward
    value[1] = vel * out_nx;
    value[2] = vel * out_ny;
    value[3] = vel * out_nz;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  Assembly_GhostUpdate();

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = parabolic \n");
    SYS_T::commPrint("                       velo_y = parabolic \n");
    SYS_T::commPrint("                       velo_z = parabolic \n");
    SYS_T::commPrint("                       flow rate 100.0 .\n");
    SYS_T::commPrint("                       direction [%e %e %e].\n", out_nx, out_ny, out_nz);
  }
}

// EOF
