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
    const FEANode * const &fNode,
    const ALocal_InflowBC * const &infbc,
    const PDNSolution * const &curr_disp,
    const std::vector<double> &curr_area,
    const std::vector<Vector_3> &curr_centroid,
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
      Init_flow_parabolic( pNode, infbc, curr_disp, curr_area, curr_centroid );
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
    const ALocal_InflowBC * const &infbc,
    const PDNSolution * const &curr_disp,
    const std::vector<double> &curr_area,
    const std::vector<Vector_3> &curr_centroid )
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
    const double vmax = 2.0 / curr_area[nbc_id];
    const double out_nx = infbc->get_outvec( nbc_id ).x();
    const double out_ny = infbc->get_outvec( nbc_id ).y();
    const double out_nz = infbc->get_outvec( nbc_id ).z();

    // If this sub-domain contains inflow nodes, set their values based on the
    // parabolic flow profile
    if( infbc->get_Num_LD( nbc_id ) > 0)
    {
      const std::vector<double> local_d = curr_disp -> GetLocalArray();

      std::vector<Vector_3> inner_points {};
      std::vector<Vector_3> ring_points {};

      std::vector<int> LD_loc_tag {};

      for(int jj=0; jj<infbc->get_num_local_node(nbc_id); ++jj)
      {   
        // Update all inlet points
        Vector_3 pt = infbc->get_local_pt_xyz(nbc_id, jj);

        const int pt_pos = infbc->get_local_node_pos(nbc_id, jj);

        pt.x() += local_d[3 * pt_pos    ];
        pt.y() += local_d[3 * pt_pos + 1];
        pt.z() += local_d[3 * pt_pos + 2];

        // pick out LD point
        if( infbc->is_inLDN( nbc_id, pNode_ptr->get_node_loc(pt_pos)) )
        {
          inner_points.push_back(pt);
          LD_loc_tag.push_back(pNode_ptr->get_node_loc(pt_pos));
        }
        else
          ring_points.push_back(pt);
      }

      // Apply inflow BC at LD points
      for(int ii=0; ii < VEC_T::get_size(LD_loc_tag); ++ii)
      {
        const int pos = LD_loc_tag[ii];
        const int location[3] = { pos, pos + 1, pos + 2 };

        const Vector_3 pt = inner_points[ii];

        const double r = get_curr_radius(ring_points, pt, curr_centroid[nbc_id]);
        const double vel = vmax * (1.0 - r*r);

        const double value[3] = { vel * out_nx, vel * out_ny, vel * out_nz };

        VecSetValues(solution, 3, location, value, INSERT_VALUES);
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
      SYS_T::commPrint("        max speed %e.\n", 2.0 / curr_area[nbc_id] );
      SYS_T::commPrint("        active area is %e.\n", infbc->get_actarea(nbc_id) );
      SYS_T::commPrint("        full area is %e (%e).\n", infbc->get_fularea(nbc_id), curr_area[nbc_id] );
      SYS_T::commPrint("        outward normal direction [%e %e %e].\n",
          infbc->get_outvec( nbc_id ).x(), infbc->get_outvec( nbc_id ).y(), infbc->get_outvec( nbc_id ).z() );
    }
  }
}

double PDNSolution_V::get_curr_radius( const std::vector<Vector_3> &outline_pt,
    const Vector_3 &target_pt,
    const Vector_3 &centroid)
{
  const double rc = Vec3::dist(target_pt, centroid);

  double rb = Vec3::dist(target_pt, outline_pt[0]);

  for(int ii=1; ii<VEC_T::get_size(outline_pt); ++ii)
  {
    const double newdist = Vec3::dist(target_pt, outline_pt[ii]);

    if(newdist < rb) rb = newdist;
  }

  return rc / (rb + rc);
}

// EOF
