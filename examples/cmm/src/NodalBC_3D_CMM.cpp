#include "NodalBC_3D_CMM.hpp"

NodalBC_3D_CMM::NodalBC_3D_CMM( const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;
  num_per_nodes = 0;

  Create_ID( nFunc );
  
  std::cout<<"===> NodalBC_3D_CMM: No nodal BC is generated. \n";
}


NodalBC_3D_CMM::NodalBC_3D_CMM( const INodalBC * const &nbc_inflow,
    const INodalBC * const &nbc_ring, const int &comp, const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  // regardless of comp, assign all interior inlet nodes as nodal/essential bc
  for(unsigned int ii=0; ii<nbc_inflow->get_num_dir_nodes(); ++ii)
    dir_nodes.push_back( nbc_inflow->get_dir_nodes(ii) );

  // obtain the type of ring nodes' specification
  const std::vector<int> cap_id = nbc_ring -> get_cap_id();

  // Record the number of nodes assigned as essential bc per dof
  std::vector<int> cap_num_dir_nodes( nbc_ring->get_num_caps(), 0 );

  switch( nbc_ring -> get_ring_bc_type() )
  {
    case 0:

      // regardless of comp, all ring nodes are added as essential bc
      for(unsigned int ii = 0; ii < nbc_ring->get_num_dir_nodes(); ++ii)
      {
        dir_nodes.push_back( nbc_ring->get_dir_nodes(ii) );
        cap_num_dir_nodes[ cap_id[ii] ] += 1;
      }
      break;

    case 1:

      // if comp is 0 (rotated x-comp, corresponding to the normal comp),
      // ring node is added to dir_nodes
      for(unsigned int ii = 0; ii < nbc_ring->get_num_dir_nodes(); ++ii)
      {
        if( comp == 0 )
        {
          dir_nodes.push_back( nbc_ring->get_dir_nodes(ii) );
          cap_num_dir_nodes[ cap_id[ii] ] += 1;
        }
      }
      break;
   
    default:

      SYS_T::print_fatal( "NodalBC_3D_CMM Error: there is no such type of essential bc for ring nodes.\n" );
      break;
  }

  // Clean up the dir_nodes and generate ID array
  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  // print data on screen
  std::cout<<"===> NodalBC_3D_CMM specified by \n";
  std::cout<<"     interior of inlet surface"<<std::endl;
  std::cout<<"     outline of inlet surface" << ": " << cap_num_dir_nodes[0] << " nodes" << std::endl;
  for(int ii=1; ii<nbc_ring -> get_num_caps(); ++ii)
    std::cout<<"     outline of outlet surface " << ii-1 << ": " << cap_num_dir_nodes[ii] << " nodes" << std::endl;
  std::cout<<"     is generated. \n";
}


NodalBC_3D_CMM::NodalBC_3D_CMM( const INodalBC * const &nbc_inflow,
        const INodalBC * const &nbc_ring, const INodalBC * const &nbc_wall,
        const int &nFunc )
{
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  dir_nodes.clear();

  // Assign the inlet nodes for this type of nodal/essential bc
  for(unsigned int ii=0; ii<nbc_inflow->get_num_dir_nodes(); ++ii)
    dir_nodes.push_back( nbc_inflow->get_dir_nodes(ii) );

  // Assign the ring nodes for this type of nodal/essential bc
  for(unsigned int ii=0; ii<nbc_ring->get_num_dir_nodes(); ++ii)
    dir_nodes.push_back( nbc_ring->get_dir_nodes(ii) );

  // Assign the wall nodes for this type of nodal/essential bc
  for(unsigned int ii=0; ii<nbc_wall->get_num_dir_nodes(); ++ii)
    dir_nodes.push_back( nbc_wall->get_dir_nodes(ii) );

  VEC_T::sort_unique_resize(dir_nodes);
  
  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  // print data on screen
  std::cout<<"===> NodalBC_3D_CMM specified by \n";
  std::cout<<"     interior of inlet surface, ring nodes, as well as the wall nodes\n";
  std::cout<<"     is generated. \n";
}


NodalBC_3D_CMM::NodalBC_3D_CMM( const INodalBC * const &nbc_wall, const int &nFunc )
{
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  dir_nodes.clear();

  // prepare wall nodes' indices
  std::vector<unsigned int> wall_gnode;
  wall_gnode.clear();
  for(unsigned int ii=0; ii<nbc_wall->get_num_dir_nodes(); ++ii)
    wall_gnode.push_back( nbc_wall->get_dir_nodes(ii) );

  // get wall nodes excluded
  for(unsigned int ii=0; ii< static_cast<unsigned int>(nFunc); ++ii)
  {
    if( !VEC_T::is_invec( wall_gnode, ii) )
      dir_nodes.push_back( ii );
  }

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  // print data on screen
  std::cout<<"===> NodalBC_3D_CMM specified by \n";
  std::cout<<"     the whole domain, with the wall nodes excluded,\n";
  std::cout<<"     is generated. \n";
}

// EOF
