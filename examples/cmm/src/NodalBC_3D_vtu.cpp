#include "NodalBC_3D_vtu.hpp"

NodalBC_3D_vtu::NodalBC_3D_vtu( const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;
  num_per_nodes = 0;

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtu: No nodal BC is generated. \n";
}


NodalBC_3D_vtu::NodalBC_3D_vtu( const INodalBC * const &nbc_inflow,
    const INodalBC * const &nbc_ring,
    const int &comp, const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  // Assign the inlet nodes for this type of nodal/essential bc
  for(unsigned int ii=0; ii<nbc_inflow->get_num_dir_nodes(); ++ii)
    dir_nodes.push_back( nbc_inflow->get_dir_nodes(ii) );

  // obtain the type of ring nodes' specification
  const int ringbc_type = nbc_ring -> get_ring_bc_type();
  const std::vector<int> cap_id = nbc_ring -> get_cap_id();
  const std::vector<int> dom_n_comp = nbc_ring -> get_dominant_n_comp();
  const std::vector<int> dom_t_comp = nbc_ring -> get_dominant_t_comp();

  std::vector<int> num_dom_n_pts( nbc_ring->get_num_caps(), 0 );
  std::vector<int> num_dom_t_pts( nbc_ring->get_num_caps(), 0 );

  switch( ringbc_type )
  {
    case 0:
      // ring nodes are added as essential bc
      for(unsigned int ii=0; ii<nbc_ring->get_num_dir_nodes(); ++ii)
      {
        dir_nodes.push_back( nbc_ring -> get_dir_nodes(ii) );
        num_dom_n_pts[ cap_id[ii] ] += 1;
        num_dom_t_pts[ cap_id[ii] ] += 1;
      }
      break;

    case 1:
      // if dom_n_comp equals comp, ring node is added to dir_nodes
      for(unsigned int ii=0; ii<nbc_ring->get_num_dir_nodes(); ++ii)
      {
        if( comp == dom_n_comp[ cap_id[ii] ] )
        {
          dir_nodes.push_back( nbc_ring -> get_dir_nodes(ii) );
          num_dom_n_pts[ cap_id[ii] ] += 1;
        }
      }
      break;

    case 2:
      // if dom_n_comp or dom_t_comp equals comp, ring node is added to dir_nodes
      for(unsigned int ii=0; ii<nbc_ring->get_num_dir_nodes(); ++ii)
      {
        if( comp == dom_n_comp[ cap_id[ii] ] )
        {
          dir_nodes.push_back( nbc_ring -> get_dir_nodes(ii) );
          num_dom_n_pts[ cap_id[ii] ] += 1;
        }
        else if( comp == dom_t_comp[ii] )
        {
          dir_nodes.push_back( nbc_ring -> get_dir_nodes(ii) );
          num_dom_t_pts[ cap_id[ii] ] += 1;
        }
      }
      break;

    default:
      SYS_T::print_fatal("Error: there is no such type of essential bc for ring nodes.\n");
      break;
  }

  // Clean up the dir_nodes and generate ID array
  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  // print data on screen
  std::cout<<"===> NodalBC_3D_vtu specified by \n";
  std::cout<<"     interior of inlet surface"<<std::endl;
  std::cout<<"     outline of inlet surface"<<": "<<num_dom_n_pts[0]<<" dom_n nodes, "<< num_dom_t_pts[0]<<" dom_t nodes"<<std::endl;
  for(int ii=1; ii<nbc_ring -> get_num_caps(); ++ii)
    std::cout<<"     outline of outlet surface "<<ii-1<<": "<<num_dom_n_pts[ii]<<" dom_n nodes, "<< num_dom_t_pts[ii]<<" dom_t nodes"<<std::endl;
  std::cout<<"     is generated. \n";
}


NodalBC_3D_vtu::NodalBC_3D_vtu( const std::string &vtufilename, 
    const int &nFunc )
{
  SYS_T::file_check( vtufilename );

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  TET_T::read_vtu_grid( vtufilename, numpts, numcels, pts, ien, gnode, gelem );

  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
  num_dir_nodes = numpts;

  dir_nodes.resize( gnode.size() );
  for(unsigned int ii=0; ii<gnode.size(); ++ii)
  {
    if(gnode[ii]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");

    dir_nodes[ii] = static_cast<unsigned int>( gnode[ii] ); 
  }

  // Generate the ID array
  Create_ID(nFunc);

  std::cout<<"===> NodalBC_3D_vtu specified by "<<vtufilename<<" is generated. \n";
}


NodalBC_3D_vtu::NodalBC_3D_vtu( const std::vector<std::string> &vtufileList,
    const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  const unsigned int num_file = vtufileList.size();

  for(unsigned int ii=0; ii<num_file; ++ii)
  {
    SYS_T::file_check( vtufileList[ii] );

    int numpts, numcels;
    std::vector<double> pts;
    std::vector<int> ien, gnode, gelem;

    TET_T::read_vtu_grid( vtufileList[ii], numpts, numcels, pts, ien, gnode, gelem );

    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      if(gnode[jj]<0) SYS_T::print_fatal("Error: there are negative nodal index! \n");

      dir_nodes.push_back( static_cast<unsigned int>( gnode[jj]) );
    }
  }

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtu specified by \n";
  for(unsigned int ii=0; ii<num_file; ++ii)
    std::cout<<"     "<<vtufileList[ii]<<std::endl;
  std::cout<<"     is generated. \n";
}


NodalBC_3D_vtu::~NodalBC_3D_vtu()
{}

// EOF
