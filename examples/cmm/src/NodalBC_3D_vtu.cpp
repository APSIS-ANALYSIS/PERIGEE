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


NodalBC_3D_vtu::NodalBC_3D_vtu( const std::string &inflow_vtu_file,
    const std::string &wall_vtu_file,
    const std::vector<std::string> &outflow_vtu_files,
    const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  // ==== WOMERSLEY CHANGES BEGIN ====
  // // Assign all inlet nodes
  // SYS_T::file_check( inflow_vtu_file );

  // TET_T::read_vtu_grid( inflow_vtu_file, numpts, numcels, pts, ien, gnode, gelem );

  // if( numpts != static_cast<int>(gnode.size()) )
  //   SYS_T::print_fatal("Error: numpts != global_node.size() for the inlet! \n");

  // for(unsigned int ii=0; ii<gnode.size(); ++ii)
  // {
  //   if(gnode[ii]<0) SYS_T::print_fatal("Error: negative nodal index on the inlet! \n");

  //   dir_nodes.push_back( static_cast<unsigned int>( gnode[ii] ) );
  // }
  // ==== WOMERSLEY CHANGES END ====

  // Assign only outline (ring) nodes on each outlet
  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  SYS_T::file_check( wall_vtu_file );

  TET_T::read_vtu_grid( wall_vtu_file, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );

  if( wall_numpts != static_cast<int>(wall_gnode.size()) )
      SYS_T::print_fatal("Error: numpts != global_node.size() for the wall! \n");

  for(unsigned int ii=0; ii<wall_gnode.size(); ++ii)
    if(wall_gnode[ii]<0) SYS_T::print_fatal("Error: negative nodal index on the wall! \n");

  const unsigned int num_outlets = outflow_vtu_files.size();
  for(unsigned int ii=0; ii<num_outlets; ++ii)
  {
    SYS_T::file_check( outflow_vtu_files[ii] );

    TET_T::read_vtu_grid( outflow_vtu_files[ii], numpts, numcels, pts, ien, gnode, gelem );

    if( numpts != static_cast<int>(gnode.size()) )
      SYS_T::print_fatal("Error: numpts != global_node.size() for outlet %d! \n", ii);

    int num_outline_pts = 0;
    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      if(gnode[jj]<0) SYS_T::print_fatal("Error: negative nodal index on outlet %d! \n", ii);

      if( VEC_T::is_invec( wall_gnode, gnode[jj]) )
      {
        dir_nodes.push_back( static_cast<unsigned int>( gnode[jj] ) );
        num_outline_pts += 1;
      }

      // Detect usage of the sv exterior surface (containing caps) as the wall surface
      if(num_outline_pts == numpts)
        SYS_T::print_fatal( "Error: Outlet %d has %d outline nodes and %d total nodes. This is likely due to an improper wall mesh.\n", ii, num_outline_pts, numpts );
    }

  } // end loop over outlets

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtu specified by \n";
  // ==== WOMERSLEY CHANGES BEGIN ====
  // std::cout<<"     "<<inflow_vtu_file<<std::endl;
  // ==== WOMERSLEY CHANGES END ====
  for(unsigned int ii=0; ii<num_outlets; ++ii)
    std::cout<<"     outline of "<<outflow_vtu_files[ii]<<std::endl;
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

  if( numpts != static_cast<int>(gnode.size()) )
    SYS_T::print_fatal("Error: the numpts != global_node.size()! \n");

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

    if( numpts != static_cast<int>(gnode.size()) )
      SYS_T::print_fatal("Error: the numpts != global_node.size()! \n");

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
