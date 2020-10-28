#include "NodalBC_3D_vtp.hpp"

NodalBC_3D_vtp::NodalBC_3D_vtp( const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;
  num_per_nodes = 0;

  Create_ID( nFunc );
  
  std::cout<<"===> NodalBC_3D_vtp: No nodal BC is generated. \n";
}


NodalBC_3D_vtp::NodalBC_3D_vtp( const std::string &inflow_vtp_file,
    const std::string &wall_vtp_file,
    const std::vector<std::string> &outflow_vtp_files,
    const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  // Assign all inlet nodes
  SYS_T::file_check( inflow_vtp_file );

  TET_T::read_vtp_grid( inflow_vtp_file, numpts, numcels, pts, ien, gnode, gelem );

  if( numpts != static_cast<int>(gnode.size()) )
    SYS_T::print_fatal("Error: numpts != global_node.size() for the inlet! \n");

  for(unsigned int ii=0; ii<gnode.size(); ++ii)
  {
    if(gnode[ii]<0) SYS_T::print_fatal("Error: negative nodal index on the inlet! \n");

    dir_nodes.push_back( static_cast<unsigned int>( gnode[ii] ) );
  }

  // Assign only outline (ring) nodes on each outlet
  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  SYS_T::file_check( wall_vtp_file );

  TET_T::read_vtp_grid( wall_vtp_file, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );

  if( wall_numpts != static_cast<int>(wall_gnode.size()) )
      SYS_T::print_fatal("Error: numpts != global_node.size() for the wall! \n");

  for(unsigned int ii=0; ii<wall_gnode.size(); ++ii)
    if(wall_gnode[ii]<0) SYS_T::print_fatal("Error: negative nodal index on the wall! \n");

  const unsigned int num_outlets = outflow_vtp_files.size();
  for(unsigned int ii=0; ii<num_outlets; ++ii)
  {
    SYS_T::file_check( outflow_vtp_files[ii] );

    TET_T::read_vtp_grid( outflow_vtp_files[ii], numpts, numcels, pts, ien, gnode, gelem );

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

  std::cout<<"===> NodalBC_3D_vtp specified by \n";
  std::cout<<"     "<<inflow_vtp_file<<std::endl;
  for(unsigned int ii=0; ii<num_outlets; ++ii)
    std::cout<<"     outline of "<<outflow_vtp_files[ii]<<std::endl;
  std::cout<<"     is generated. \n";

}


NodalBC_3D_vtp::NodalBC_3D_vtp( const std::string &vtpfileName,
    const int &nFunc )
{
  SYS_T::file_check( vtpfileName );

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  TET_T::read_vtp_grid( vtpfileName, numpts, numcels, pts, ien, gnode, gelem );

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

  std::cout<<"===> NodalBC_3D_vtp specified by "<<vtpfileName<<" is generated. \n";
}


NodalBC_3D_vtp::NodalBC_3D_vtp( const std::vector<std::string> &vtpfileList, 
    const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  const unsigned int num_file = vtpfileList.size();

  for(unsigned int ii=0; ii<num_file; ++ii)
  {
    SYS_T::file_check( vtpfileList[ii] );

    int numpts, numcels;
    std::vector<double> pts;
    std::vector<int> ien, gnode, gelem;

    TET_T::read_vtp_grid( vtpfileList[ii], numpts, numcels, pts, ien, gnode, gelem );

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

  std::cout<<"===> NodalBC_3D_vtp specified by \n";
  for(unsigned int ii=0; ii<num_file; ++ii)
    std::cout<<"     "<<vtpfileList[ii]<<std::endl;
  std::cout<<"     is generated. \n";
}


NodalBC_3D_vtp::~NodalBC_3D_vtp()
{}

// EOF
