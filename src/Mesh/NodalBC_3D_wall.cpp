#include "NodalBC_3D_wall.hpp"

NodalBC_3D_wall::NodalBC_3D_wall( const std::string &inflow_file,
    const std::string &wall_file,
    const std::vector<std::string> &outflow_files,
    const int &nFunc, const int &elemtype )
{
  // No periodic nodes
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  dir_nodes.clear();

  // Aggregate inlet and outlet data
  std::vector<std::string> cap_files = outflow_files;
  cap_files.insert( cap_files.begin(), inflow_file );

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  std::vector<int> ring_gnode {};

  if( elemtype == 501 )
  {
    SYS_T::file_check(wall_file);

    TET_T::read_vtp_grid( wall_file, wall_numpts, wall_numcels, wall_pts,
        wall_ien, wall_gnode, wall_gelem );

    for( const auto &capfile : cap_files )
    {
      SYS_T::file_check( capfile );

      TET_T::read_vtp_grid( capfile, numpts, numcels, pts, ien, gnode, gelem );
    
      VEC_T::insert_end( ring_gnode, gnode );
    }
  }
  else if( elemtype == 502 )
  {
    SYS_T::file_check(wall_file);

    TET_T::read_vtu_grid( wall_file, wall_numpts, wall_numcels, wall_pts,
        wall_ien, wall_gnode, wall_gelem );

    for( const auto &capfile : cap_files )
    {
      SYS_T::file_check( capfile );

      TET_T::read_vtu_grid( capfile, numpts, numcels, pts, ien, gnode, gelem );

      VEC_T::insert_end( ring_gnode, gnode );
    }
  }
  else
    SYS_T::print_fatal("Error: Nodal_3D_ring unknown file type.\n");

  VEC_T::sort_unique_resize( ring_gnode );

  // exclude the ring nodes
  for(unsigned int ii=0; ii<wall_gnode.size(); ++ii)
  {
    if( !VEC_T::is_invec( ring_gnode, wall_gnode[ii] ) )
      dir_nodes.push_back( wall_gnode[ii] );
  }

  num_dir_nodes = dir_nodes.size();

  // Generate the ID array
  Create_ID( nFunc );

  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_wall specified by "<<wall_file;
  std::cout<<", with ring nodes on\n";
  for( const auto  &capfile : cap_files )
    std::cout << "     on the outline of " << capfile << std::endl;
  std::cout<<"     excluded, is generated.\n";
}


NodalBC_3D_wall::~NodalBC_3D_wall()
{}

// EOF
