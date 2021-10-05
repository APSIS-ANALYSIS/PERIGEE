#include "NodalBC_3D_wall.hpp"

NodalBC_3D_wall::NodalBC_3D_wall(
    const std::vector<std::string> &inflow_files,
    const std::string &wall_file,
    const std::vector<std::string> &outflow_files,
    const int &nFunc, const int &elemtype )
: num_nbc(0)
{
  // No periodic nodes
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes.clear();

  dir_nodes.resize(     num_nbc );
  num_dir_nodes.resize( num_nbc );

  // Aggregate inlet and outlet data
  std::vector<std::string> cap_files = inflow_files;
  for( unsigned int ii=0; ii<outflow_files.size(); ++ii )
    cap_files.push_back( outflow_files[ii] );

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
    SYS_T::print_fatal("Error: Nodal_3D_wall unknown file type.\n");

  VEC_T::sort_unique_resize( ring_gnode );

  const int nbc_id = 0;
  dir_nodes[nbc_id].clear();

  // exclude the ring nodes
  for(unsigned int ii=0; ii<wall_gnode.size(); ++ii)
  {
    if( !VEC_T::is_invec( ring_gnode, wall_gnode[ii] ) )
      dir_nodes[nbc_id].push_back( wall_gnode[ii] );
  }

  num_dir_nodes[nbc_id] = dir_nodes[nbc_id].size();

  // Generate the ID array
  Create_ID( nFunc );

  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_wall specified by "<<wall_file;
  std::cout<<", with ring nodes on\n";
  for( const auto  &capfile : cap_files )
    std::cout << "     on the outline of " << capfile << std::endl;
  std::cout<<"     excluded, is generated.\n";
}

// EOF
