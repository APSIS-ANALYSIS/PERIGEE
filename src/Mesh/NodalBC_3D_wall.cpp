#include "NodalBC_3D_wall.hpp"

NodalBC_3D_wall::NodalBC_3D_wall(
    const std::vector<std::string> &inflow_files,
    const std::string &wall_file,
    const std::vector<std::string> &outflow_files,
    const int &nFunc, const int &elemtype )
{
  // No periodic nodes
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  // Aggregate inlet and outlet data
  std::vector<std::string> cap_files = inflow_files;
  for( unsigned int ii=0; ii<outflow_files.size(); ++ii )
    cap_files.push_back( outflow_files[ii] );

  std::vector<int> ring_gnode {};

  for( const auto &capfile : cap_files )
  {
    SYS_T::file_check( capfile );

    const auto gnode = VTK_T::read_int_PointData(capfile, "GlobalNodeID"); 

    VEC_T::insert_end( ring_gnode, gnode );
  }

  const auto wall_gnode = VTK_T::read_int_PointData(wall_file, "GlobalNodeID");

  VEC_T::sort_unique_resize( ring_gnode );

  dir_nodes.clear();

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

// EOF
