#include "NodalBC_3D_wall.hpp"

NodalBC_3D_wall::NodalBC_3D_wall( const std::string &inflow_file,
    const std::string &wall_file,
    const std::vector<std::string> &outflow_files,
    const int &nFunc, const int &elemtype = 501 )
{
  // No periodic nodes
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  dir_nodes.clear();




}


NodalBC_3D_wall::~NodalBC_3D_wall()
{}

// EOF
