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

  // Aggregate inlet and outlet data
  std::vector<std::string> cap_files = outflow_files;
  cap_files.insert( cap_files.begin(), inflow_file );

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  if( elemtype == 501 )
  {
    SYS_T::file_check(wallfile);

    TET_T::read_vtp_grid( wallfile, wall_numpts, wall_numcels, wall_pts,
        wall_ien, wall_gnode, wall_gelem );

    std::vector<int> ring_gnode {};

    for(unsigned int ii=0; ii<cap_files.size(); ++ii)
    {
      SYS_T::file_check( cap_files[ii] );

      TET_T::read_vtp_grid( cap_files[ii], numpts, numcels, pts, ien, gnode, gelem );
    }

  }
  else if( elemtype == 502 )
  {
  }
  else
    SYS_T::print_fatal("Error: Nodal_3D_ring unknown file type.\n");


}


NodalBC_3D_wall::~NodalBC_3D_wall()
{}

// EOF
