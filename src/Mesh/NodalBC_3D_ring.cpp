#include "NodalBC_3D_ring.hpp"

NodalBC_3D_ring::NodalBC_3D_ring(const int &nFunc)
{
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
  dir_nodes.clear();
  num_dir_nodes = 0;

  Create_ID( nFunc );
  
  num_caps = 0;
  cap_id.clear();
  dominant_n_comp.clear();
  dominant_t_comp.clear();
  outnormal.clear();
  tangential.clear();

  std::cout<<"===> NodalBC_3D_ring::empty is generated. \n";
}


NodalBC_3D_ring::NodalBC_3D_ring( const std::string &inflow_file,
    const std::vector<double> &inflow_outward_vec,
    const std::string &wallfile,
    const std::vector<std::string> &outflow_files,
    const std::vector< std::vector<double> > &outflow_outward_vec,
    const int &nFunc, const int &elemtype )
{
  // No periodic nodes
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  // Aggregate inlet & outlet data
  std::vector<std::string> cap_files = outflow_files;
  cap_files.insert( cap_files.begin(), inflow_file ); 

  num_caps = cap_files.size();
  dominant_n_comp.resize(num_caps);

  outnormal = inflow_outward_vec;
  Vector_3 outvec = Vector_3( inflow_outward_vec[0], inflow_outward_vec[1], inflow_outward_vec[2] );
  dominant_n_comp[0] = outvec.get_dominant_comp();

  for(unsigned int ii=0; ii<outflow_outward_vec.size(); ++ii)
  {
    outnormal.push_back(outflow_outward_vec[ii][0]);
    outnormal.push_back(outflow_outward_vec[ii][1]);
    outnormal.push_back(outflow_outward_vec[ii][2]);

    outvec = Vector_3( outflow_outward_vec[ii][0], outflow_outward_vec[ii][1], outflow_outward_vec[ii][2] );
    dominant_n_comp[ii + 1] = outvec.get_dominant_comp();
  }

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  // Generate the dir-node list with all ring nodes.
  dir_nodes.clear();

  if( elemtype == 501 )
  { 
    SYS_T::file_check(wallfile);

    TET_T::read_vtp_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );

    for(int ii=0; ii<num_caps; ++ii)
    {
      SYS_T::file_check( cap_files[ii] );

      TET_T::read_vtp_grid( cap_files[ii], numpts, numcels, pts, ien, gnode, gelem );
     
      Vector_3 centroid;
      compute_cap_centroid( pts, centroid );

      int num_outline_pts = 0;
      for(unsigned int jj=0; jj<gnode.size(); ++jj)
      {
        if( gnode[jj]<0 ) SYS_T::print_fatal("Error: negative nodal index on cap %d! \n", ii);

        if( VEC_T::is_invec( wall_gnode, gnode[jj]) )
        {
          dir_nodes.push_back( static_cast<unsigned int>( gnode[jj] ) );
          cap_id.push_back( ii );
          compute_tangential( ii, centroid, pts[3*jj], pts[3*jj + 1], pts[3*jj + 2]);
          num_outline_pts += 1;
        }
      }

      // Detect usage of the sv exterior surface (containing caps) as the wall surface
      if(num_outline_pts == numpts)
        SYS_T::print_fatal( "Error: Cap %d has %d outline nodes and %d total nodes. This is likely due to an improper wall mesh.\n", ii, num_outline_pts, numpts );
    }
  }
  else if( elemtype == 502 )
  {
    SYS_T::file_check(wallfile);

    TET_T::read_vtu_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );

    for(int ii=0; ii<num_caps; ++ii)
    {
      SYS_T::file_check( cap_files[ii] );

      TET_T::read_vtu_grid( cap_files[ii], numpts, numcels, pts, ien, gnode, gelem );

      Vector_3 centroid;
      compute_cap_centroid( pts, centroid );

      int num_outline_pts = 0;
      for(unsigned int jj=0; jj<gnode.size(); ++jj)
      {
        if( gnode[jj]<0 ) SYS_T::print_fatal("Error: negative nodal index on cap %d! \n", ii);

        if( VEC_T::is_invec( wall_gnode, gnode[jj]) )
        {
          dir_nodes.push_back( static_cast<unsigned int>( gnode[jj] ) );
          cap_id.push_back( ii );
          compute_tangential( ii, centroid, pts[3*jj], pts[3*jj + 1], pts[3*jj + 2]);
          num_outline_pts += 1;
        }
      }

      // Detect usage of the sv exterior surface (containing caps) as the wall surface
      if(num_outline_pts == numpts)
        SYS_T::print_fatal( "Error: Cap %d has %d outline nodes and %d total nodes. This is likely due to an improper wall mesh.\n", ii, num_outline_pts, numpts );
    }
  }
  else
    SYS_T::print_fatal("Error: Nodal_3D_ring unknown file type.\n");

  num_dir_nodes = dir_nodes.size(); 

  // Generate ID array
  Create_ID( nFunc );

  // Finish and print info on screen
  std::cout << "===> NodalBC_3D_ring specified by \n";
  for(int ii=0; ii<num_caps; ++ii)
    std::cout << "     outline of " << cap_files[ii] << std::endl;
  std::cout<<"     is generated. \n";
}


void NodalBC_3D_ring::compute_cap_centroid( const std::vector<double> &pts, Vector_3 &centroid ) const
{
  const int num_node = static_cast<int>( pts.size() / 3 );

  for(int ii=0; ii<num_node; ++ii)
  {
    centroid(0) += pts[3*ii + 0];
    centroid(1) += pts[3*ii + 1];
    centroid(2) += pts[3*ii + 2];
  }

  centroid(0) /= (double) num_node;
  centroid(1) /= (double) num_node;
  centroid(2) /= (double) num_node;
}


void NodalBC_3D_ring::compute_tangential( const int &cap_id, const Vector_3 &centroid,
    const double &pt_x, const double &pt_y, const double &pt_z )
{
  // Generate radial vector using nodal & centroidal coordinates
  Vector_3 radial_vec = Vector_3( pt_x, pt_y, pt_z );
  radial_vec -= centroid;

  Vector_3 normal_vec = Vector_3( outnormal[3*cap_id], outnormal[3*cap_id+1], outnormal[3*cap_id+2] );
  Vector_3 tan_vec    = cross_product( normal_vec, radial_vec );
  tan_vec.normalize();

  for(int ii=0; ii<3; ++ii) tangential.push_back( tan_vec(ii) );

  // Ensure dominant component indices in the normal and tangential vectors
  // aren't equal 
  if( tan_vec.get_dominant_comp() == dominant_n_comp[ cap_id ] )
  {
    Vector_3 temp( tan_vec );
    temp( tan_vec.get_dominant_comp() ) = 0.0;
    dominant_t_comp.push_back( temp.get_dominant_comp() );
  }
  else
    dominant_t_comp.push_back( tan_vec.get_dominant_comp() );
}

// EOF
