#include "NodalBC_3D_ring.hpp"

NodalBC_3D_ring::NodalBC_3D_ring(
    const std::vector<std::string> &inflow_files,
    const std::vector< Vector_3 > &inlet_outnormal,
    const std::string &wallfile,
    const std::vector<std::string> &outflow_files,
    const std::vector< Vector_3 > &outlet_outnormal,
    const int &nFunc, const int &in_ring_bc_type,
    const int &elemtype ) : ring_bc_type(in_ring_bc_type)
{
  // Aggregate inlet & outlet data
  std::vector<std::string> cap_files = inflow_files;
  for(unsigned int ii=0; ii<outflow_files.size(); ++ii)
    cap_files.push_back( outflow_files[ii] );

  num_caps = cap_files.size();
  Q.resize(9 * num_caps);

  outnormal.clear();

  for(unsigned int ii=0; ii<inlet_outnormal.size(); ++ii)
  {
    outnormal.push_back( inlet_outnormal[ii](0) );
    outnormal.push_back( inlet_outnormal[ii](1) );
    outnormal.push_back( inlet_outnormal[ii](2) );
  }

  for(unsigned int ii=0; ii<outlet_outnormal.size(); ++ii)
  {
    outnormal.push_back( outlet_outnormal[ii](0) );
    outnormal.push_back( outlet_outnormal[ii](1) );
    outnormal.push_back( outlet_outnormal[ii](2) );
  }

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  // Generate the dir-node list with all ring nodes.
  dir_nodes.clear();
  cap_id.clear();

  if( elemtype == 501 )
  { 
    SYS_T::file_check(wallfile);

    TET_T::read_vtp_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );

    for(int ii=0; ii<num_caps; ++ii)
    {
      SYS_T::file_check( cap_files[ii] );

      TET_T::read_vtp_grid( cap_files[ii], numpts, numcels, pts, ien, gnode, gelem );
     
      const Vector_3 centroid = compute_cap_centroid( pts );

      int num_outline_pts = 0;
      for(unsigned int jj=0; jj<gnode.size(); ++jj)
      {
        if( gnode[jj]<0 ) SYS_T::print_fatal("Error: negative nodal index on cap %d! \n", ii);

        if( VEC_T::is_invec( wall_gnode, gnode[jj]) )
        {
          dir_nodes.push_back( static_cast<unsigned int>( gnode[jj] ) );
          cap_id.push_back( ii );

          // Compute the cap's skew bc transformation matrix with the first ring node
          if( num_outline_pts == 0 )
          {
            Vector_3 radial_vec = Vector_3( pts[3*jj], pts[3*jj + 1], pts[3*jj + 2] );
            radial_vec -= centroid;
            radial_vec.normalize();

            Vector_3 normal_vec = Vector_3( outnormal[3*ii], outnormal[3*ii+1], outnormal[3*ii+2] );
            Vector_3 tan_vec    = cross_product( normal_vec, radial_vec );
            tan_vec.normalize();

            Q[9*ii+0] = normal_vec.x(); Q[9*ii+1] = normal_vec.y(); Q[9*ii+2] = normal_vec.z();
            Q[9*ii+3] = radial_vec.x(); Q[9*ii+4] = radial_vec.y(); Q[9*ii+5] = radial_vec.z();
            Q[9*ii+6] = tan_vec.x();    Q[9*ii+7] = tan_vec.y();    Q[9*ii+8] = tan_vec.z();
          }

          num_outline_pts += 1;
        }
      }

      // Detect usage of the sv exterior surface (containing caps) as the wall surface
      if( num_outline_pts == numpts )
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

      const Vector_3 centroid = compute_cap_centroid( pts );

      int num_outline_pts = 0;
      for(unsigned int jj=0; jj<gnode.size(); ++jj)
      {
        if( gnode[jj]<0 ) SYS_T::print_fatal("Error: negative nodal index on cap %d! \n", ii);

        if( VEC_T::is_invec( wall_gnode, gnode[jj]) )
        {
          dir_nodes.push_back( static_cast<unsigned int>( gnode[jj] ) );
          cap_id.push_back( ii );

          // Compute the cap's skew bc transformation matrix with the first ring node
          if( num_outline_pts == 0 )
          {
            Vector_3 radial_vec = Vector_3( pts[3*jj], pts[3*jj + 1], pts[3*jj + 2] );
            radial_vec -= centroid;

            Vector_3 normal_vec = Vector_3( outnormal[3*ii], outnormal[3*ii+1], outnormal[3*ii+2] );
            Vector_3 tan_vec    = cross_product( normal_vec, radial_vec );
            tan_vec.normalize();

            Q[9*ii+0] = normal_vec.x(); Q[9*ii+1] = normal_vec.y(); Q[9*ii+2] = normal_vec.z();
            Q[9*ii+3] = radial_vec.x(); Q[9*ii+4] = radial_vec.y(); Q[9*ii+5] = radial_vec.z();
            Q[9*ii+6] = tan_vec.x();    Q[9*ii+7] = tan_vec.y();    Q[9*ii+8] = tan_vec.z();
          }

          num_outline_pts += 1;
        }
      }

      // Detect usage of the sv exterior surface (containing caps) as the wall surface
      if( num_outline_pts == numpts )
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
  std::cout<<"     is generated ";
  if(ring_bc_type == 0) std::cout<<"for fully clamped case (ring_bc_type = 0).\n";
  else if(ring_bc_type == 1) std::cout<<"for in-plane motion (ring_bc_type = 1).\n";
  else SYS_T::print_fatal("Error: NodalBC_3D_ring does not allow this ring_bc_type!\n");
}

Vector_3 NodalBC_3D_ring::compute_cap_centroid( const std::vector<double> &pts ) const
{
  const int num_node = static_cast<int>( pts.size() / 3 );

  Vector_3 centroid(0.0, 0.0, 0.0);

  for(int ii=0; ii<num_node; ++ii)
  {
    centroid.x() += pts[3*ii + 0];
    centroid.y() += pts[3*ii + 1];
    centroid.z() += pts[3*ii + 2];
  }

  centroid *= (1.0 / static_cast<double>( num_node ));

  return centroid;
}

// EOF