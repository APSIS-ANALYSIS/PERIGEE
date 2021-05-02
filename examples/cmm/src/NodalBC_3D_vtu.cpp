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

  // Assign all inlet nodes
  SYS_T::file_check( inflow_vtu_file );

  TET_T::read_vtu_grid( inflow_vtu_file, numpts, numcels, pts, ien, gnode, gelem );

  for(unsigned int ii=0; ii<gnode.size(); ++ii)
  {
    if(gnode[ii]<0) SYS_T::print_fatal("Error: negative nodal index on the inlet! \n");

    dir_nodes.push_back( static_cast<unsigned int>( gnode[ii] ) );
  }

  // Assign only outline (ring) nodes on each outlet
  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  SYS_T::file_check( wall_vtu_file );

  TET_T::read_vtu_grid( wall_vtu_file, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );

  for(unsigned int ii=0; ii<wall_gnode.size(); ++ii)
    if(wall_gnode[ii]<0) SYS_T::print_fatal("Error: negative nodal index on the wall! \n");

  const unsigned int num_outlets = outflow_vtu_files.size();
  for(unsigned int ii=0; ii<num_outlets; ++ii)
  {
    SYS_T::file_check( outflow_vtu_files[ii] );

    TET_T::read_vtu_grid( outflow_vtu_files[ii], numpts, numcels, pts, ien, gnode, gelem );

    int num_outline_pts = 0;
    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      if(gnode[jj]<0) SYS_T::print_fatal("Error: negative nodal index on outlet %d! \n", ii);

      if( VEC_T::is_invec( wall_gnode, gnode[jj]) )
      {
        dir_nodes.push_back( static_cast<unsigned int>( gnode[jj] ) );
        num_outline_pts += 1;
      }
    }

    // Detect usage of the sv exterior surface (containing caps) as the wall surface
    if(num_outline_pts == numpts)
      SYS_T::print_fatal( "Error: Outlet %d has %d outline nodes and %d total nodes. This is likely due to an improper wall mesh.\n", ii, num_outline_pts, numpts );

  } // end loop over outlets

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtu specified by \n";
  std::cout<<"     "<<inflow_vtu_file<<std::endl;
  for(unsigned int ii=0; ii<num_outlets; ++ii)
    std::cout<<"     outline of "<<outflow_vtu_files[ii]<<std::endl;
  std::cout<<"     is generated. \n";
}


NodalBC_3D_vtu::NodalBC_3D_vtu( const std::string &inflow_vtu_file,
        const std::vector<double> &inflow_outward_vec,
        const std::string &wall_vtu_file,
        const std::vector<std::string> &outflow_vtu_files,
        const std::vector< std::vector<double> > &outflow_outward_vec,
        const int &type, const int &comp, const int &nFunc )
{
  if( type != 0 && type != 1 )
    SYS_T::print_fatal("Error: NodalBC_3D_vtu: No such nodal BC type.\n");

  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  const unsigned int num_outlets = outflow_vtu_files.size();

  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  SYS_T::file_check( wall_vtu_file );

  TET_T::read_vtu_grid( wall_vtu_file, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );

  for(unsigned int ii=0; ii<wall_gnode.size(); ++ii)
    if(wall_gnode[ii]<0) SYS_T::print_fatal("Error: negative nodal index on the wall! \n");

  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;
  Vector_3 outvec, centroid;
  int dom_n_comp, dom_t_comp;

  // Inlet: Assign all interior nodes
  SYS_T::file_check( inflow_vtu_file );

  TET_T::read_vtu_grid( inflow_vtu_file, numpts, numcels, pts, ien, gnode, gelem );

  outvec = Vector_3( inflow_outward_vec[0], inflow_outward_vec[1], inflow_outward_vec[2] );
  dom_n_comp = outvec.get_dominant_comp();

  compute_cap_centroid( pts, centroid );

  std::vector<int> num_dom_n_pts, num_dom_t_pts;
  num_dom_n_pts.resize( 1 + num_outlets ); num_dom_t_pts.resize( 1 + num_outlets );

  int num_outline_pts = 0;
  for(unsigned int ii=0; ii<gnode.size(); ++ii)
  {
    if(gnode[ii]<0) SYS_T::print_fatal("Error: negative nodal index on the inlet! \n");

    if( VEC_T::is_invec( wall_gnode, gnode[ii]) )
    {
      // Assign dominant component index of unit normal vector
      if( dom_n_comp == comp )
      {
        dir_nodes.push_back( static_cast<unsigned int>( gnode[ii] ) );
        num_dom_n_pts[0] += 1;
      }

      // For purely radial motion, assign dominant component index of unit tangential vector.
      // This is guaranteed by compute_tangential() to be a different index than dom_n_comp. 
      else if( type == 1 )
      {
        dom_t_comp = compute_tangential( outvec, centroid, pts[3*ii], pts[3*ii + 1], pts[3*ii + 2] );
        if( dom_t_comp == comp )
        {
          dir_nodes.push_back( static_cast<unsigned int>( gnode[ii] ) );
          num_dom_t_pts[0] += 1;
        }
      }

      num_outline_pts += 1;
    }
    else
      dir_nodes.push_back( static_cast<unsigned int>( gnode[ii] ) );
  }

  // Detect usage of the sv exterior surface (containing caps) as the wall surface
  if(num_outline_pts == numpts)
    SYS_T::print_fatal( "Error: Inlet has %d outline nodes and %d total nodes. This is likely due to an improper wall mesh.\n", num_outline_pts, numpts );

  // Outlets
  for(unsigned int ii=0; ii<num_outlets; ++ii)
  {
    SYS_T::file_check( outflow_vtu_files[ii] );

    TET_T::read_vtu_grid( outflow_vtu_files[ii], numpts, numcels, pts, ien, gnode, gelem );

    outvec = Vector_3( outflow_outward_vec[ii][0], outflow_outward_vec[ii][1], outflow_outward_vec[ii][2] );
    dom_n_comp = outvec.get_dominant_comp();

    compute_cap_centroid( pts, centroid );

    num_outline_pts = 0;
    for(unsigned int jj=0; jj<gnode.size(); ++jj)
    {
      if(gnode[jj]<0) SYS_T::print_fatal("Error: negative nodal index on outlet %d! \n", ii);

      if( VEC_T::is_invec( wall_gnode, gnode[jj]) )
      {
        // Assign dominant component index of unit normal vector
        if( dom_n_comp == comp )
        {
          dir_nodes.push_back( static_cast<unsigned int>( gnode[jj] ) );
          num_dom_n_pts[1+ii] += 1;
        }

        // For purely radial motion, assign dominant component index of unit tangential vector.
        // This is guaranteed by compute_tangential() to be a different index than dom_n_comp. 
        else if( type == 1 )
        {
          dom_t_comp = compute_tangential( outvec, centroid, pts[3*jj], pts[3*jj + 1], pts[3*jj + 2] );
          if( dom_t_comp == comp )
          {
            dir_nodes.push_back( static_cast<unsigned int>( gnode[jj] ) );
            num_dom_t_pts[1+ii] += 1;
          }
        }

        num_outline_pts += 1;
      }
    }

    // Detect usage of the sv exterior surface (containing caps) as the wall surface
    SYS_T::print_fatal_if( num_outline_pts == numpts, "Error: Outlet %d has %d outline nodes and %d total nodes. This is likely due to an improper wall mesh.\n",
        ii, num_outline_pts, numpts );

  } // end loop over outlets

  VEC_T::sort_unique_resize(dir_nodes);

  num_dir_nodes = dir_nodes.size();

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_vtu specified by \n";
  std::cout<<"     interior of "<<inflow_vtu_file<<std::endl;
  std::cout<<"     outline of "<<inflow_vtu_file<<": "<<num_dom_n_pts[0]<<" dom_n nodes, "<< num_dom_t_pts[0]<<" dom_t nodes"<<std::endl;
  for(unsigned int ii=0; ii<num_outlets; ++ii)
    std::cout<<"     outline of "<<outflow_vtu_files[ii+1]<<": "<<num_dom_n_pts[ii+1]<<" dom_n nodes, "<< num_dom_t_pts[ii+1]<<" dom_t nodes"<<std::endl;
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


void NodalBC_3D_vtu::compute_cap_centroid( const std::vector<double> &pts, Vector_3 &centroid ) const
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


int NodalBC_3D_vtu::compute_tangential( const Vector_3 &outvec, const Vector_3 &centroid,
    const int &pt_x, const int &pt_y, const int &pt_z )
{
  // Generate radial vector using nodal & centroidal coordinates
  Vector_3 radial_vec = Vector_3( pt_x, pt_y, pt_z );
  radial_vec -= centroid;

  Vector_3 tan_vec = cross_product( outvec, radial_vec );
  tan_vec.normalize();

  // Ensure dominant component indices in the normal and tangential vectors
  // aren't equal 
  if( tan_vec.get_dominant_comp() == outvec.get_dominant_comp() )
  {
    Vector_3 temp( tan_vec );
    temp( tan_vec.get_dominant_comp() ) = 0.0;
    return temp.get_dominant_comp();
  }
  else
    return tan_vec.get_dominant_comp();
}

// EOF
