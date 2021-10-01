#include "NodalBC_3D_inflow.hpp"

NodalBC_3D_inflow::NodalBC_3D_inflow(const int &nFunc)
: num_nbc( 0 )
{
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes.clear();

  dir_nodes.clear();
  num_dir_nodes.clear();

  Create_ID( nFunc );
  
  inf_active_area.clear();

  centroid.clear();
  outnormal.clear();
  
  num_out_bc_pts.clear();
  outline_pts.clear();
  intNA.clear();

  std::cout<<"===> NodalBC_3D_inflow::empty is generated. \n";
}


NodalBC_3D_inflow::NodalBC_3D_inflow( const std::string &inffile,
    const std::string &wallfile, const int &nFunc,
    const Vector_3 &in_outnormal,
    const int &elemtype )
: num_nbc( 1 )
{
  SYS_T::file_check(inffile);
  SYS_T::file_check(wallfile);

  // 1. No periodic nodes
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes.clear();

  // 2. Analyze the file type and read in the data
  num_node.resize( num_nbc );
  num_cell.resize( num_nbc );
  nLocBas.resize(  num_nbc );

  tri_ien.resize(num_nbc );
  pt_xyz.resize( num_nbc );

  global_node.resize(  num_nbc );
  global_cell.resize(  num_nbc );

  dir_nodes.resize(     num_nbc );
  num_dir_nodes.resize( num_nbc );

  centroid.resize(  num_nbc );
  outnormal.resize( num_nbc );

  outline_pts.resize(    num_nbc );
  num_out_bc_pts.resize( num_nbc );

  inf_active_area.resize( num_nbc );
  face_area.resize(       num_nbc );

  intNA.resize( num_nbc ); 

  // Read the files
  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  if( elemtype == 501 )
  { 
    nLocBas[0] = 3;

    TET_T::read_vtp_grid( inffile, num_node[0], num_cell[0], pt_xyz[0], tri_ien[0], global_node[0], global_cell[0] );

    TET_T::read_vtp_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );
  }
  else if( elemtype == 502 )
  {
    nLocBas[0] = 6;

    TET_T::read_vtu_grid( inffile, num_node[0], num_cell[0], pt_xyz[0], tri_ien[0], global_node[0], global_cell[0] );

    TET_T::read_vtu_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );
  }
  else SYS_T::print_fatal("Error: unknown element type.\n");

  // Generate the dir-node list. Nodes belonging to the wall are excluded.
  dir_nodes[0].clear();
  
  for(unsigned int ii=0; ii<global_node[0].size(); ++ii)
  {
    SYS_T::print_fatal_if( global_node[0][ii]<0, "Error: negative nodal index! \n");

    if( !VEC_T::is_invec( wall_gnode, global_node[0][ii]) )
      dir_nodes[0].push_back( global_node[0][ii] );
  }

  num_dir_nodes[0] = dir_nodes[0].size(); 

  // Generate ID array
  Create_ID( nFunc );

  // Calculate the centroid of the surface
  centroid[0].gen_zero();
  for(int ii=0; ii<num_node[0]; ++ii)
  {
    centroid[0](0) += pt_xyz[0][3*ii+0];
    centroid[0](1) += pt_xyz[0][3*ii+1];
    centroid[0](2) += pt_xyz[0][3*ii+2];
  }
  centroid[0].scale( 1.0 / (double) num_node[0] );

  // assign outward normal vector from the input
  outnormal[0] = in_outnormal;

  // Collect nodes that belong to the wall, and set up a vector that
  // is 1 on the interior nodes and 0 on the wall bc nodes.
  outline_pts[0].clear();

  num_out_bc_pts[0] = 0;

  double * temp_sol = new double [num_node[0]];  

  for(int ii=0; ii<num_node[0]; ++ii)
  {
    // If the node is not on the wall, it is an interior node, so set
    // the element to 1.
    if( !VEC_T::is_invec(wall_gnode, global_node[0][ii]) ) 
      temp_sol[ii] = 1.0;

    // otherwise, the node is on the wall surface, so set element to 0.
    else 
    {
      temp_sol[ii] = 0.0;

      // Also store the point's coordinates in outline points
      num_out_bc_pts[0] += 1;
      outline_pts[0].push_back( pt_xyz[0][3*ii+0] );
      outline_pts[0].push_back( pt_xyz[0][3*ii+1] );
      outline_pts[0].push_back( pt_xyz[0][3*ii+2] );
    }
  }

  // If the number of surface nodes matches num_out_bc_pts, the wall
  // mesh contains the inlet surface. This is a common error when
  // adopting sv files, where the user uses the combined exterior surface
  // as the wall mesh. We will throw an error message if detected.
  if( num_out_bc_pts[0] == num_node[0] ) SYS_T::print_fatal( "Error: the number of outline points is %d and the number of total points on the surface is %d. This is likely due to an improper wall mesh. \n", num_out_bc_pts[0], num_node[0] );

  inf_active_area[0] = 0.0;
  face_area[0] = 0.0;

  intNA[0].resize( num_node[0] );

  // zero the container
  for(int ii=0; ii<num_node[0]; ++ii) intNA[0][ii] = 0.0;

  if( elemtype == 501 )
  {
    double eptx[3]; double epty[3]; double eptz[3];
    int node_idx[3]; double R[3];

    const int nqp_tri = 3;                       // num qua points
    QuadPts_Gauss_Triangle quad( nqp_tri );      // quadrature rule
    FEAElement_Triangle3_3D_der0 ele( nqp_tri ); // element

    for(int ee=0; ee<num_cell[0]; ++ee)
    {
      for(int ii=0; ii<3; ++ii)
      {
        node_idx[ii] = tri_ien[0][3*ee+ii];
        eptx[ii] = pt_xyz[0][ 3*node_idx[ii]+0 ];
        epty[ii] = pt_xyz[0][ 3*node_idx[ii]+1 ];
        eptz[ii] = pt_xyz[0][ 3*node_idx[ii]+2 ];
      }

      ele.buildBasis(&quad, eptx, epty, eptz);

      for(int qua=0; qua<3; ++qua)
      {
        ele.get_R( qua, R );

        const double gwts = ele.get_detJac(qua) * quad.get_qw(qua);

        for(int ii=0; ii<3; ++ii)
        {
          inf_active_area[0] += gwts * R[ii] * temp_sol[ tri_ien[0][3*ee+ii] ];
          face_area[0] += gwts * R[ii];

          intNA[0][node_idx[ii]] += gwts * R[ii];
        }
      } // end qua-loop
    } // end ee-loop
  }
  else if( elemtype == 502 )
  {
    double eptx[6]; double epty[6]; double eptz[6];
    int node_idx[6]; double R[6];

    const int nqp_tri = 6;                       // num qua points
    QuadPts_Gauss_Triangle quad( nqp_tri );      // quadrature rule
    FEAElement_Triangle6_3D_der0 ele( nqp_tri ); // element

    for(int ee=0; ee<num_cell[0]; ++ee)
    {
      for(int ii=0; ii<6; ++ii)
      {
        node_idx[ii] = tri_ien[0][6*ee+ii];
        eptx[ii] = pt_xyz[0][ 3*node_idx[ii]+0 ];
        epty[ii] = pt_xyz[0][ 3*node_idx[ii]+1 ];
        eptz[ii] = pt_xyz[0][ 3*node_idx[ii]+2 ];
      }

      ele.buildBasis(&quad, eptx, epty, eptz);

      for(int qua=0; qua<6; ++qua)
      {
        ele.get_R( qua, R );

        const double gwts = ele.get_detJac(qua) * quad.get_qw(qua);

        for(int ii=0; ii<6; ++ii)
        {
          inf_active_area[0] += gwts * R[ii] * temp_sol[ tri_ien[0][6*ee+ii] ];
          face_area[0] += gwts * R[ii];

          intNA[0][node_idx[ii]] += gwts * R[ii];
        }
      } // end qua-loop
    } // end ee-loop
  }
  else SYS_T::print_fatal("Error: unknown element type.\n");

  delete [] temp_sol;

  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_inflow specified by "<<inffile<<", with nodes on \n";
  std::cout<<"     "<<wallfile<<" excluded, is generated. \n";
  std::cout<<"     num_node: "<<num_node[0]<<", num_cell: "<<num_cell[0]<<'\n';
  std::cout<<"     centroid: "<<centroid[0](0)<<'\t'<<centroid[0](1)<<'\t'<<centroid[0](2)<<'\n';
  std::cout<<"     number of outline points is "<<num_out_bc_pts[0]<<'\n';
  std::cout<<"     outward normal is ["<<outnormal[0](0)<<'\t'<<outnormal[0](1)<<'\t'<<outnormal[0](2)<<"]. \n";
  std::cout<<"     area is "<<face_area[0]<<", and active area is "<<inf_active_area[0]<<'\n';
}


NodalBC_3D_inflow::NodalBC_3D_inflow( const std::vector<std::string> &inffileList,
    const std::string &wallfile,
    const int &nFunc,
    const std::vector<Vector_3> &in_outnormal,
    const int &elemtype )
: num_nbc( static_cast<int>( inffileList.size() ) )
{
  // 1. No periodic nodes
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes.clear();

  // 2. Analyze the file type and read in the data
  num_node.resize( num_nbc );
  num_cell.resize( num_nbc );
  nLocBas.resize(  num_nbc );

  tri_ien.resize(num_nbc );
  pt_xyz.resize( num_nbc );

  global_node.resize(  num_nbc );
  global_cell.resize(  num_nbc );

  dir_nodes.resize(     num_nbc );
  num_dir_nodes.resize( num_nbc );

  centroid.resize(  num_nbc );
  outnormal.resize( num_nbc );

  outline_pts.resize(    num_nbc );
  num_out_bc_pts.resize( num_nbc );

  inf_active_area.resize( num_nbc );
  face_area.resize(       num_nbc );

  intNA.resize( num_nbc ); 

  SYS_T::file_check(wallfile);

  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  for( int ii=0; ii<num_nbc; ++ii )
  {
    SYS_T::file_check( inffileList[ii] );

    if( elemtype == 501 )
    {
      nLocBas[ii] = 3;

      TET_T::read_vtp_grid( inffileList[ii], num_node[ii], num_cell[ii],
        pt_xyz[ii], tri_ien[ii], global_node[ii], global_cell[ii] );

      TET_T::read_vtp_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );

    }
    else if( elemtype == 502 )
    {
      nLocBas[ii] = 6;

      TET_T::read_vtu_grid( inffileList[ii], num_node[ii], num_cell[ii],
        pt_xyz[ii], tri_ien[ii], global_node[ii], global_cell[ii] );

      TET_T::read_vtu_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );

    }

    else SYS_T::print_fatal("Error: unknown element type.\n");

    // Generate the dir-node list. Nodes belonging to the wall are excluded.
    dir_nodes[ii].clear();
    
    for(unsigned int jj=0; jj<global_node[ii].size(); ++jj)
    {
      SYS_T::print_fatal_if( global_node[ii][jj]<0, "Error: negative nodal index! \n");

      if( !VEC_T::is_invec( wall_gnode, global_node[ii][jj]) )
        dir_nodes[ii].push_back( global_node[ii][jj] );
    }

    num_dir_nodes[ii] = dir_nodes[ii].size(); 

    // Calculate the centroid of the surface
    centroid[ii].gen_zero();
    for(int jj=0; jj<num_node[ii]; ++jj)
    {
      centroid[ii](0) += pt_xyz[ii][3*jj+0];
      centroid[ii](1) += pt_xyz[ii][3*jj+1];
      centroid[ii](2) += pt_xyz[ii][3*jj+2];
    }
    centroid[ii].scale( 1.0 / (double) num_node[ii] );

    // assign outward normal vector from the input
    outnormal[ii] = in_outnormal[ii];

    // Collect nodes that belong to the wall, and set up a vector that
    // is 1 on the interior nodes and 0 on the wall bc nodes.
    outline_pts[ii].clear();

    num_out_bc_pts[ii] = 0;

    double * temp_sol = new double [num_node[ii]];  

    for(int jj=0; jj<num_node[ii]; ++jj)
    {
      // If the node is not on the wall, it is an interior node, so set
      // the element to 1.
      if( !VEC_T::is_invec(wall_gnode, global_node[ii][jj]) ) 
        temp_sol[jj] = 1.0;

      // otherwise, the node is on the wall surface, so set element to 0.
      else 
      {
        temp_sol[jj] = 0.0;

        // Also store the point's coordinates in outline points
        num_out_bc_pts[ii] += 1;
        outline_pts[ii].push_back( pt_xyz[ii][3*jj+0] );
        outline_pts[ii].push_back( pt_xyz[ii][3*jj+1] );
        outline_pts[ii].push_back( pt_xyz[ii][3*jj+2] );
      }
    }

    // If the number of surface nodes matches num_out_bc_pts, the wall
    // mesh contains the inlet surface. This is a common error when
    // adopting sv files, where the user uses the combined exterior surface
    // as the wall mesh. We will throw an error message if detected.
    if( num_out_bc_pts[ii] == num_node[ii] ) SYS_T::print_fatal( "Error: the number of outline points is %d and the number of total points on the surface is %d. This is likely due to an improper wall mesh. \n", num_out_bc_pts[ii], num_node[ii] );

    inf_active_area[ii] = 0.0;
    face_area[ii] = 0.0;

    intNA[ii].resize( num_node[ii] );

    // zero the container
    for(int jj=0; jj<num_node[ii]; ++jj) intNA[ii][jj] = 0.0;

    if( elemtype == 501 )
    {
      double eptx[3]; double epty[3]; double eptz[3];
      int node_idx[3]; double R[3];

      const int nqp_tri = 3;                       // num qua points
      QuadPts_Gauss_Triangle quad( nqp_tri );      // quadrature rule
      FEAElement_Triangle3_3D_der0 ele( nqp_tri ); // element

      for(int ee=0; ee<num_cell[ii]; ++ee)
      {
        for(int jj=0; jj<3; ++jj)
        {
          node_idx[jj] = tri_ien[ii][3*ee+jj];
          eptx[jj] = pt_xyz[ii][ 3*node_idx[jj]+0 ];
          epty[jj] = pt_xyz[ii][ 3*node_idx[jj]+1 ];
          eptz[jj] = pt_xyz[ii][ 3*node_idx[jj]+2 ];
        }

        ele.buildBasis(&quad, eptx, epty, eptz);

        for(int qua=0; qua<nqp_tri; ++qua)
        {
          ele.get_R( qua, R );

          const double gwts = ele.get_detJac(qua) * quad.get_qw(qua);

          for(int jj=0; jj<3; ++jj)
          {
            inf_active_area[ii] += gwts * R[jj] * temp_sol[ tri_ien[ii][3*ee+jj] ];
            face_area[ii] += gwts * R[jj];

            intNA[ii][node_idx[jj]] += gwts * R[jj];
          }
        } // end qua-loop
      } // end ee-loop
    }
    else if( elemtype == 502 )
    {
      double eptx[6]; double epty[6]; double eptz[6]; 
      int node_idx[6]; double R[6];

      const int nqp_tri = 6;                       // num qua points
      QuadPts_Gauss_Triangle quad( nqp_tri );      // quadrature rule
      FEAElement_Triangle6_3D_der0 ele( nqp_tri ); // element

      for(int ee=0; ee<num_cell[ii]; ++ee)
      {
        for(int jj=0; jj<6; ++jj)
        {
          node_idx[jj] = tri_ien[ii][6*ee+jj];
          eptx[jj] = pt_xyz[ii][ 3*node_idx[jj]+0 ];
          epty[jj] = pt_xyz[ii][ 3*node_idx[jj]+1 ];
          eptz[jj] = pt_xyz[ii][ 3*node_idx[jj]+2 ];
        }

        ele.buildBasis(&quad, eptx, epty, eptz);

        for(int qua=0; qua<nqp_tri; ++qua)
        {
          ele.get_R( qua, R );

          const double gwts = ele.get_detJac(qua) * quad.get_qw(qua);

          for(int jj=0; jj<6; ++jj)
          {
            inf_active_area[ii] += gwts * R[jj] * temp_sol[ tri_ien[ii][6*ee+jj] ];
            face_area[ii] += gwts * R[jj];

            intNA[ii][node_idx[jj]] += gwts * R[jj];
          }
        } // end qua-loop
      } // end ee-loop
    }
    else SYS_T::print_fatal("Error: unknown element type.\n");

    delete [] temp_sol;

  } // end ii-loop

  // Generate ID array
  Create_ID( nFunc );

  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_inflow specified by\n";
  for(int ii=0; ii<num_nbc; ++ii)
  {
    std::cout<<"     nbc_id = "<<ii<<": "<<inffileList[ii]<<" with nodes on ";
    std::cout<<"     "<<wallfile<<" excluded.\n";
    std::cout<<"          num_node: "<<num_node[ii]<<", num_cell: "<<num_cell[ii]<<'\n';
    std::cout<<"          centroid: "<<centroid[ii](0)<<'\t'<<centroid[ii](1)<<'\t'<<centroid[ii](2)<<'\n';
    std::cout<<"          number of outline points is "<<num_out_bc_pts[ii]<<'\n';
    std::cout<<"          outward normal is ["<<outnormal[ii](0)<<'\t'<<outnormal[ii](1)<<'\t'<<outnormal[ii](2)<<"]. \n";
    std::cout<<"          area is "<<face_area[ii]<<", and active area is "<<inf_active_area[ii]<<'\n';
  }
}

// EOF
