#include "NodalBC_3D_inflow.hpp"

NodalBC_3D_inflow::NodalBC_3D_inflow(const int &nFunc) : num_nbc( 0 )
{
  dir_nodes.clear();
  num_dir_nodes = 0;

  dir_nodes_on_inlet.clear();
  num_dir_nodes_on_inlet.clear();

  Create_ID( nFunc );

  inf_active_area.clear(); face_area.clear();

  centroid.clear(); outnormal.clear();
  num_out_bc_pts.clear(); outline_pts.clear();
  intNA.clear();

  num_node.clear(); num_cell.clear(); nLocBas.clear();
  tri_ien.clear(); pt_xyz.clear(); global_node.clear(); global_cell.clear();

  std::cout<<"===> NodalBC_3D_inflow::empty is generated. \n";
}

NodalBC_3D_inflow::NodalBC_3D_inflow( const std::string &inffile,
    const std::string &wallfile, const int &nFunc,
    const Vector_3 &in_outnormal,
    const int &elemtype ) : num_nbc( 1 )
{
  SYS_T::file_check(inffile);
  SYS_T::file_check(wallfile);

  const int nbc_id = 0;

  // 1. Clear the container for Dirichlet nodes
  dir_nodes.clear();
  num_dir_nodes = 0;

  dir_nodes_on_inlet.resize( num_nbc );
  for(int ii=0; ii<num_nbc; ++ii) dir_nodes_on_inlet[ii].clear();

  num_dir_nodes_on_inlet.resize( num_nbc );

  // 2. Analyze the file type and read in the data
  num_node.resize( num_nbc );
  num_cell.resize( num_nbc );
  nLocBas.resize(  num_nbc );

  tri_ien.resize(num_nbc );
  pt_xyz.resize( num_nbc );

  global_node.resize(  num_nbc );
  global_cell.resize(  num_nbc );

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
    nLocBas[nbc_id] = 3;

    TET_T::read_vtp_grid( inffile, num_node[nbc_id], num_cell[nbc_id], pt_xyz[nbc_id], tri_ien[nbc_id], global_node[nbc_id], global_cell[nbc_id] );

    TET_T::read_vtp_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );
  }
  else if( elemtype == 502 )
  {
    nLocBas[nbc_id] = 6;

    TET_T::read_vtu_grid( inffile, num_node[nbc_id], num_cell[nbc_id], pt_xyz[nbc_id], tri_ien[nbc_id], global_node[nbc_id], global_cell[nbc_id] );

    TET_T::read_vtu_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );
  }
  else SYS_T::print_fatal("Error: unknown element type.\n");

  // Generate the dir-node list. Nodes belonging to the wall are excluded.
  for(unsigned int ii=0; ii<global_node[nbc_id].size(); ++ii)
  {
    SYS_T::print_fatal_if( global_node[nbc_id][ii]<0, "Error: negative nodal index! \n");

    if( !VEC_T::is_invec( wall_gnode, global_node[nbc_id][ii]) )
      dir_nodes.push_back( global_node[nbc_id][ii] );
  }

  num_dir_nodes = dir_nodes.size(); // Now assign the value of the number of Dirichlet nodes 

  // Generate ID array
  Create_ID( nFunc );

  // Calculate the centroid of the surface
  centroid[nbc_id].gen_zero();
  for(int ii=0; ii<num_node[nbc_id]; ++ii)
  {
    centroid[nbc_id](0) += pt_xyz[nbc_id][3*ii+0];
    centroid[nbc_id](1) += pt_xyz[nbc_id][3*ii+1];
    centroid[nbc_id](2) += pt_xyz[nbc_id][3*ii+2];
  }
  centroid[nbc_id].scale( 1.0 / (double) num_node[nbc_id] );

  // assign outward normal vector from the input
  outnormal[nbc_id] = in_outnormal;

  // Collect nodes that belong to the wall, and set up a vector that
  // is 1 on the interior nodes and 0 on the wall bc nodes.
  outline_pts[nbc_id].clear();

  num_out_bc_pts[nbc_id] = 0;

  double * temp_sol = new double [num_node[nbc_id]];  

  for(int ii=0; ii<num_node[nbc_id]; ++ii)
  {
    // If the node is not on the wall, it is an interior node, so set
    // the element to 1.
    if( !VEC_T::is_invec(wall_gnode, global_node[nbc_id][ii]) ) 
      temp_sol[ii] = 1.0;

    // otherwise, the node is on the wall surface, so set element to 0.
    else 
    {
      temp_sol[ii] = 0.0;

      // Also store the point's coordinates in outline points
      num_out_bc_pts[nbc_id] += 1;
      outline_pts[nbc_id].push_back( pt_xyz[nbc_id][3*ii+0] );
      outline_pts[nbc_id].push_back( pt_xyz[nbc_id][3*ii+1] );
      outline_pts[nbc_id].push_back( pt_xyz[nbc_id][3*ii+2] );
    }
  }

  // If the number of surface nodes matches num_out_bc_pts, the wall
  // mesh contains the inlet surface. This is a common error when
  // adopting sv files, where the user uses the combined exterior surface
  // as the wall mesh. We will throw an error message if detected.
  if( num_out_bc_pts[nbc_id] == num_node[nbc_id] ) SYS_T::print_fatal( "Error: the number of outline points is %d and the number of total points on the surface is %d. This is likely due to an improper wall mesh. \n", num_out_bc_pts[nbc_id], num_node[nbc_id] );

  inf_active_area[nbc_id] = 0.0;
  face_area[nbc_id] = 0.0;

  intNA[nbc_id].resize( num_node[nbc_id] );

  // zero the container
  for(int ii=0; ii<num_node[nbc_id]; ++ii) intNA[nbc_id][ii] = 0.0;

  if( elemtype == 501 )
  {
    double eptx[3]; double epty[3]; double eptz[3];
    int node_idx[3]; double R[3];

    const int nqp_tri = 3;                       // num qua points
    QuadPts_Gauss_Triangle quad( nqp_tri );      // quadrature rule
    FEAElement_Triangle3_3D_der0 ele( nqp_tri ); // element

    for(int ee=0; ee<num_cell[nbc_id]; ++ee)
    {
      for(int ii=0; ii<3; ++ii)
      {
        node_idx[ii] = tri_ien[nbc_id][3*ee+ii];
        eptx[ii] = pt_xyz[nbc_id][ 3*node_idx[ii]+0 ];
        epty[ii] = pt_xyz[nbc_id][ 3*node_idx[ii]+1 ];
        eptz[ii] = pt_xyz[nbc_id][ 3*node_idx[ii]+2 ];
      }

      ele.buildBasis(&quad, eptx, epty, eptz);

      for(int qua=0; qua<3; ++qua)
      {
        ele.get_R( qua, R );

        const double gwts = ele.get_detJac(qua) * quad.get_qw(qua);

        for(int ii=0; ii<3; ++ii)
        {
          inf_active_area[nbc_id] += gwts * R[ii] * temp_sol[ tri_ien[nbc_id][3*ee+ii] ];
          face_area[nbc_id] += gwts * R[ii];

          intNA[nbc_id][node_idx[ii]] += gwts * R[ii];
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

    for(int ee=0; ee<num_cell[nbc_id]; ++ee)
    {
      for(int ii=0; ii<6; ++ii)
      {
        node_idx[ii] = tri_ien[nbc_id][6*ee+ii];
        eptx[ii] = pt_xyz[nbc_id][ 3*node_idx[ii]+0 ];
        epty[ii] = pt_xyz[nbc_id][ 3*node_idx[ii]+1 ];
        eptz[ii] = pt_xyz[nbc_id][ 3*node_idx[ii]+2 ];
      }

      ele.buildBasis(&quad, eptx, epty, eptz);

      for(int qua=0; qua<6; ++qua)
      {
        ele.get_R( qua, R );

        const double gwts = ele.get_detJac(qua) * quad.get_qw(qua);

        for(int ii=0; ii<6; ++ii)
        {
          inf_active_area[nbc_id] += gwts * R[ii] * temp_sol[ tri_ien[nbc_id][6*ee+ii] ];
          face_area[nbc_id] += gwts * R[ii];

          intNA[nbc_id][node_idx[ii]] += gwts * R[ii];
        }
      } // end qua-loop
    } // end ee-loop
  }
  else SYS_T::print_fatal("Error: unknown element type.\n");

  delete [] temp_sol; temp_sol = nullptr;

  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_inflow specified by "<<inffile<<", with nodes on \n";
  std::cout<<"     "<<wallfile<<" excluded, is generated. \n";
  std::cout<<"     num_node: "<<num_node[nbc_id]<<", num_cell: "<<num_cell[nbc_id]<<'\n';
  std::cout<<"     centroid: "<<centroid[nbc_id](0)<<'\t'<<centroid[nbc_id](1)<<'\t'<<centroid[nbc_id](2)<<'\n';
  std::cout<<"     number of outline points is "<<num_out_bc_pts[nbc_id]<<'\n';
  std::cout<<"     outward normal is ["<<outnormal[nbc_id](0)<<'\t'<<outnormal[nbc_id](1)<<'\t'<<outnormal[nbc_id](2)<<"]. \n";
  std::cout<<"     area is "<<face_area[nbc_id]<<", and active area is "<<inf_active_area[nbc_id]<<'\n';
}


NodalBC_3D_inflow::NodalBC_3D_inflow( const std::vector<std::string> &inffileList,
    const std::string &wallfile,
    const int &nFunc,
    const std::vector<Vector_3> &in_outnormal,
    const int &elemtype ) : num_nbc( static_cast<int>( inffileList.size() ) )
{
  init( inffileList, wallfile, nFunc, in_outnormal, elemtype );
}

void NodalBC_3D_inflow::init( const std::vector<std::string> &inffileList,
    const std::string &wallfile,
    const int &nFunc,
    const std::vector<Vector_3> &in_outnormal,
    const int &elemtype )
{ 
  SYS_T::file_check(wallfile);

  // 1. Clear the container for Dirichlet nodes
  dir_nodes.clear();
  num_dir_nodes = 0;

  dir_nodes_on_inlet.resize( num_nbc );
  for(int ii=0; ii<num_nbc; ++ii) dir_nodes_on_inlet[ii].clear();

  num_dir_nodes_on_inlet.resize( num_nbc );

  // 2. Analyze the file type and read in the data
  num_node.resize( num_nbc );
  num_cell.resize( num_nbc );
  nLocBas.resize(  num_nbc );

  tri_ien.resize(num_nbc );
  pt_xyz.resize( num_nbc );

  global_node.resize(  num_nbc );
  global_cell.resize(  num_nbc );

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
    for(unsigned int jj=0; jj<global_node[ii].size(); ++jj)
    {
      SYS_T::print_fatal_if( global_node[ii][jj]<0, "Error: negative nodal index! \n");

      if( !VEC_T::is_invec( wall_gnode, global_node[ii][jj]) )
      {
        dir_nodes.push_back( global_node[ii][jj] );
        dir_nodes_on_inlet[ii].push_back( global_node[ii][jj] );
      }
    }

    num_dir_nodes_on_inlet[ii] = dir_nodes_on_inlet[ii].size();

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

    delete [] temp_sol; temp_sol = nullptr;

  } // end ii-loop

  num_dir_nodes = dir_nodes.size();

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
