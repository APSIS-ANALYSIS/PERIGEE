#include "NodalBC_3D_inflow.hpp"

NodalBC_3D_inflow::NodalBC_3D_inflow(const int &nFunc)
{
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
  dir_nodes.clear();
  num_dir_nodes = 0;

  Create_ID( nFunc );
  
  inf_active_area = 0;

  centroid[0] = 0.0; centroid[1] = 0.0; centroid[2] = 0.0;

  outnormal.gen_zero();
  
  num_out_bc_pts = 0;
  VEC_T::clean(outline_pts);
  VEC_T::clean(intNA);

  std::cout<<"===> NodalBC_3D_inflow::empty is generated. \n";
}


NodalBC_3D_inflow::NodalBC_3D_inflow( const std::string &inffile,
    const std::string &wallfile, const int &nFunc,
    const Vector_3 &in_outnormal,
    const int &elemtype )
{
  SYS_T::file_check(inffile);
  SYS_T::file_check(wallfile);

  // 1. There is no periodic nodes
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  // 2. Analyze the file type and read in the data
  if( elemtype == 501 ) nLocBas = 3;
  else if( elemtype == 502 ) nLocBas = 6;
  else SYS_T::print_fatal("Error: unknown element type.\n");

  // Read the files
  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  std::string fend; fend.assign( inffile.end()-4 , inffile.end() );

  if( fend.compare(".vtp") == 0 )
  { 
    TET_T::read_vtp_grid( inffile, num_node, num_cell, pt_xyz, tri_ien, global_node, global_cell );

    TET_T::read_vtp_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );
  }
  else if( fend.compare(".vtu") == 0 )
  {
    TET_T::read_vtu_grid( inffile, num_node, num_cell, pt_xyz, tri_ien, global_node, global_cell );

    TET_T::read_vtu_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );
  }
  else
    SYS_T::print_fatal("Error: Nodal_3D_inflow unknown file type.\n");

  // Generate the dir-node list. The nodes belonging to the wall are excluded.
  dir_nodes.clear();
  for(unsigned int ii=0; ii<global_node.size(); ++ii)
  {
    SYS_T::print_fatal_if( global_node[ii]<0, "Error: there are negative nodal index! \n");

    if( !VEC_T::is_invec( wall_gnode, global_node[ii]) )
      dir_nodes.push_back( global_node[ii] );
  }

  num_dir_nodes = dir_nodes.size(); 

  // Generate ID array
  Create_ID( nFunc );

  // Calculate the centroid of the surface
  centroid[0] = 0.0; centroid[1] = 0.0; centroid[2] = 0.0;
  for(int ii=0; ii<num_node; ++ii)
  {
    centroid[0] += pt_xyz[3*ii+0];
    centroid[1] += pt_xyz[3*ii+1];
    centroid[2] += pt_xyz[3*ii+2];
  }
  centroid[0] = centroid[0] / (double) num_node;
  centroid[1] = centroid[1] / (double) num_node;
  centroid[2] = centroid[2] / (double) num_node;

  // Collect the nodes that belong to the wall, and setup a vector that
  // is 1 on the interior nodes and 0 on the wall bc nodes.
  outline_pts.clear();
  num_out_bc_pts = 0;
  double * temp_sol = new double [num_node];  
  for(int ii=0; ii<num_node; ++ii)
  {
    // If the node is not in the wall, it is an interior node and set
    // the vector to be 1.
    if( !VEC_T::is_invec(wall_gnode, global_node[ii]) ) 
      temp_sol[ii] = 1.0;
    else 
    {
      // otherwise, the node is on the wall surface, set the vector to
      // be 0.
      temp_sol[ii] = 0.0;
      // Also store the point's coordinates into outline points
      num_out_bc_pts += 1;
      outline_pts.push_back( pt_xyz[3*ii+0] );
      outline_pts.push_back( pt_xyz[3*ii+1] );
      outline_pts.push_back( pt_xyz[3*ii+2] );
    }
  }

  // If the number of the surface nodes matches with the num_out_bc_pts, this
  // means the wall mesh contains the inlet surface. This is a common error
  // happens from adopting sv files, where the user uses the combined outer
  // surface as the wall mesh. We will throw an error message if this is
  // detected.
  if( num_out_bc_pts == num_node ) SYS_T::print_fatal( "Error: the number of outline points is %d and the number of total points on the surface is %d. This is likely due to an improper wall mesh. \n", num_out_bc_pts, num_node );

  inf_active_area = 0.0;
  face_area = 0.0;

  if( elemtype == 501 )
  {
    double eptx[3]; double epty[3]; double eptz[3]; double R[3];

    QuadPts_Gauss_Triangle quad(3); // quadrature rule
    FEAElement_Triangle3_3D_der0 ele(3); // element

    for(int ee=0; ee<num_cell; ++ee)
    {
      for(int ii=0; ii<3; ++ii)
      {
        const int nodidx = tri_ien[3*ee+ii];
        eptx[ii] = pt_xyz[ 3*nodidx ];
        epty[ii] = pt_xyz[ 3*nodidx+1 ];
        eptz[ii] = pt_xyz[ 3*nodidx+2 ];
      }

      ele.buildBasis(&quad, eptx, epty, eptz);

      for(int qua=0; qua<3; ++qua)
      {
        ele.get_R( qua, R );

        for(int ii=0; ii<3; ++ii)
        {
          inf_active_area += ele.get_detJac(qua) * quad.get_qw(qua) 
            * R[ii] * temp_sol[ tri_ien[3*ee+ii] ];

          face_area += ele.get_detJac(qua) * quad.get_qw(qua) * R[ii];
        }
      }
    }
  }
  else if( elemtype == 502 )
  {
    double eptx[6]; double epty[6]; double eptz[6]; double R[6];

    QuadPts_Gauss_Triangle quad(6); // quadrature rule
    FEAElement_Triangle6_3D_der0 ele(6); // element

    for(int ee=0; ee<num_cell; ++ee)
    {
      for(int ii=0; ii<6; ++ii)
      {
        const int nodidx = tri_ien[6*ee+ii];
        eptx[ii] = pt_xyz[ 3*nodidx ];
        epty[ii] = pt_xyz[ 3*nodidx+1 ];
        eptz[ii] = pt_xyz[ 3*nodidx+2 ];
      }

      ele.buildBasis(&quad, eptx, epty, eptz);

      for(int qua=0; qua<6; ++qua)
      {
        ele.get_R( qua, R );

        for(int ii=0; ii<6; ++ii)
        {
          inf_active_area += ele.get_detJac(qua) * quad.get_qw(qua) 
            * R[ii] * temp_sol[ tri_ien[6*ee+ii] ];

          face_area += ele.get_detJac(qua) * quad.get_qw(qua) * R[ii];
        }
      }
    }
  }
  else SYS_T::print_fatal("Error: unknown element type.\n");

  delete [] temp_sol;

  // assign outward normal vector from the input
  outnormal = in_outnormal;

  // Perform surface integral
  intNA.resize( num_node );

  // zero the container
  for(int ii=0; ii<num_node; ++ii) intNA[ii] = 0.0;

  if( elemtype == 501 )
  {
    const int nqp_tri = 3; // number of quadrature points
    IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );
    FEAElement * elems = new FEAElement_Triangle3_3D_der0( nqp_tri );

    double ectrl_x[3]; double ectrl_y[3]; double ectrl_z[3];
    int node_idx[3]; double R[3];

    // Calculate the surface integral of basis functions
    for( int ee = 0; ee<num_cell; ++ee )
    {
      for(int ii=0; ii<3; ++ii)
      {
        node_idx[ii] = tri_ien[3*ee+ii];
        ectrl_x[ii] = pt_xyz[3*node_idx[ii] + 0];
        ectrl_y[ii] = pt_xyz[3*node_idx[ii] + 1];
        ectrl_z[ii] = pt_xyz[3*node_idx[ii] + 2];
      }

      elems -> buildBasis(quads, ectrl_x, ectrl_y, ectrl_z);

      for(int qua=0; qua<nqp_tri; ++qua)
      {
        elems -> get_R(qua, R);

        const double gwts = elems->get_detJac( qua ) * quads -> get_qw( qua );

        for(int ii=0; ii<3; ++ii) intNA[node_idx[ii]] += gwts * R[ii];

      } // loop over quadrature points
    } // loop over linear triangle elements

    delete quads; delete elems;
  }
  else if(elemtype == 502 )
  {
    const int nqp_tri = 6; // number of quadrature points
    IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );
    FEAElement * elems = new FEAElement_Triangle6_3D_der0( nqp_tri );

    double ectrl_x[6]; double ectrl_y[6]; double ectrl_z[6];
    int node_idx[6]; double R[6];

    // Calculate the surface integral for quadratic triangle element
    for( int ee=0; ee<num_cell; ++ee )
    {
      for(int ii=0; ii<6; ++ii)
      {
        node_idx[ii] = tri_ien[6*ee+ii];
        ectrl_x[ii] = pt_xyz[3*node_idx[ii] + 0];
        ectrl_y[ii] = pt_xyz[3*node_idx[ii] + 1];
        ectrl_z[ii] = pt_xyz[3*node_idx[ii] + 2];
      }

      elems -> buildBasis(quads, ectrl_x, ectrl_y, ectrl_z);

      for(int qua=0; qua<nqp_tri; ++qua)
      {
        elems -> get_R(qua, R);

        const double gwts = elems->get_detJac( qua ) * quads -> get_qw( qua );

        for(int ii=0; ii<6; ++ii) intNA[node_idx[ii]] += gwts * R[ii];

      } // loop over quadrature points
    } // loop over quadratic triangle elements

    delete quads; delete elems;
  }
  else SYS_T::print_fatal("Error: unknown element type.\n");

  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_inflow specified by "<<inffile<<", with nodes on \n";
  std::cout<<"     "<<wallfile<<" excluded, is generated. \n";
  std::cout<<"     num_node: "<<num_node<<", num_cell: "<<num_cell<<'\n';
  std::cout<<"     centroid: "<<centroid[0]<<'\t'<<centroid[1]<<'\t'<<centroid[2]<<'\n';
  std::cout<<"     number of outline points is "<<num_out_bc_pts<<'\n';
  std::cout<<"     outward normal is ["<<outnormal(0)<<'\t'<<outnormal(1)<<'\t'<<outnormal(2)<<"]. \n";
  std::cout<<"     area is "<<face_area<<", and active area is "<<inf_active_area<<'\n';
}

// EOF
