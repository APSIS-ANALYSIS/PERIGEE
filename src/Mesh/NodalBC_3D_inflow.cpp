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

  outnormal.resize(3);
  outnormal[0] = 0.0; outnormal[1] = 0.0; outnormal[2] = 0.0;

  num_out_bc_pts = 0;
  VEC_T::clean(outline_pts);
  VEC_T::clean(intNA);

  std::cout<<"===> NodalBC_3D_inflow::empty is generated. \n";
}


NodalBC_3D_inflow::NodalBC_3D_inflow( const std::string &inffile,
    const std::string &wallfile, const int &nFunc,
    const std::vector<double> &in_outnormal,
    const int &elemtype )
{
  SYS_T::file_check(inffile);
  SYS_T::file_check(wallfile);

  // 1. There is no periodic nodes
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  // 2. Analyze the file type and read in the data
  // Read the files
  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  int wall_numpts, wall_numcels;
  std::vector<double> wall_pts;
  std::vector<int> wall_ien, wall_gnode, wall_gelem;

  std::string fend; fend.assign( inffile.end()-4 , inffile.end() );
 
  if( fend.compare(".vtp") == 0 )
  { 
    TET_T::read_vtp_grid( inffile, numpts, numcels, pts, ien, gnode, gelem );

    TET_T::read_vtp_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );
  }
  else if( fend.compare(".vtu") == 0 )
  {
    TET_T::read_vtu_grid( inffile, numpts, numcels, pts, ien, gnode, gelem );

    TET_T::read_vtu_grid( wallfile, wall_numpts, wall_numcels, wall_pts, 
        wall_ien, wall_gnode, wall_gelem );
  }
  else
    SYS_T::print_fatal("Error: get_out_normal unknown file type.\n");

  // Generate the dir-node list. The nodes belonging to the wall are excluded.
  dir_nodes.clear();
  for(unsigned int ii=0; ii<gnode.size(); ++ii)
  {
    SYS_T::print_fatal_if( gnode[ii]<0,
        "Error: there are negative nodal index! \n");

    if( !VEC_T::is_invec( wall_gnode, gnode[ii]) )
      dir_nodes.push_back( gnode[ii] );
  }

  num_dir_nodes = dir_nodes.size(); 

  // Generate ID array
  Create_ID( nFunc );

  // Calculate the centroid of the surface
  centroid[0] = 0.0; centroid[1] = 0.0; centroid[2] = 0.0;
  for(int ii=0; ii<numpts; ++ii)
  {
    centroid[0] += pts[3*ii+0];
    centroid[1] += pts[3*ii+1];
    centroid[2] += pts[3*ii+2];
  }
  centroid[0] = centroid[0] / (double) numpts;
  centroid[1] = centroid[1] / (double) numpts;
  centroid[2] = centroid[2] / (double) numpts;

  // Collect the nodes that belong to the wall, and setup a vector that
  // is 1 on the interior nodes and 0 on the wall bc nodes.
  outline_pts.clear();
  num_out_bc_pts = 0;
  double * temp_sol = new double [numpts];  
  for(int ii=0; ii<numpts; ++ii)
  {
    // If the node is not in the wall, it is an interior node and set
    // the vector to be 1.
    if( !VEC_T::is_invec(wall_gnode, gnode[ii]) ) 
      temp_sol[ii] = 1.0;
    else 
    {
      // otherwise, the node is on the wall surface, set the vector to
      // be 0.
      temp_sol[ii] = 0.0;
      // Also store the point's coordinates into outline points
      num_out_bc_pts += 1;
      outline_pts.push_back( pts[3*ii+0] );
      outline_pts.push_back( pts[3*ii+1] );
      outline_pts.push_back( pts[3*ii+2] );
    }
  }

  inf_active_area = 0.0;
  face_area = 0.0;

  double eptx[3]; double epty[3]; double eptz[3]; double R[3];
  double nx, ny, nz, ele_surface_area;
  
  QuadPts_Gauss_Triangle quad(3); // quadrature rule
  FEAElement_Triangle3_3D_der0 ele(3); // element

  for(int ee=0; ee<numcels; ++ee)
  {
    for(int ii=0; ii<3; ++ii)
    {
      const int nodidx = ien[3*ee+ii];
      eptx[ii] = pts[ 3*nodidx ];
      epty[ii] = pts[ 3*nodidx+1 ];
      eptz[ii] = pts[ 3*nodidx+2 ];
    }

    ele.buildBasis(&quad, eptx, epty, eptz);

    for(int qua=0; qua<3; ++qua)
    {
      ele.get_R( qua, R );
      ele.get_2d_normal_out(qua, nx, ny, nz, ele_surface_area);

      for(int ii=0; ii<3; ++ii)
      {
        inf_active_area += ele_surface_area * quad.get_qw(qua) 
          * R[ii] * temp_sol[ ien[3*ee+ii] ];

        face_area += ele_surface_area * quad.get_qw(qua) * R[ii];
      }
    }
  }
  delete [] temp_sol;

  // assign outward normal vector from the input
  outnormal = in_outnormal;

  // Perform surface integral
  intNA.resize( gnode.size() );
 
  IQuadPts * quads = nullptr;
  FEAElement * elems = nullptr;
   
  if( elemtype == 501 )
  {
    quads = new QuadPts_Gauss_Triangle( 3 );
    elems = new FEAElement_Triangle3_3D_der0( 3 );
  }
  else if(elemtype == 502 )
  {
    quads = new QuadPts_Gauss_Triangle( 6 );
    elems = new FEAElement_Triangle6_3D_der0( 6 );
  }
  else SYS_T::print_fatal("Error: unknown element type.\n");




  if(quads != nullptr) delete quads;

  if(elems != nullptr) delete elems;

  // Finish and print info on screen
  std::cout<<"===> NodalBC_3D_inflow specified by "<<inffile
    <<" and wall file "<<wallfile<<" is generated. \n";
  std::cout<<"     centroid: "<<centroid[0]<<'\t'<<centroid[1]<<'\t'<<centroid[2]<<'\n';
  std::cout<<"     number of outline points is "<<num_out_bc_pts<<'\n';
  std::cout<<"     outward normal is ["<<outnormal[0]<<'\t'<<outnormal[1]<<'\t'<<outnormal[2]<<"]. \n";
  std::cout<<"     area is "<<face_area<<", and active area is "<<inf_active_area<<'\n';
}

// EOF
