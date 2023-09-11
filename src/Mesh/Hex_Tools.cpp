#include "Hex_Tools.hpp"

void HEX_T::gen_hex_grid( vtkUnstructuredGrid * const &grid_w,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array )
{
  // Check the input data compatibility
  SYS_T::print_exit_if( int(pt.size()) != 3*numpts,
    "Error: HEX_T::gen_hex_grid point, vector size does not match the number of points. \n" );

  // detect the element type
  int nlocbas {-1};
  if ( int(ien_array.size()) == 8*numcels )
    nlocbas = 8;
  else if ( int(ien_array.size()) == 27*numcels ) 
    nlocbas = 27;
  else
    SYS_T::print_exit("Error: HEX_T::gen_hex_grid, ien array size does not match the number of cells. \n");
  // Now we donnot consider the serendipity elements i.e. 8-node quadrilateral and 20-node hexahedron.

  // Check the connectivity array
  std::vector<int> temp = ien_array;
  VEC_T::sort_unique_resize(temp);
  SYS_T::print_exit_if( int(temp.size()) != numpts,
  "Error: HEX_T::gen_hex_grid, numpts does not match the number of unique points in the ien array. Please re-organize the input. \n" );
  VEC_T::clean(temp);

  // 1. nodal points coordinates
  vtkPoints * ppt = vtkPoints::New(); 
  ppt->SetDataTypeToDouble();
  for(int ii=0; ii<numpts; ++ii)
  {
    const double coor[3] = { pt[3*ii], pt[3*ii+1], pt[3*ii+2] };
    ppt -> InsertPoint(ii, coor);
  }

  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  // 2. Cell
  // Assign an additional data type for element aspect ratio
  vtkDoubleArray * edge_aspect_ratio = vtkDoubleArray::New();
  edge_aspect_ratio -> SetName("Aspect_ratio");
  edge_aspect_ratio -> SetNumberOfComponents(1);

  if(nlocbas == 8) 
  {
    vtkCell * cl = vtkHexahedron::New();

    for(int ii=0; ii<numcels; ++ii)
    {
      for(int lnode=0; lnode<8; ++lnode)
        cl->GetPointIds()->SetId( lnode, ien_array[8*ii + lnode] );

      grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

      std::vector<double> cell_node; cell_node.clear();
      for(int lnode=0; lnode<8; ++lnode)
      {
        int node_offset = 3 * ien_array[8*ii + lnode];
        cell_node.push_back(pt[node_offset]);
        cell_node.push_back(pt[node_offset+1]);
        cell_node.push_back(pt[node_offset+2]);
      }
      edge_aspect_ratio -> InsertNextValue( HEX_T::get_aspect_ratio(cell_node) );
    }

    cl -> Delete();
  }
  else if(nlocbas == 27)
  {
    vtkCell * cl = vtkTriQuadraticHexahedron::New();

    for(int ii=0; ii<numcels; ++ii)
    { 
      std::vector<int> local_ien_msh(ien_array.begin()+27*ii, ien_array.begin() + 27*ii + 27);
      std::vector<int> local_ien_vtk = HEX_T::reset_node(local_ien_msh);
      for(int lnode=0; lnode<27; ++lnode)
        cl->GetPointIds()->SetId( lnode, local_ien_vtk[lnode] );

      grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

      std::vector<double> cell_node {};
      for(int lnode=0; lnode<8; ++lnode)
      {
        int node_offset = 3 * local_ien_vtk[lnode];
        cell_node.push_back(pt[node_offset]);
        cell_node.push_back(pt[node_offset+1]);
        cell_node.push_back(pt[node_offset+2]);
      }
      edge_aspect_ratio -> InsertNextValue( HEX_T::get_aspect_ratio(cell_node) );
    }

    cl -> Delete();
  }
  else
    SYS_T::print_exit("Error: HEX_T::gen_hex_grid, unknown local basis number.\n");

  // Add the asepct-ratio to grid_w
  grid_w -> GetCellData() -> AddArray( edge_aspect_ratio );
  edge_aspect_ratio -> Delete();
}

void HEX_T::write_hex_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, const std::vector<int> &ien_array,
      const std::vector< DataVecStr<int> > &IOdata, const bool &isXML )
{
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // Generate the mesh and compute aspect ratios
  gen_hex_grid( grid_w, numpts, numcels, pt, ien_array );

  // We need to make sure there are no data in IOdata that have the same name
  std::vector<std::string> name_list {};
  for( auto data : IOdata ) name_list.push_back( data.get_name() );

  VEC_T::sort_unique_resize( name_list );

  SYS_T::print_exit_if( name_list.size() != IOdata.size(),
    "Error: In HEX_T::write_hex_grid, there are %d data in the IOdata that have the same name.\n", IOdata.size() - name_list.size() + 1 );
  
  // We need to assign GlobalNodeID if the user does not provide it explicitly
  // in IOdata
  if( !VEC_T::is_invec(name_list, static_cast<std::string>("GlobalNodeID")) )
  {
    std::vector<int> node_idx(numpts, -1);
    for(int ii=0; ii<numpts; ++ii) node_idx[ii] = ii;

    VTK_T::add_int_PointData( grid_w, node_idx, "GlobalNodeID" );
    VEC_T::clean( node_idx );
  }

  // We need to assign GlobalElementID if the user does not provide it
  // explicitly in IOdata
  if( !VEC_T::is_invec(name_list, static_cast<std::string>("GlobalElementID")) )
  {
    std::vector<int> elem_idx(numcels, -1);
    for(int ii=0; ii<numcels; ++ii) elem_idx[ii] = ii;

    VTK_T::add_int_CellData( grid_w, elem_idx, "GlobalElementID" );
    VEC_T::clean( elem_idx );
  }

  // We add the IOdata for VTK
  for( auto data : IOdata )
  {
    if( data.get_object() == AssociateObject::Node )
      VTK_T::add_int_PointData( grid_w, data.get_data(), data.get_name() );
    else if( data.get_object() == AssociateObject::Cell )
      VTK_T::add_int_CellData( grid_w, data.get_data(), data.get_name() );
    else
      SYS_T::print_exit( "Error: In HEX_T::write_hex_grid, there is an unknown object type in DataVecStr %s", data.get_name().c_str() );
  }

  // Write the prepared grid_w to a vtu or vtk file on disk
  VTK_T::write_vtkPointSet(filename, grid_w, isXML);

  grid_w->Delete();
}

void HEX_T::gen_quad_grid( vtkPolyData * const &grid_w,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array )
{
  // check the input data compatibility
  SYS_T::print_exit_if(pt.size() != 3*numpts, 
    "Error: HEX_T::gen_quad_grid, point vector size does not match the number of points. \n");

  SYS_T::print_exit_if(ien_array.size() != 4*numcels,
    "Error: HEX_T::gen_quad_grid, ien array size does not match the number of cells. \n");

  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New();
  ppt->SetDataTypeToDouble();
  for(int ii=0; ii<numpts; ++ii)
  {
    const double coor[3] = { pt[3*ii], pt[3*ii+1], pt[3*ii+2] };
    ppt -> InsertPoint(ii, coor);
  }

  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  // 2. Cell
  vtkCellArray * cl = vtkCellArray::New();
  for(int ii=0; ii<numcels; ++ii)
  {
    vtkQuad * tr = vtkQuad::New();

    tr->GetPointIds()->SetId( 0, ien_array[4*ii] );
    tr->GetPointIds()->SetId( 1, ien_array[4*ii+1] );
    tr->GetPointIds()->SetId( 2, ien_array[4*ii+2] );
    tr->GetPointIds()->SetId( 3, ien_array[4*ii+3] );
    cl -> InsertNextCell(tr);
    tr -> Delete();
  }
  grid_w->SetPolys(cl);
  cl->Delete();
}

void HEX_T::write_quad_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata )
{
  // Setup the VTK objects
  vtkPolyData * grid_w = vtkPolyData::New();

  // Generate the mesh
  gen_quad_grid( grid_w, numpts, numcels, pt, ien_array );

  // We need to make sure there are no data in IOdata that have the same name
  std::vector<std::string> name_list {};
  for( auto data : IOdata ) name_list.push_back( data.get_name() );

  VEC_T::sort_unique_resize( name_list );

  SYS_T::print_exit_if( name_list.size() != IOdata.size(),
    "Error: In HEX_T::write_quad_grid, there are %d data in the IOdata that have the same name.\n", IOdata.size() - name_list.size() + 1 );

// We add the IOdata for VTK
  for( auto data : IOdata )
  {
    if( data.get_object() == AssociateObject::Node )
      VTK_T::add_int_PointData( grid_w, data.get_data(), data.get_name() );
    else if( data.get_object() == AssociateObject::Cell )
      VTK_T::add_int_CellData( grid_w, data.get_data(), data.get_name() );
    else
      SYS_T::print_exit( "Error: In HEX_T::write_quad_grid, there is an unknown object type in DataVecStr %s", data.get_name().c_str() );
  }

  // Write the prepared grid_w to a vtu or vtk file on disk
  VTK_T::write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}

void HEX_T::gen_quadratic_quad_grid( vtkUnstructuredGrid * const &grid_w,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array )
{
  // check the input data compatibility
  SYS_T::print_exit_if(pt.size() != 3*numpts,
    "Error: HEX_T::gen_quadratic_quad_grid point, vector size does not match the number of points. \n");

  SYS_T::print_exit_if(ien_array.size() != 9*numcels, 
    "Error: HEX_T::gen_quadratic_quad_grid, ien array size does not match the number of cells. \n");

  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New(); 
  ppt->SetDataTypeToDouble();
  for(int ii=0; ii<numpts; ++ii)
  {
    const double coor[3] = { pt[3*ii], pt[3*ii+1], pt[3*ii+2] };
    ppt -> InsertPoint(ii, coor);
  }

  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  // 2. Cell
  vtkCellArray * cl = vtkCellArray::New();

  for(int ii=0; ii<numcels; ++ii)
  {
    vtkBiQuadraticQuad * tr = vtkBiQuadraticQuad::New();

    tr->GetPointIds()->SetId( 0, ien_array[9*ii] );
    tr->GetPointIds()->SetId( 1, ien_array[9*ii+1] );
    tr->GetPointIds()->SetId( 2, ien_array[9*ii+2] );
    tr->GetPointIds()->SetId( 3, ien_array[9*ii+3] );
    tr->GetPointIds()->SetId( 4, ien_array[9*ii+4] );
    tr->GetPointIds()->SetId( 5, ien_array[9*ii+5] );
    tr->GetPointIds()->SetId( 6, ien_array[9*ii+6] );
    tr->GetPointIds()->SetId( 7, ien_array[9*ii+7] );
    tr->GetPointIds()->SetId( 8, ien_array[9*ii+8] );
    cl -> InsertNextCell(tr);
    tr -> Delete();
  }

  grid_w -> SetCells(28, cl);
  cl -> Delete();
}

void HEX_T::write_quadratic_quad_grid( const std::string &filename,
      const int &numpts, const int &numcels,
      const std::vector<double> &pt, 
      const std::vector<int> &ien_array,
      const std::vector<DataVecStr<int>> &IOdata )
{
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // Generate the mesh
  gen_quadratic_quad_grid( grid_w, numpts, numcels, pt, ien_array );

  // We need to make sure there are no data in IOdata that have the same name
  std::vector<std::string> name_list {};
  for( auto data : IOdata ) name_list.push_back( data.get_name() );

  VEC_T::sort_unique_resize( name_list );

  SYS_T::print_exit_if( name_list.size() != IOdata.size(),
    "Error: In HEX_T::write_quadratic_quad_grid, there are %d data in the IOdata that have the same name.\n", IOdata.size() - name_list.size() + 1 );

// We add the IOdata for VTK
  for( auto data : IOdata )
  {
    if( data.get_object() == AssociateObject::Node )
      VTK_T::add_int_PointData( grid_w, data.get_data(), data.get_name() );
    else if( data.get_object() == AssociateObject::Cell )
      VTK_T::add_int_CellData( grid_w, data.get_data(), data.get_name() );
    else
      SYS_T::print_exit( "Error: In HEX_T::write_quadratic_quad_grid, there is an unknown object type in DataVecStr %s", data.get_name().c_str() );
  }

  // Write the prepared grid_w to a vtu file on disk
  VTK_T::write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}

double HEX_T::get_aspect_ratio( const std::vector<double> &pt )
{
  double edge[12];
  //e01
  edge[0] = MATH_T::norm2(pt[3]-pt[0], pt[4]-pt[1], pt[5]-pt[2]);
  //e12
  edge[1] = MATH_T::norm2(pt[6]-pt[3], pt[7]-pt[4], pt[8]-pt[5]);
  //e23
  edge[2] = MATH_T::norm2(pt[9]-pt[6], pt[10]-pt[7], pt[11]-pt[8]);
  //e03
  edge[3] = MATH_T::norm2(pt[9]-pt[0], pt[10]-pt[1], pt[11]-pt[2]);
  //e45
  edge[4] = MATH_T::norm2(pt[15]-pt[12], pt[16]-pt[13], pt[17]-pt[14]);
  //e56
  edge[5] = MATH_T::norm2(pt[18]-pt[15], pt[19]-pt[16], pt[20]-pt[17]);
  //e67
  edge[6] = MATH_T::norm2(pt[21]-pt[18], pt[22]-pt[19], pt[23]-pt[20]);
  //e47
  edge[7] = MATH_T::norm2(pt[21]-pt[12], pt[22]-pt[13], pt[23]-pt[14]);
  //e04
  edge[8] = MATH_T::norm2(pt[12]-pt[0], pt[13]-pt[1], pt[14]-pt[2]);
  //e15
  edge[9] = MATH_T::norm2(pt[15]-pt[3], pt[16]-pt[4], pt[17]-pt[5]);
  //e26
  edge[10] = MATH_T::norm2(pt[18]-pt[6], pt[19]-pt[7], pt[20]-pt[8]);
  //e37
  edge[11] = MATH_T::norm2(pt[21]-pt[9], pt[22]-pt[10], pt[23]-pt[11]);

  const double emax = *std::max_element(edge, edge+12);
  const double emin = *std::min_element(edge, edge+12);

  return emax / emin;
}

Vector_3 HEX_T::get_out_normal( const std::string &file,
      const std::vector<double> &vol_ctrlPts,
      const IIEN * const &vol_ien )
{
  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien;

  // Analyze the file type
  std::string fend; fend.assign( file.end()-4 , file.end() );

  VTK_T::read_grid( file, numpts, numcels, pts, ien);

  const std::vector<int> gnode = VTK_T::read_int_PointData(file, "GlobalNodeID");
  const std::vector<int> gelem = VTK_T::read_int_CellData(file, "GlobalElementID");

  // quad nodes' global indices
  const std::vector<int> qun = { gnode[ ien[0] ], gnode[ ien[1] ], gnode[ ien[2] ], gnode[ ien[3] ] };

  const int hexe0 = gelem[0]; // quad's associated hex element indices

  SYS_T::print_exit_if(hexe0 == -1, "Error: HEX_T::get_out_normal the vtp file %s contains -1 for GlobalElementID.\n", file.c_str());

  const std::vector<int> hen { vol_ien->get_IEN(hexe0, 0), vol_ien->get_IEN(hexe0, 1), 
                               vol_ien->get_IEN(hexe0, 2), vol_ien->get_IEN(hexe0, 3),
                               vol_ien->get_IEN(hexe0, 4), vol_ien->get_IEN(hexe0, 5),
                               vol_ien->get_IEN(hexe0, 6), vol_ien->get_IEN(hexe0, 7)};

  // the coordinates of the center of the hex element, defined as the average of the coordinates of 8 vertices
  Vector_3 cen_hex( 0.0, 0.0, 0.0 );
  for (int vertex : hen)
    cen_hex += 0.125 * Vector_3( vol_ctrlPts[3 * vertex], vol_ctrlPts[3 * vertex + 1], vol_ctrlPts[3 * vertex + 2] );

  // the coordinates of the center of the surface quad element, defined as the average of the coordinates of 4 vertices
  Vector_3 cen_quad(0.0, 0.0, 0.0);
  for (int vertex : qun)
    cen_quad += 0.25 * Vector_3( vol_ctrlPts[3 * vertex], vol_ctrlPts[3 * vertex + 1], vol_ctrlPts[3 * vertex + 2] );

  // obtain the inward vector
  const Vector_3 inw = cen_hex - cen_quad;

  // make cross line-0-1 and line-0-3
  const Vector_3 l01 ( vol_ctrlPts[3*qun[1]    ] - vol_ctrlPts[3*qun[0]    ],
                       vol_ctrlPts[3*qun[1] + 1] - vol_ctrlPts[3*qun[0] + 1],
                       vol_ctrlPts[3*qun[1] + 2] - vol_ctrlPts[3*qun[0] + 2] );
  
  const Vector_3 l03 ( vol_ctrlPts[3*qun[3]    ] - vol_ctrlPts[3*qun[0]    ],
                       vol_ctrlPts[3*qun[3] + 1] - vol_ctrlPts[3*qun[0] + 1],
                       vol_ctrlPts[3*qun[3] + 2] - vol_ctrlPts[3*qun[0] + 2] );

  Vector_3 outVec = cross_product( l01, l03 );

  outVec.normalize();

  // inner product outward with inward, and correct outVec
  if(inw.dot_product(outVec) > 0.0) outVec *= -1.0;

  return outVec;
}

std::vector<int> HEX_T::reset_node (const std::vector<int> &ien)
{
  if (ien.size() == 27)
  {
   std::vector<int> local_ien_vtk {ien[0], ien[1], ien[2], ien[3], ien[4], ien[5],
     ien[6], ien[7], ien[8], ien[11], ien[13], ien[9], ien[16], ien[18], ien[19],
     ien[17], ien[10], ien[12], ien[14], ien[15], ien[22], ien[23], ien[21], ien[24],
     ien[20], ien[25], ien[26]};
   return local_ien_vtk;
  }
  else if (ien.size() == 20)
  {
    std::vector<int> local_ien_vtk {ien[0], ien[1], ien[2], ien[3], ien[4], ien[5],
     ien[6], ien[7], ien[8], ien[11], ien[13], ien[9], ien[16], ien[18], ien[19],
     ien[17], ien[10], ien[12], ien[14], ien[15]};
   return local_ien_vtk;
  }
  else
    SYS_T::print_exit("Error: In HEX_T::reset_node, undefined hexahedron type.");
}

namespace HEX_T
{
  Hex8::Hex8()
  {
    pts.resize(24);
    pts[0]  = -1.0; pts[1]  = -1.0; pts[2]  = -1.0;
    pts[3]  =  1.0; pts[4]  = -1.0; pts[5]  = -1.0;
    pts[6]  =  1.0; pts[7]  =  1.0; pts[8]  = -1.0;
    pts[9]  = -1.0; pts[10] =  1.0; pts[11] = -1.0;
    pts[12] = -1.0; pts[13] = -1.0; pts[14] =  1.0;
    pts[15] =  1.0; pts[16] = -1.0; pts[17] =  1.0;
    pts[18] =  1.0; pts[19] =  1.0; pts[20] =  1.0;
    pts[21] = -1.0; pts[22] =  1.0; pts[23] =  1.0;

    gindex[0] = 0; gindex[1] = 1;
    gindex[2] = 2; gindex[3] = 3;
    gindex[4] = 4; gindex[5] = 5;
    gindex[6] = 6; gindex[7] = 7;
  }

  Hex8::Hex8( const std::vector<double> &in_nodes )
  {
    pts.resize(24);
    SYS_T::print_exit_if( in_nodes.size() != 24,
      "Error: input nodal list shall have 8 nodes with xyz-coordinates. \n");
    
    for(int ii=0; ii<24; ++ii) pts[ii] = in_nodes[ii];

    gindex[0] = 0; gindex[1] = 1;
    gindex[2] = 2; gindex[3] = 3;
    gindex[4] = 4; gindex[5] = 5;
    gindex[6] = 6; gindex[7] = 7;
  }

  Hex8::Hex8( const std::vector<double> &ctrlPts,
          const int &ien0, const int &ien1, const int &ien2,
          const int &ien3, const int &ien4, const int &ien5,
          const int &ien6, const int &ien7 )
  {
    pts.resize(24);
    gindex[0] = ien0; gindex[1] = ien1;
    gindex[2] = ien2; gindex[3] = ien3;
    gindex[4] = ien4; gindex[5] = ien5;
    gindex[6] = ien6; gindex[7] = ien7;

    for(int ii=0; ii<8; ++ii)
    {
      pts[ii*3]   = ctrlPts[gindex[ii]*3];
      pts[ii*3+1] = ctrlPts[gindex[ii]*3+1];
      pts[ii*3+2] = ctrlPts[gindex[ii]*3+2];
    }
  }

  Hex8::~Hex8()
  {
  }

  void Hex8::reset( const std::vector<double> &in_nodes )
  {
    SYS_T::print_exit_if( in_nodes.size() != 24,
      "Error: input nodal list shall have 8 nodes with xyz-coordinates. \n");
    
    for(int ii=0; ii<24; ++ii) pts[ii] = in_nodes[ii];

    gindex[0] = 0; gindex[1] = 1;
    gindex[2] = 2; gindex[3] = 3;
    gindex[4] = 4; gindex[5] = 5;
    gindex[6] = 6; gindex[7] = 7;
  }

  void Hex8::reset( const std::vector<double> &ctrlPts,
          const int &ien0, const int &ien1, const int &ien2,
          const int &ien3, const int &ien4, const int &ien5,
          const int &ien6, const int &ien7 )
  {
    gindex[0] = ien0; gindex[1] = ien1;
    gindex[2] = ien2; gindex[3] = ien3;
    gindex[4] = ien4; gindex[5] = ien5;
    gindex[6] = ien6; gindex[7] = ien7;

    for(int ii=0; ii<8; ++ii)
    {
      pts[ii*3]   = ctrlPts[gindex[ii]*3];
      pts[ii*3+1] = ctrlPts[gindex[ii]*3+1];
      pts[ii*3+2] = ctrlPts[gindex[ii]*3+2];
    }
  }

  void Hex8::reset( const int &ien0, const int &ien1,
          const int &ien2, const int &ien3, const int &ien4,
          const int &ien5, const int &ien6, const int &ien7 )
  {
    gindex[0] = ien0; gindex[1] = ien1;
    gindex[2] = ien2; gindex[3] = ien3;
    gindex[4] = ien4; gindex[5] = ien5;
    gindex[6] = ien6; gindex[7] = ien7;
  }

  void Hex8::reset( const std::vector<double> &ctrlPts,
          const IIEN * const &ien_ptr, const int &ee )
  {
    gindex[0] = ien_ptr->get_IEN(ee, 0);
    gindex[1] = ien_ptr->get_IEN(ee, 1);
    gindex[2] = ien_ptr->get_IEN(ee, 2);
    gindex[3] = ien_ptr->get_IEN(ee, 3);
    gindex[4] = ien_ptr->get_IEN(ee, 4);
    gindex[5] = ien_ptr->get_IEN(ee, 5);
    gindex[6] = ien_ptr->get_IEN(ee, 6);
    gindex[7] = ien_ptr->get_IEN(ee, 7);

    for(int ii=0; ii<8; ++ii)
    {
      pts[ii*3]   = ctrlPts[gindex[ii]*3];
      pts[ii*3+1] = ctrlPts[gindex[ii]*3+1];
      pts[ii*3+2] = ctrlPts[gindex[ii]*3+2];
    }
  }

  void Hex8::print_info() const
  {
    std::cout<<"Hex8 object : \n";
    std::cout<<" -- global indices: "<<'\n'<<gindex[0]<<'\t';
    std::cout<<gindex[1]<<'\t'<<gindex[2]<<'\t'<<gindex[3]<<'\t';
    std::cout<<gindex[4]<<'\t'<<gindex[5]<<'\t'<<gindex[6]<<'\t'<<gindex[7]<<'\n';
    std::cout<<" -- pt0 : "<<pts[0]<<'\t'<<pts[1]<<'\t'<<pts[2]<<'\n';
    std::cout<<" -- pt1 : "<<pts[3]<<'\t'<<pts[4]<<'\t'<<pts[5]<<'\n';
    std::cout<<" -- pt2 : "<<pts[6]<<'\t'<<pts[7]<<'\t'<<pts[8]<<'\n';
    std::cout<<" -- pt3 : "<<pts[9]<<'\t'<<pts[10]<<'\t'<<pts[11]<<'\n';
    std::cout<<" -- pt4 : "<<pts[12]<<'\t'<<pts[13]<<'\t'<<pts[14]<<'\n';
    std::cout<<" -- pt5 : "<<pts[15]<<'\t'<<pts[16]<<'\t'<<pts[17]<<'\n';
    std::cout<<" -- pt6 : "<<pts[18]<<'\t'<<pts[19]<<'\t'<<pts[20]<<'\n';
    std::cout<<" -- pt7 : "<<pts[21]<<'\t'<<pts[22]<<'\t'<<pts[23]<<std::endl;
  }

  double Hex8::get_aspect_ratio() const
  {
    double edge[12];
    //e01
    edge[0] = MATH_T::norm2(pts[3]-pts[0], pts[4]-pts[1], pts[5]-pts[2]);
    //e12
    edge[1] = MATH_T::norm2(pts[6]-pts[3], pts[7]-pts[4], pts[8]-pts[5]);
    //e23
    edge[2] = MATH_T::norm2(pts[9]-pts[6], pts[10]-pts[7], pts[11]-pts[8]);
    //e03
    edge[3] = MATH_T::norm2(pts[9]-pts[0], pts[10]-pts[1], pts[11]-pts[2]);
    //e45
    edge[4] = MATH_T::norm2(pts[15]-pts[12], pts[16]-pts[13], pts[17]-pts[14]);
    //e56
    edge[5] = MATH_T::norm2(pts[18]-pts[15], pts[19]-pts[16], pts[20]-pts[17]);
    //e67
    edge[6] = MATH_T::norm2(pts[21]-pts[18], pts[22]-pts[19], pts[23]-pts[20]);
    //e47
    edge[7] = MATH_T::norm2(pts[21]-pts[12], pts[22]-pts[13], pts[23]-pts[14]);
    //e04
    edge[8] = MATH_T::norm2(pts[12]-pts[0], pts[13]-pts[1], pts[14]-pts[2]);
    //e15
    edge[9] = MATH_T::norm2(pts[15]-pts[3], pts[16]-pts[4], pts[17]-pts[5]);
    //e26
    edge[10] = MATH_T::norm2(pts[18]-pts[6], pts[19]-pts[7], pts[20]-pts[8]);
    //e37
    edge[11] = MATH_T::norm2(pts[21]-pts[9], pts[22]-pts[10], pts[23]-pts[11]);

    const double emax = *std::max_element(edge, edge+12);
    const double emin = *std::min_element(edge, edge+12);

    return emax / emin;
  }

  double Hex8::get_volume() const
  {
    TET_T::Tet4 tet1(pts, 0, 1, 3, 7);
    TET_T::Tet4 tet2(pts, 0, 4, 1, 7);
    TET_T::Tet4 tet3(pts, 1, 4, 5, 7);
    TET_T::Tet4 tet4(pts, 1, 2, 3, 7);
    TET_T::Tet4 tet5(pts, 1, 5, 6, 7);
    TET_T::Tet4 tet6(pts, 1, 6, 2, 7);

    return tet1.get_volume() + tet2.get_volume() + tet3.get_volume() + tet4.get_volume() + tet5.get_volume() + tet6.get_volume();
  }

}

// EOF
