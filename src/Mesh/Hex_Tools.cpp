#include "Hex_Tools.hpp"

void HEX_T::gen_hex_grid( vtkUnstructuredGrid * const &grid_w,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array )
{
  // Check the input data compatibility
  SYS_T::print_exit_if( int(pt.size()) != 3*numpts,
    "Error: HEX_T::gen_hex_grid point vector size does not match the number of points. \n" );

  // detect the element type
  int nlocbas {-1};
  if ( int(ien_array.size()) == 8*numcels )
    nlocbas = 8;
  else if ( int(ien_array.size()) == 27*numcels ) 
    nlocbas = 27;
  else
    SYS_T::print_exit("Error: HEX_T::gen_hex_grid ien array size does not match the number of cells. \n");
  // Now we donnot consider the serendipity elements i.e. 8-node quadrangle and 20-node hexahedron.

  // Check the connectivity array
  std::vector<int> temp = ien_array;
  VEC_T::sort_unique_resize(temp);
  SYS_T::print_exit_if( int(temp.size()) != numpts,
  "Error: HEX_T::gen_hex_grid numpts does not match the number of unique points in the ien array. Please re-organize the input. \n" );
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
    vtkCell * cl = vtkQuadraticHexahedron::New();

    for(int ii=0; ii<numcels; ++ii)
    {
      for(int lnode=0; lnode<27; ++lnode)
        cl->GetPointIds()->SetId( lnode, ien_array[27*ii + lnode] );

      grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

      std::vector<double> cell_node; cell_node.clear();
      for(int lnode=0; lnode<8; ++lnode)
      {
        int node_offset = 3 * ien_array[27*ii + lnode];
        cell_node.push_back(pt[node_offset]);
        cell_node.push_back(pt[node_offset+1]);
        cell_node.push_back(pt[node_offset+2]);
      }
      edge_aspect_ratio -> InsertNextValue( HEX_T::get_aspect_ratio(cell_node) );
    }

    cl -> Delete();
  }
  else
    SYS_T::print_exit("Error: HEX_T::gen_hex_grid unknown local basis number.\n");

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