#include "VTK_Tools.hpp"

void VTK_T::read_vtu_grid( const std::string &filename,
    int &numpts, int &numcels,
    std::vector<double> &pt, std::vector<int> &ien_array )
{
  vtkXMLUnstructuredGridReader * reader = vtkXMLUnstructuredGridReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  reader -> GlobalWarningDisplayOff();
  vtkUnstructuredGrid * vtkugrid = reader -> GetOutput();

  // Number of grid points in the mesh
  numpts  = static_cast<int>( vtkugrid -> GetNumberOfPoints() );
  
  SYS_T::print_fatal_if(numpts <= 0, "Error: the file %s contains no point. \n", filename.c_str());
  
  // Number of cells in the mesh
  numcels = static_cast<int>( vtkugrid -> GetNumberOfCells() );

  SYS_T::print_fatal_if(numcels <= 0, "Error: the file %s contains no cell. \n", filename.c_str());
  
  // xyz coordinates of the points
  pt.clear();
  for(int ii=0; ii<numpts; ++ii)
  {
    double pt_xyz[3];
    vtkugrid -> GetPoint(ii, pt_xyz);
    pt.push_back(pt_xyz[0]);
    pt.push_back(pt_xyz[1]);
    pt.push_back(pt_xyz[2]);
  }
 
  // Connectivity of the mesh 
  ien_array.clear();
  
  for(int ii=0; ii<numcels; ++ii)
  {
    vtkCell * cell = vtkugrid -> GetCell(ii);

    if( cell->GetCellType() == 10 ) 
    {
      // cell type 10 is four-node tet
      ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(3) ) );
    }
    else if( cell-> GetCellType() == 22 )
    {
      // cell type 22 is six-node triangle
      ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(3) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(4) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(5) ) );
    }
    else if( cell-> GetCellType() == 24 )
    {
      // cell type 24 is ten-node tet
      ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(3) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(4) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(5) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(6) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(7) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(8) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(9) ) );
    }
    else SYS_T::print_fatal("Error: VTK_T::read_vtu_grid read a mesh with VTK cell type 10, 24, or 22. \n"); 
  }

  reader->Delete();
}


std::vector<int> VTK_T::read_int_CellData( const std::string &filename,
    const std::string &dataname )
{
  vtkXMLGenericDataObjectReader * reader = vtkXMLGenericDataObjectReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  
  vtkCellData * celldata = nullptr;
  int numcels = -1;

  // Downcasting will return null if fails
  if(dynamic_cast<vtkPolyData*>(reader->GetOutput()))
  {
    vtkPolyData * vtkgrid = reader -> GetPolyDataOutput (); 
    celldata = vtkgrid->GetCellData();
    numcels = static_cast<int>( vtkgrid -> GetNumberOfCells() );
  }
  else if(dynamic_cast<vtkUnstructuredGrid*>(reader->GetOutput()))
  {
    vtkUnstructuredGrid * vtkgrid = reader -> GetUnstructuredGridOutput();
    celldata = vtkgrid->GetCellData();
    numcels = static_cast<int>( vtkgrid -> GetNumberOfCells() );
  }
  else
    SYS_T::print_fatal("VTK_T::read_int_CellData unknown vtk object type.\n");

  vtkDataArray * cd = celldata->GetScalars( dataname.c_str() );

  std::vector<int> data( numcels );
  for(int ii=0; ii<numcels; ++ii)
    data[ii] = static_cast<int>( cd->GetComponent(ii, 0) );

  reader -> Delete();

  return data;
}


std::vector<double> VTK_T::read_double_CellData( const std::string &filename,
    const std::string &dataname )
{
  vtkXMLGenericDataObjectReader * reader = vtkXMLGenericDataObjectReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  
  vtkCellData * celldata = nullptr;
  int numcels = -1;

  // Downcasting will return null if fails
  if(dynamic_cast<vtkPolyData*>(reader->GetOutput()))
  {
    vtkPolyData * vtkgrid = reader -> GetPolyDataOutput (); 
    celldata = vtkgrid->GetCellData();
    numcels = static_cast<int>( vtkgrid -> GetNumberOfCells() );
  }
  else if(dynamic_cast<vtkUnstructuredGrid*>(reader->GetOutput()))
  {
    vtkUnstructuredGrid * vtkgrid = reader -> GetUnstructuredGridOutput();
    celldata = vtkgrid->GetCellData();
    numcels = static_cast<int>( vtkgrid -> GetNumberOfCells() );
  }
  else
    SYS_T::print_fatal("VTK_T::read_double_CellData unknown vtk object type.\n");

  vtkDataArray * cd = celldata->GetScalars( dataname.c_str() );

  std::vector<double> data( numcels );
  for(int ii=0; ii<numcels; ++ii)
    data[ii] = static_cast<double>( cd->GetComponent(ii, 0) );

  reader -> Delete();

  return data;
}


std::vector<int> VTK_T::read_int_PointData( const std::string &filename,
    const std::string &dataname )
{
  vtkXMLGenericDataObjectReader * reader = vtkXMLGenericDataObjectReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  
  vtkPointData * pointdata = nullptr;
  int numpts = -1;

  // Downcasting will return null if fails
  if(dynamic_cast<vtkPolyData*>(reader->GetOutput()))
  {
    vtkPolyData * vtkgrid = reader -> GetPolyDataOutput (); 
    pointdata = vtkgrid->GetPointData();
    numpts = static_cast<int>( vtkgrid -> GetNumberOfPoints() );
  }
  else if(dynamic_cast<vtkUnstructuredGrid*>(reader->GetOutput()))
  {
    vtkUnstructuredGrid * vtkgrid = reader -> GetUnstructuredGridOutput();
    pointdata = vtkgrid->GetPointData();
    numpts = static_cast<int>( vtkgrid -> GetNumberOfPoints() );
  }
  else
    SYS_T::print_fatal("VTK_T::read_int_PointData unknown vtk object type.\n");

  vtkDataArray * pd = pointdata->GetScalars( dataname.c_str() );

  std::vector<int> data( numpts );
  for(int ii=0; ii<numpts; ++ii)
    data[ii] = static_cast<int>( pd->GetComponent(ii, 0) );

  reader -> Delete();

  return data;
}


std::vector<double> VTK_T::read_double_PointData( const std::string &filename,
    const std::string &dataname )
{
  vtkXMLGenericDataObjectReader * reader = vtkXMLGenericDataObjectReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  
  vtkPointData * pointdata = nullptr;
  int numpts = -1;

  // Downcasting will return null if fails
  if(dynamic_cast<vtkPolyData*>(reader->GetOutput()))
  {
    vtkPolyData * vtkgrid = reader -> GetPolyDataOutput (); 
    pointdata = vtkgrid->GetPointData();
    numpts = static_cast<int>( vtkgrid -> GetNumberOfPoints() );
  }
  else if(dynamic_cast<vtkUnstructuredGrid*>(reader->GetOutput()))
  {
    vtkUnstructuredGrid * vtkgrid = reader -> GetUnstructuredGridOutput();
    pointdata = vtkgrid->GetPointData();
    numpts = static_cast<int>( vtkgrid -> GetNumberOfPoints() );
  }
  else
    SYS_T::print_fatal("VTK_T::read_double_PointData unknown vtk object type.\n");

  vtkDataArray * pd = pointdata->GetScalars( dataname.c_str() );

  std::vector<double> data( numpts );
  for(int ii=0; ii<numpts; ++ii)
    data[ii] = static_cast<double>( pd->GetComponent(ii, 0) );

  reader -> Delete();

  return data;
}


void VTK_T::read_vtu_grid( const std::string &filename,
    int &numpts, int &numcels,
    std::vector<double> &pt, std::vector<int> &ien_array,
    std::vector<int> &phy_tag )
{
  read_vtu_grid(filename, numpts, numcels, pt, ien_array);
  
  phy_tag = read_int_CellData(filename, "Physics_tag");
}


void VTK_T::read_vtu_grid( const std::string &filename,
    int &numpts, int &numcels,
    std::vector<double> &pt, std::vector<int> &ien_array,
    std::vector<int> &global_node_index,
    std::vector<int> &global_elem_index )
{
  read_vtu_grid(filename, numpts, numcels, pt, ien_array);
  
  global_node_index = read_int_PointData(filename, "GlobalNodeID"); 

  global_elem_index = read_int_CellData(filename, "GlobalElementID");
}


void VTK_T::read_vtp_grid( const std::string &filename,
    int &numpts, int &numcels,
    std::vector<double> &pt, std::vector<int> &ien_array )
{
  vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  vtkPolyData * polydata = reader -> GetOutput();
  
  // Number of grid points in the mesh
  numpts = static_cast<int>( polydata -> GetNumberOfPoints() );
  
  SYS_T::print_fatal_if(numpts <= 0, "Error: the file %s contains no point. \n", filename.c_str());
  
  // Number of cells in the mesh
  numcels = static_cast<int>( polydata -> GetNumberOfPolys() );

  SYS_T::print_fatal_if(numcels <= 0, "Error: the file %s contains no cell. \n", filename.c_str());
  
  // xyz coordinates of the points
  pt.clear();
  for(int ii=0; ii<numpts; ++ii)
  {
    double pt_xyz[3];
    polydata -> GetPoint(ii, pt_xyz);
    pt.push_back(pt_xyz[0]);
    pt.push_back(pt_xyz[1]);
    pt.push_back(pt_xyz[2]);
  }

  // xyz coordinates of the points
  ien_array.clear();
  for(int ii=0; ii<numcels; ++ii)
  {
    vtkCell * cell = polydata -> GetCell(ii);

    if( cell->GetCellType() == 5 )
    {
      ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
    }
    else SYS_T::print_fatal("Error: read_vtp_grid read a mesh with non triangle elements. \n");
  }

  reader->Delete();
}


void VTK_T::read_vtp_grid( const std::string &filename,
    int &numpts, int &numcels,
    std::vector<double> &pt, std::vector<int> &ien_array,
    std::vector<int> &global_node_index,
    std::vector<int> &global_elem_index )
{
  read_vtp_grid(filename, numpts, numcels, pt, ien_array);
  
  global_node_index = read_int_PointData(filename, "GlobalNodeID"); 

  global_elem_index = read_int_CellData(filename, "GlobalElementID");
}


void VTK_T::gen_tet_grid( vtkUnstructuredGrid * const &grid_w,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array )
{
  // Check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: VTK_T::write_tet_grid point vector size does not match the number of points. \n");

  // detect the element type
  int nlocbas = -1;
  if( int(ien_array.size()) == 4*numcels ) nlocbas = 4;
  else if( int(ien_array.size()) == 10*numcels ) nlocbas = 10;
  else SYS_T::print_fatal("Error: VTK_T::write_tet_grid ien array size does not match the number of cells. \n");

  // Check the connectivity array
  std::vector<int> temp = ien_array;
  VEC_T::sort_unique_resize(temp);
  if( int(temp.size()) != numpts ) SYS_T::print_fatal("Error: VTK_T::write_tet_grid numpts does not match the number of unique points in the ien array. Please re-organize the input. \n");
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

  if(nlocbas == 4) 
  {
    vtkCell * cl = vtkTetra::New();

    for(int ii=0; ii<numcels; ++ii)
    {
      for(int lnode=0; lnode<4; ++lnode)
        cl->GetPointIds()->SetId( lnode, ien_array[4*ii + lnode] );

      grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

      std::vector<double> cell_node; cell_node.clear();
      for(int lnode=0; lnode<4; ++lnode)
      {
        int node_offset = 3 * ien_array[4*ii + lnode];
        cell_node.push_back(pt[node_offset]);
        cell_node.push_back(pt[node_offset+1]);
        cell_node.push_back(pt[node_offset+2]);
      }
      edge_aspect_ratio -> InsertNextValue( VTK_T::get_aspect_ratio(cell_node) );
    }

    cl -> Delete();
  }
  else if(nlocbas == 10)
  {
    vtkCell * cl = vtkQuadraticTetra::New();

    for(int ii=0; ii<numcels; ++ii)
    {
      for(int lnode=0; lnode<10; ++lnode)
        cl->GetPointIds()->SetId( lnode, ien_array[10*ii + lnode] );

      grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

      std::vector<double> cell_node; cell_node.clear();
      for(int lnode=0; lnode<4; ++lnode)
      {
        int node_offset = 3 * ien_array[10*ii + lnode];
        cell_node.push_back(pt[node_offset]);
        cell_node.push_back(pt[node_offset+1]);
        cell_node.push_back(pt[node_offset+2]);
      }
      edge_aspect_ratio -> InsertNextValue( VTK_T::get_aspect_ratio(cell_node) );
    }

    cl -> Delete();
  }
  else SYS_T::print_fatal("Error: VTK_T::write_tet_grid unknown local basis number.\n");

  // Add the asepct-ratio to grid_w
  grid_w -> GetCellData() -> AddArray( edge_aspect_ratio );
  edge_aspect_ratio -> Delete();
}


void VTK_T::write_tet_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt, const std::vector<int> &ien_array )
{
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // Generate the mesh and compute aspect ratios
  gen_tet_grid( grid_w, numpts, numcels, pt, ien_array );

  // nodal indices (natural numbering)
  std::vector<int> node_idx(numpts);
  for(int ii=0; ii<numpts; ++ii) node_idx[ii] = ii;

  add_int_PointData( grid_w, node_idx, "GlobalNodeID" );

  // cell indices (natural numbering)
  std::vector<int> elem_idx(numcels);
  for(int ii=0; ii<numcels; ++ii) elem_idx[ii] = ii;

  add_int_CellData( grid_w, elem_idx, "GlobalElementID" );

  // write vtu (by default of the writer function)
  write_vtkPointSet(filename, grid_w);

  VEC_T::clean( node_idx );
  VEC_T::clean( elem_idx );

  grid_w->Delete();
}


void VTK_T::write_tet_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt, const std::vector<int> &ien_array,
    const std::vector<int> &node_idx, const std::vector<int> &elem_idx )
{
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // Generate the mesh and compute aspect ratios
  gen_tet_grid( grid_w, numpts, numcels, pt, ien_array );

  // nodal indices
  add_int_PointData( grid_w, node_idx, "GlobalNodeID" );

  // cell indices
  add_int_CellData( grid_w, elem_idx, "GlobalElementID" );

  // write vtu (by default of the writer function)
  write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}


void VTK_T::write_tet_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt, const std::vector<int> &ien_array,
    const std::vector<int> &phytag, const bool &isXML,
    const int &start_cell_index )
{
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // Generate the mesh and compute aspect ratios
  gen_tet_grid( grid_w, numpts, numcels, pt, ien_array );

  // nodal indices (natural numbering)
  std::vector<int> node_idx(numpts);
  for(int ii=0; ii<numpts; ++ii) node_idx[ii] = ii;

  add_int_PointData( grid_w, node_idx, "GlobalNodeID" );

  // cell indices (natural numbering)
  std::vector<int> elem_idx(numcels);
  for(int ii=0; ii<numcels; ++ii) elem_idx[ii] = ii + start_cell_index;

  add_int_CellData( grid_w, elem_idx, "GlobalElementID" );

  // physics tags
  add_int_CellData( grid_w, phytag, "Physics_tag" );

  // write vtu or vtk
  write_vtkPointSet(filename, grid_w, isXML);

  VEC_T::clean( node_idx );
  VEC_T::clean( elem_idx );

  grid_w->Delete();
}


void VTK_T::write_tet_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt, const std::vector<int> &ien_array,
    const std::vector<int> &node_idx, const std::vector<int> &elem_idx, 
    const std::vector<int> &phytag, const bool &isXML )
{
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // Generate the mesh and compute aspect ratios
  gen_tet_grid( grid_w, numpts, numcels, pt, ien_array );

  // nodal indices
  add_int_PointData( grid_w, node_idx, "GlobalNodeID" );

  // cell indices
  add_int_CellData( grid_w, elem_idx, "GlobalElementID" );

  // physics tags
  add_int_CellData( grid_w, phytag, "Physics_tag" );

  // write vtu or vtk
  write_vtkPointSet(filename, grid_w, isXML);

  grid_w->Delete();
}


void VTK_T::write_triangle_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<int> &node_index,
    const std::vector<int> &ele_index )
{
  // Setup the VTK objects
  vtkPolyData * grid_w = vtkPolyData::New();

  // Generate the mesh
  gen_triangle_grid( grid_w, numpts, numcels, pt, ien_array );

  // nodal indices
  add_int_PointData( grid_w, node_index, "GlobalNodeID" );

  // cell indices
  add_int_CellData( grid_w, ele_index, "GlobalElementID");

  // write grid_w to vtp file
  write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}


void VTK_T::gen_triangle_grid( vtkPolyData * const &grid_w,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array )
{
  // check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: VTK_T::gen_triangle_grid point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 3*numcels) SYS_T::print_fatal("Error: VTK_T::gen_triangle_grid ien array size does not match the number of cells. \n");

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
    vtkTriangle * tr = vtkTriangle::New();

    tr->GetPointIds()->SetId( 0, ien_array[3*ii] );
    tr->GetPointIds()->SetId( 1, ien_array[3*ii+1] );
    tr->GetPointIds()->SetId( 2, ien_array[3*ii+2] );
    cl -> InsertNextCell(tr);
    tr -> Delete();
  }
  grid_w->SetPolys(cl);
  cl->Delete();
}


void VTK_T::write_quadratic_triangle_grid( 
    const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<int> &node_index,
    const std::vector<int> &ele_index )
{
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // Generate the mesh
  gen_quadratic_triangle_grid( grid_w, numpts, numcels, pt, ien_array );

  // nodal indices
  add_int_PointData( grid_w, node_index, "GlobalNodeID" );

  // cell indices
  add_int_CellData( grid_w, ele_index, "GlobalElementID" );

  // write vtu (by default of the writer function)
  write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}


void VTK_T::gen_quadratic_triangle_grid( vtkUnstructuredGrid * const &grid_w,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array )
{
  // check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: VTK_T::gen_quadratic_triangle_grid point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 6*numcels) SYS_T::print_fatal("Error: VTK_T::gen_quadratic_quadratic_triangle_grid ien array size does not match the number of cells. \n");

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
    vtkQuadraticTriangle * tr = vtkQuadraticTriangle::New();

    tr->GetPointIds()->SetId( 0, ien_array[6*ii] );
    tr->GetPointIds()->SetId( 1, ien_array[6*ii+1] );
    tr->GetPointIds()->SetId( 2, ien_array[6*ii+2] );
    tr->GetPointIds()->SetId( 3, ien_array[6*ii+3] );
    tr->GetPointIds()->SetId( 4, ien_array[6*ii+4] );
    tr->GetPointIds()->SetId( 5, ien_array[6*ii+5] );
    cl -> InsertNextCell(tr);
    tr -> Delete();
  }

  grid_w -> SetCells(22, cl);
  cl -> Delete();
}

void VTK_T::write_triangle_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<int> &node_index,
    const std::vector<int> &ele_index_1,
    const std::vector<int> &ele_index_2 )
{
  // Setup the VTK objects
  vtkPolyData * grid_w = vtkPolyData::New();

  // Generate the mesh
  gen_triangle_grid( grid_w, numpts, numcels, pt, ien_array );

  // nodal indices
  add_int_PointData( grid_w, node_index, "GlobalNodeID" );

  // cell indices
  add_int_CellData( grid_w, ele_index_1, "GlobalElementID_1" );

  add_int_CellData( grid_w, ele_index_2, "GlobalElementID_2" );

  // write vtp
  write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}


void VTK_T::write_quadratic_triangle_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<int> &node_index,
    const std::vector<int> &ele_index_1,
    const std::vector<int> &ele_index_2 )
{
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // Generate the mesh
  gen_quadratic_triangle_grid( grid_w, numpts, numcels, pt, ien_array );

  // nodal indices
  add_int_PointData( grid_w, node_index, "GlobalNodeID" );

  // cell indices
  add_int_CellData( grid_w, ele_index_1, "GlobalElementID_1" );

  add_int_CellData( grid_w, ele_index_2, "GlobalElementID_2" );

  // write vtu (by default of write function)
  write_vtkPointSet(filename, grid_w);

  grid_w->Delete();
}


void VTK_T::add_int_PointData( vtkPointSet * const &grid_w,
    const std::vector<int> &ptdata, const std::string &dataname )
{
  SYS_T::print_fatal_if( ptdata.size() != static_cast<unsigned int>( grid_w -> GetNumberOfPoints() ), "Error: add_int_PointData data size does not match with the number of points.\n" );

  vtkIntArray * data = vtkIntArray::New();
  data -> SetNumberOfComponents(1);
  data -> SetName(dataname.c_str());

  for(unsigned int ii=0; ii<ptdata.size(); ++ii)
    data -> InsertComponent(ii, 0, ptdata[ii]);

  grid_w -> GetPointData() -> AddArray( data );
  data -> Delete();
}


void VTK_T::add_double_PointData( vtkPointSet * const &grid_w,
    const std::vector<double> &ptdata, const std::string &dataname )
{
  SYS_T::print_fatal_if( ptdata.size() != static_cast<unsigned int>( grid_w -> GetNumberOfPoints() ), "Error: add_double_PointData data size does not match with the number of points.\n" );

  vtkDoubleArray * data = vtkDoubleArray::New();
  data -> SetNumberOfComponents(1);
  data -> SetName(dataname.c_str());

  for(unsigned int ii=0; ii<ptdata.size(); ++ii)
    data -> InsertComponent(ii, 0, ptdata[ii]);

  grid_w -> GetPointData() -> AddArray( data );
  data -> Delete();
}


void VTK_T::add_Vector3_PointData( vtkPointSet * const &grid_w,
    const std::vector<Vector_3> &ptdata, const std::string &dataname )
{
  SYS_T::print_fatal_if( ptdata.size() != static_cast<unsigned int>( grid_w -> GetNumberOfPoints() ), "Error: add_Vector3_PointData data size does not match with the number of points.\n" );

  vtkDoubleArray * data = vtkDoubleArray::New();
  data -> SetNumberOfComponents(3);
  data -> SetName(dataname.c_str());
  
  for(unsigned int ii=0; ii<ptdata.size(); ++ii)
  {
    data -> InsertComponent(ii, 0, ptdata[ii].x());
    data -> InsertComponent(ii, 1, ptdata[ii].y());
    data -> InsertComponent(ii, 2, ptdata[ii].z());
  }

  grid_w -> GetPointData() -> AddArray( data );
  data -> Delete();
}


void VTK_T::add_int_CellData( vtkPointSet * const &grid_w,
    const std::vector<int> &cldata, const std::string &dataname )
{
  SYS_T::print_fatal_if( cldata.size() != static_cast<unsigned int>( grid_w -> GetNumberOfCells() ), "Error: add_int_CellData data size does not match with the number of cells.\n" );

  vtkIntArray * data = vtkIntArray::New();
  data -> SetNumberOfComponents(1);
  data -> SetName(dataname.c_str());

  for(unsigned int ii=0; ii<cldata.size(); ++ii)
    data -> InsertComponent(ii, 0, cldata[ii]);

  grid_w -> GetCellData() -> AddArray( data );
  data -> Delete();
}


void VTK_T::add_double_CellData( vtkPointSet * const &grid_w,
    const std::vector<double> &cldata, const std::string &dataname )
{
  SYS_T::print_fatal_if( cldata.size() != static_cast<unsigned int>( grid_w -> GetNumberOfCells() ), "Error: add_double_CellData data size does not match with the number of cells.\n" );

  vtkDoubleArray * data = vtkDoubleArray::New();
  data -> SetNumberOfComponents(1);
  data -> SetName(dataname.c_str());

  for(unsigned int ii=0; ii<cldata.size(); ++ii)
    data -> InsertComponent(ii, 0, cldata[ii]);

  grid_w -> GetCellData() -> AddArray( data );
  data -> Delete();
}


void VTK_T::write_vtkPointSet( const std::string &filename,
    vtkPointSet * const &grid_w, const bool &isXML )
{
  if( grid_w -> GetDataObjectType() == VTK_UNSTRUCTURED_GRID )
  {
    if ( isXML )
    {
      vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
      std::string name_to_write(filename);
      name_to_write.append(".vtu");
      writer -> SetFileName( name_to_write.c_str() );

      writer->SetInputData(grid_w);
      writer->Write();
      writer->Delete();
    }
    else
    {
      vtkUnstructuredGridWriter * writer = vtkUnstructuredGridWriter::New();
      std::string name_to_write(filename);
      name_to_write.append(".vtk");
      writer -> SetFileName( name_to_write.c_str() );

      writer->SetInputData(grid_w);
      writer->Write();
      writer->Delete();
    }
  }
  else if( grid_w -> GetDataObjectType() == VTK_POLY_DATA )
  {
    vtkXMLPolyDataWriter * writer = vtkXMLPolyDataWriter::New();
    std::string name_to_write(filename);
    name_to_write.append(".vtp");
    writer -> SetFileName( name_to_write.c_str() );
    writer->SetInputData(grid_w);
    writer->Write();
    writer->Delete();
  }
  else
    SYS_T::print_fatal("Error: VTK_T::write_vtkPointSet unknown vtkPointSet data. \n");
}

// EOF
