#include "Tet_Tools.hpp"

void TET_T::read_vtu_grid( const std::string &filename,
    int &numpts, int &numcels,
    std::vector<double> &pt, std::vector<int> &ien_array )
{
  vtkXMLUnstructuredGridReader * reader = vtkXMLUnstructuredGridReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  vtkUnstructuredGrid * vtkugrid = reader -> GetOutput();

  // Number of point in the mesh
  numpts  = static_cast<int>( vtkugrid -> GetNumberOfPoints() );
  
  // Number of cell in the mesh
  numcels = static_cast<int>( vtkugrid -> GetNumberOfCells() );

  // Coordinates of the points
  double pt_xyz[3];
  pt.clear();
  for(int ii=0; ii<numpts; ++ii)
  {
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

    if( cell->GetCellType() != 10 ) SYS_T::print_fatal("Error: read_vtu_grid read a mesh with non-tet elements. \n"); 

    ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(3) ) );
  }

  // Close the mesh
  reader->Delete();
}


void TET_T::read_vtu_grid( const std::string &filename,
    int &numpts, int &numcels,
    std::vector<double> &pt, std::vector<int> &ien_array,
    std::vector<int> &phy_tag )
{
  vtkXMLUnstructuredGridReader * reader = vtkXMLUnstructuredGridReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  vtkUnstructuredGrid * vtkugrid = reader -> GetOutput();

  numpts  = static_cast<int>( vtkugrid -> GetNumberOfPoints() );
  numcels = static_cast<int>( vtkugrid -> GetNumberOfCells() );
  
  vtkCellData * celldata = vtkugrid->GetCellData();
  vtkDataArray * cd = celldata->GetScalars("Physics_tag");

  double pt_xyz[3];
  pt.clear();
  for(int ii=0; ii<numpts; ++ii)
  {
    vtkugrid -> GetPoint(ii, pt_xyz);
    pt.push_back(pt_xyz[0]);
    pt.push_back(pt_xyz[1]);
    pt.push_back(pt_xyz[2]);
  }
  
  ien_array.clear();
  phy_tag.clear();
  for(int ii=0; ii<numcels; ++ii)
  {
    vtkCell * cell = vtkugrid -> GetCell(ii);

    if( cell->GetCellType() != 10 ) SYS_T::print_fatal("Error: read_vtu_grid read a mesh with non-tet elements. \n"); 

    ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(3) ) );

    phy_tag.push_back( static_cast<int>( cd->GetComponent(ii, 0) ) );
  }

  reader->Delete();
}


void TET_T::read_vtp_grid( const std::string &filename,
    int &numpts, int &numcels,
    std::vector<double> &pt, std::vector<int> &ien_array,
    std::vector<int> &global_node_index,
    std::vector<int> &global_ele_index )
{
  vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();

  vtkPolyData * polydata = reader -> GetOutput();
  numpts = static_cast<int>( polydata -> GetNumberOfPoints() );
  numcels = static_cast<int>( polydata -> GetNumberOfPolys() );

  // read the cell data
  vtkCellData * celldata = polydata->GetCellData();
  vtkDataArray * cd = celldata->GetScalars("ElemIndex");

  // read the point nodal data
  vtkPointData * pointdata = polydata->GetPointData();
  vtkDataArray * pd = pointdata->GetScalars("NodalIndex");

  pt.clear();
  global_node_index.clear();
  double pt_xyz[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    polydata -> GetPoint(ii, pt_xyz);
    pt.push_back(pt_xyz[0]);
    pt.push_back(pt_xyz[1]);
    pt.push_back(pt_xyz[2]);

    global_node_index.push_back( static_cast<int>(pd->GetComponent(ii,0)) );
  }

  ien_array.clear();
  global_ele_index.clear();
  for(int ii=0; ii<numcels; ++ii)
  {
    vtkCell * cell = polydata -> GetCell(ii);

    if( cell->GetCellType() != 5 ) SYS_T::print_fatal("Error: read_vtp_grid read a mesh with non triangle elements. \n");

    ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
  
    global_ele_index.push_back( static_cast<int>( cd->GetComponent(ii, 0)));
  }

  reader->Delete();
}


void TET_T::read_vtp_grid( const std::string &filename,
    int &numpts, int &numcels,
    std::vector<double> &pt, std::vector<int> &ien_array )
{
  vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();

  vtkPolyData * polydata = reader -> GetOutput();
  numpts = static_cast<int>( polydata -> GetNumberOfPoints() );
  numcels = static_cast<int>( polydata -> GetNumberOfPolys() );
 
  pt.clear();
  double pt_xyz[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    polydata -> GetPoint(ii, pt_xyz);
    pt.push_back(pt_xyz[0]);
    pt.push_back(pt_xyz[1]);
    pt.push_back(pt_xyz[2]);
  }

  ien_array.clear();
  for(int ii=0; ii<numcels; ++ii)
  {
    vtkCell * cell = polydata -> GetCell(ii);

    if( cell->GetCellType() != 5 ) SYS_T::print_fatal("Error: read_vtp_grid read a mesh with non triangle elements. \n");

    ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
    ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
  }

  reader->Delete();
}


void TET_T::write_tet_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt, const std::vector<int> &ien_array )
{
  // Check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: TET_T::write_tet_grid point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 4*numcels) SYS_T::print_fatal("Error: TET_T::write_tet_grid ien array size does not match the numver of cells. \n");
 
  std::vector<int> temp = ien_array;
  VEC_T::sort_unique_resize(temp);
  if( int(temp.size()) != numpts ) SYS_T::print_fatal("Error: TET_T::write_tet_grid numpts is more than the point needed for the cells. Please re-organize the input. \n");
  VEC_T::clean(temp); 

  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();
  
  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New(); 
  ppt->SetDataTypeToDouble();
  double coor[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    coor[0] = pt[3*ii];
    coor[1] = pt[3*ii+1];
    coor[2] = pt[3*ii+2];
    ppt -> InsertPoint(ii, coor);
  }
  
  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  // 2. Cell
  vtkDoubleArray * edge_aspect_ratio = vtkDoubleArray::New();
  edge_aspect_ratio -> SetName("Aspect_ratio");
  edge_aspect_ratio -> SetNumberOfComponents(1);

  vtkCell * cl = vtkTetra::New();
  for(int ii=0; ii<numcels; ++ii)
  {
    cl->GetPointIds()->SetId( 0, ien_array[4*ii] );
    cl->GetPointIds()->SetId( 1, ien_array[4*ii+1] );
    cl->GetPointIds()->SetId( 2, ien_array[4*ii+2] );
    cl->GetPointIds()->SetId( 3, ien_array[4*ii+3] );
    
    grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

    std::vector<double> cell_node; cell_node.clear();
    for(int lnode=0; lnode<4; ++lnode)
    {
      int node_offset = 3 * ien_array[4*ii + lnode];
      cell_node.push_back(pt[node_offset]);
      cell_node.push_back(pt[node_offset+1]);
      cell_node.push_back(pt[node_offset+2]);
    }
    edge_aspect_ratio -> InsertNextValue( TET_T::get_aspect_ratio(cell_node) );
  }
  cl -> Delete();
  grid_w -> GetCellData() -> AddArray( edge_aspect_ratio );
  edge_aspect_ratio -> Delete();

  // write vtk
  vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
  std::string name_to_write(filename);
  name_to_write.append(".vtu");
  writer -> SetFileName( name_to_write.c_str() );

  writer->SetInputData(grid_w);
  writer->Write();
  writer->Delete();
  grid_w->Delete();
}


void TET_T::write_tet_grid_node_elem_index( 
    const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt, const std::vector<int> &ien_array )
{
  // Check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: TET_T::write_tet_grid point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 4*numcels) SYS_T::print_fatal("Error: TET_T::write_tet_grid ien array size does not match the numver of cells. \n");
   
  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();
  
  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New(); 
  ppt->SetDataTypeToDouble();
  double coor[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    coor[0] = pt[3*ii];
    coor[1] = pt[3*ii+1];
    coor[2] = pt[3*ii+2];
    ppt -> InsertPoint(ii, coor);
  }
  
  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  vtkIntArray * ptindex = vtkIntArray::New();
  ptindex -> SetNumberOfComponents(1);
  ptindex -> SetName("NodalIndex");
  for(int ii=0; ii<numpts; ++ii)
    ptindex -> InsertComponent(ii, 0, ii );
  
  grid_w -> GetPointData() -> AddArray( ptindex );
  ptindex->Delete();
  
  // 2. Cell
  vtkDoubleArray * edge_aspect_ratio = vtkDoubleArray::New();
  edge_aspect_ratio -> SetName("Aspect_ratio");
  edge_aspect_ratio -> SetNumberOfComponents(1);

  vtkIntArray * cellindex = vtkIntArray::New();
  cellindex -> SetName("ElemIndex");
  cellindex -> SetNumberOfComponents(1);

  vtkCell * cl = vtkTetra::New();
  for(int ii=0; ii<numcels; ++ii)
  {
    cellindex -> InsertNextValue(ii);

    cl->GetPointIds()->SetId( 0, ien_array[4*ii] );
    cl->GetPointIds()->SetId( 1, ien_array[4*ii+1] );
    cl->GetPointIds()->SetId( 2, ien_array[4*ii+2] );
    cl->GetPointIds()->SetId( 3, ien_array[4*ii+3] );
  
    grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

    std::vector<double> cell_node; cell_node.clear();
    for(int lnode=0; lnode<4; ++lnode)
    {
      int node_offset = 3 * ien_array[4*ii + lnode];
      cell_node.push_back(pt[node_offset]);
      cell_node.push_back(pt[node_offset+1]);
      cell_node.push_back(pt[node_offset+2]);
    }
    edge_aspect_ratio -> InsertNextValue( TET_T::get_aspect_ratio(cell_node) );
  }
  cl -> Delete();
  grid_w -> GetCellData() -> AddArray( edge_aspect_ratio );
  grid_w -> GetCellData() -> AddArray( cellindex );
  edge_aspect_ratio -> Delete();
  cellindex -> Delete();

  // write vtk
  vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
  std::string name_to_write(filename);
  name_to_write.append(".vtu");
  writer -> SetFileName( name_to_write.c_str() );

  writer->SetInputData(grid_w);
  writer->Write();
  writer->Delete();
  grid_w->Delete();
}


void TET_T::write_tet_grid_node_elem_index( 
    const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt, const std::vector<int> &ien_array,
    const std::vector<int> &node_idx, const std::vector<int> &elem_idx )
{
  // Check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: TET_T::write_tet_grid point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 4*numcels) SYS_T::print_fatal("Error: TET_T::write_tet_grid ien array size does not match the numver of cells. \n");
 
  if(int(node_idx.size()) != numpts) SYS_T::print_fatal("Error: TET_T::write_tet_grid node_idx size does not match the number of points. \n");

  if(int(elem_idx.size()) != numcels) SYS_T::print_fatal("Error: TET_T::write_tet_grid elem_idx size does not match the number of cells. \n"); 

  // Setup the VTK objects
  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();
  
  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New(); 
  ppt->SetDataTypeToDouble();
  double coor[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    coor[0] = pt[3*ii];
    coor[1] = pt[3*ii+1];
    coor[2] = pt[3*ii+2];
    ppt -> InsertPoint(ii, coor);
  }
  
  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  vtkIntArray * ptindex = vtkIntArray::New();
  ptindex -> SetNumberOfComponents(1);
  ptindex -> SetName("NodalIndex");
  for(int ii=0; ii<numpts; ++ii)
    ptindex -> InsertComponent(ii, 0, node_idx[ii] );
  
  grid_w -> GetPointData() -> AddArray( ptindex );
  ptindex->Delete();
  
  // 2. Cell
  vtkDoubleArray * edge_aspect_ratio = vtkDoubleArray::New();
  edge_aspect_ratio -> SetName("Aspect_ratio");
  edge_aspect_ratio -> SetNumberOfComponents(1);

  vtkIntArray * cellindex = vtkIntArray::New();
  cellindex -> SetName("ElemIndex");
  cellindex -> SetNumberOfComponents(1);

  vtkCell * cl = vtkTetra::New();
  for(int ii=0; ii<numcels; ++ii)
  {
    cellindex -> InsertNextValue( elem_idx[ii] );

    cl->GetPointIds()->SetId( 0, ien_array[4*ii] );
    cl->GetPointIds()->SetId( 1, ien_array[4*ii+1] );
    cl->GetPointIds()->SetId( 2, ien_array[4*ii+2] );
    cl->GetPointIds()->SetId( 3, ien_array[4*ii+3] );
  
    grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

    std::vector<double> cell_node; cell_node.clear();
    for(int lnode=0; lnode<4; ++lnode)
    {
      int node_offset = 3 * ien_array[4*ii + lnode];
      cell_node.push_back(pt[node_offset]);
      cell_node.push_back(pt[node_offset+1]);
      cell_node.push_back(pt[node_offset+2]);
    }
    edge_aspect_ratio -> InsertNextValue( TET_T::get_aspect_ratio(cell_node) );
  }
  cl -> Delete();
  grid_w -> GetCellData() -> AddArray( edge_aspect_ratio );
  grid_w -> GetCellData() -> AddArray( cellindex );
  edge_aspect_ratio -> Delete();
  cellindex -> Delete();

  // write vtk
  vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
  std::string name_to_write(filename);
  name_to_write.append(".vtu");
  writer -> SetFileName( name_to_write.c_str() );

  writer->SetInputData(grid_w);
  writer->Write();
  writer->Delete();
  grid_w->Delete();
}


void TET_T::write_tet_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt, const std::vector<int> &ien_array,
    const std::vector<int> &phytag, const bool &isXML,
    const int &start_cell_index )
{
  // Check the compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: TET_T::write_tet_grid: point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 4*numcels) SYS_T::print_fatal("Error: TET_T::write_tet_grid: ien array size does not match the numver of cells. \n");

  if(int(phytag.size()) != numcels) SYS_T::print_fatal("Error: TET_T::write_tet_grid: phytag length does not match the number of cells. \n");

  std::vector<int> temp = ien_array;
  VEC_T::sort_unique_resize(temp);
  if( int(temp.size()) != numpts ) SYS_T::print_fatal("Error: TET_T::write_tet_grid numpts is more than the point needed for the cells. Please re-organize the input. \n");
  
  if( temp.back()+1 != numpts ) SYS_T::print_fatal("Error: TET_T::write_tet_grid IEN array is not ranged from zero to numpts. \n");
  VEC_T::clean(temp); 

  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New();
  ppt->SetDataTypeToDouble();
  double coor[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    coor[0] = pt[3*ii];
    coor[1] = pt[3*ii+1];
    coor[2] = pt[3*ii+2];
    ppt -> InsertPoint(ii, coor);
  }
  grid_w -> SetPoints(ppt);
  ppt -> Delete();
  
  vtkIntArray * ptindex = vtkIntArray::New();
  ptindex -> SetNumberOfComponents(1);
  ptindex -> SetName("NodalIndex");
  for(int ii=0; ii<numpts; ++ii)
    ptindex -> InsertComponent(ii, 0, ii);

  grid_w -> GetPointData() -> AddArray( ptindex );
  ptindex->Delete();

  // 2. Cell
  vtkDoubleArray * edge_aspect_ratio = vtkDoubleArray::New();
  edge_aspect_ratio -> SetName("Aspect_ratio");
  edge_aspect_ratio -> SetNumberOfComponents(1);

  vtkIntArray * phy_tag = vtkIntArray::New();
  phy_tag -> SetName("Physics_tag");
  phy_tag -> SetNumberOfComponents(1);

  vtkIntArray * cellindex = vtkIntArray::New();
  cellindex -> SetName("ElemIndex");
  cellindex -> SetNumberOfComponents(1);

  vtkCell * cl = vtkTetra::New();
  for(int ii=0; ii<numcels; ++ii)
  {
    // cell geometry
    cl->GetPointIds()->SetId( 0, ien_array[4*ii] );
    cl->GetPointIds()->SetId( 1, ien_array[4*ii+1] );
    cl->GetPointIds()->SetId( 2, ien_array[4*ii+2] );
    cl->GetPointIds()->SetId( 3, ien_array[4*ii+3] );
    grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

    // Obtain the cell node coordinates and calculate teh aspect ratio
    std::vector<double> cell_node; cell_node.clear();
    for(int lnode=0; lnode<4; ++lnode)
    {
      int node_offset = 3 * ien_array[4*ii + lnode];
      cell_node.push_back(pt[node_offset]);
      cell_node.push_back(pt[node_offset+1]);
      cell_node.push_back(pt[node_offset+2]);
    }
    edge_aspect_ratio -> InsertNextValue( TET_T::get_aspect_ratio(cell_node) );

    cellindex -> InsertNextValue( ii + start_cell_index );

    phy_tag -> InsertNextValue( phytag[ii] );
  }
  cl -> Delete();
  
  grid_w -> GetCellData() -> AddArray( edge_aspect_ratio );
  grid_w -> GetCellData() -> AddArray( cellindex );
  grid_w -> GetCellData() -> AddArray( phy_tag );
  
  edge_aspect_ratio -> Delete();
  cellindex -> Delete();
  phy_tag -> Delete();

  if( isXML )
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
  grid_w->Delete();
}


void TET_T::write_tet_grid( const std::string &filename,
    const bool &isXML, const int &numpts, const int &numcels,
    const std::vector<double> &pt, const std::vector<int> &ien_array,
    const std::vector<int> &phytag, 
    const std::vector<int> &node_index,
    const std::vector<int> &cell_index )
{
  // Check the compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: TET_T::write_tet_grid: point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 4*numcels) SYS_T::print_fatal("Error: TET_T::write_tet_grid: ien array size does not match the numver of cells. \n");

  if(int(phytag.size()) != numcels) SYS_T::print_fatal("Error: TET_T::write_tet_grid: phytag length does not match the number of cells. \n");

  if(int(node_index.size()) != numpts) SYS_T::print_fatal("Error: TET_T::write_tet_grid: node index array does not match the number of points. \n");

  if(int(cell_index.size()) != numcels) SYS_T::print_fatal("Error: TET_T::write_tet_grid: cell index array does not match the number of cells. \n");

  std::vector<int> temp = ien_array;
  VEC_T::sort_unique_resize(temp);
  if( int(temp.size()) != numpts ) SYS_T::print_fatal("Error: TET_T::write_tet_grid numpts is more than the point needed for the cells. Please re-organize the input. \n");

  if( temp.back()+1 != numpts ) SYS_T::print_fatal("Error: TET_T::write_tet_grid IEN array is not ranged from zero to numpts. \n");
  VEC_T::clean(temp); 

  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New();
  ppt->SetDataTypeToDouble();
  double coor[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    coor[0] = pt[3*ii];
    coor[1] = pt[3*ii+1];
    coor[2] = pt[3*ii+2];
    ppt -> InsertPoint(ii, coor);
  }
  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  vtkIntArray * ptindex = vtkIntArray::New();
  ptindex -> SetNumberOfComponents(1);
  ptindex -> SetName("NodalIndex");
  for(int ii=0; ii<numpts; ++ii)
    ptindex -> InsertComponent( ii, 0, node_index[ii] );

  grid_w -> GetPointData() -> AddArray( ptindex );
  ptindex->Delete();

  // 2. Cell
  vtkDoubleArray * edge_aspect_ratio = vtkDoubleArray::New();
  edge_aspect_ratio -> SetName("Aspect_ratio");
  edge_aspect_ratio -> SetNumberOfComponents(1);

  vtkIntArray * phy_tag = vtkIntArray::New();
  phy_tag -> SetName("Physics_tag");
  phy_tag -> SetNumberOfComponents(1);

  vtkIntArray * cellindex = vtkIntArray::New();
  cellindex -> SetName("ElemIndex");
  cellindex -> SetNumberOfComponents(1);

  vtkCell * cl = vtkTetra::New();
  for(int ii=0; ii<numcels; ++ii)
  {
    // cell geometry
    cl->GetPointIds()->SetId( 0, ien_array[4*ii] );
    cl->GetPointIds()->SetId( 1, ien_array[4*ii+1] );
    cl->GetPointIds()->SetId( 2, ien_array[4*ii+2] );
    cl->GetPointIds()->SetId( 3, ien_array[4*ii+3] );
    grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

    // Obtain the cell node coordinates and calculate teh aspect ratio
    std::vector<double> cell_node; cell_node.clear();
    for(int lnode=0; lnode<4; ++lnode)
    {
      int node_offset = 3 * ien_array[4*ii + lnode];
      cell_node.push_back(pt[node_offset]);
      cell_node.push_back(pt[node_offset+1]);
      cell_node.push_back(pt[node_offset+2]);
    }
    edge_aspect_ratio -> InsertNextValue( TET_T::get_aspect_ratio(cell_node) );

    cellindex -> InsertNextValue( cell_index[ii] );

    phy_tag -> InsertNextValue( phytag[ii] );
  }
  cl -> Delete();

  grid_w -> GetCellData() -> AddArray( edge_aspect_ratio );
  grid_w -> GetCellData() -> AddArray( cellindex );
  grid_w -> GetCellData() -> AddArray( phy_tag );

  edge_aspect_ratio -> Delete();
  cellindex -> Delete();
  phy_tag -> Delete();

  if( isXML )
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
  grid_w->Delete();
}


void TET_T::write_triangle_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<int> &node_index,
    const std::vector<int> &ele_index )
{
  // check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: TET_T::write_triangle_grid point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 3*numcels) SYS_T::print_fatal("Error: TET_T::write_triangle_grid ien array size does not match the number of cells. \n");

  if(int(node_index.size()) != numpts) SYS_T::print_fatal("Error: TET_T::write_triangle_grid node_index size does not match the number of points. \n"); 

  if(int(ele_index.size()) != numcels) SYS_T::print_fatal("Error: TET_T::write_triangle_grid ele_index size does not match the number of cells. \n");

  // Setup the VTK objects
  vtkPolyData * grid_w = vtkPolyData::New();

  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New(); 
  ppt->SetDataTypeToDouble();
  double coor[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    coor[0] = pt[3*ii];
    coor[1] = pt[3*ii+1];
    coor[2] = pt[3*ii+2];
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

  // 3. nodal indices
  vtkIntArray * ptindex = vtkIntArray::New();
  ptindex -> SetNumberOfComponents(1);
  ptindex -> SetName("NodalIndex");
  for(int ii=0; ii<numpts; ++ii)
  {
    ptindex -> InsertComponent(ii, 0, node_index[ii]);
  }
  grid_w -> GetPointData() -> AddArray( ptindex );
  ptindex->Delete();

  // 4. cell indices
  vtkIntArray * clindex = vtkIntArray::New();
  clindex -> SetName("ElemIndex");
  clindex -> SetNumberOfComponents(1);
  for(int ii=0; ii<numcels; ++ii)
    clindex -> InsertNextValue( ele_index[ii] );
  grid_w -> GetCellData() -> AddArray( clindex );
  clindex -> Delete();

  // write vtk
  vtkXMLPolyDataWriter * writer = vtkXMLPolyDataWriter::New();
  std::string name_to_write(filename);
  name_to_write.append(".vtp");
  writer -> SetFileName( name_to_write.c_str() );
  writer->SetInputData(grid_w);
  writer->Write();

  writer->Delete();
  grid_w->Delete();
}


void TET_T::write_triangle_grid( const std::string &filename,
    const int &numpts, const int &numcels,
    const std::vector<double> &pt,
    const std::vector<int> &ien_array,
    const std::vector<int> &node_index,
    const std::vector<int> &ele_index_1,
    const std::vector<int> &ele_index_2 )
{
  // check the input data compatibility
  if(int(pt.size()) != 3*numpts) SYS_T::print_fatal("Error: TET_T::write_triangle_grid point vector size does not match the number of points. \n");

  if(int(ien_array.size()) != 3*numcels) SYS_T::print_fatal("Error: TET_T::write_triangle_grid ien array size does not match the number of cells. \n");

  if(int(node_index.size()) != numpts) SYS_T::print_fatal("Error: TET_T::write_triangle_grid node_index size does not match the number of points. \n"); 

  if(int(ele_index_1.size()) != numcels) SYS_T::print_fatal("Error: TET_T::write_triangle_grid ele_index_1 size does not match the number of cells. \n");

  if(int(ele_index_2.size()) != numcels) SYS_T::print_fatal("Error: TET_T::write_triangle_grid ele_index_2 size does not match the number of cells. \n");

  // Setup the VTK objects
  vtkPolyData * grid_w = vtkPolyData::New();

  // 1. nodal points
  vtkPoints * ppt = vtkPoints::New(); 
  ppt->SetDataTypeToDouble();
  double coor[3];
  for(int ii=0; ii<numpts; ++ii)
  {
    coor[0] = pt[3*ii];
    coor[1] = pt[3*ii+1];
    coor[2] = pt[3*ii+2];
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

  // 3. nodal indices
  vtkIntArray * ptindex = vtkIntArray::New();
  ptindex -> SetNumberOfComponents(1);
  ptindex -> SetName("NodalIndex");
  for(int ii=0; ii<numpts; ++ii)
  {
    ptindex -> InsertComponent(ii, 0, node_index[ii]);
  }
  grid_w -> GetPointData() -> AddArray( ptindex );
  ptindex->Delete();

  // 4. cell indices
  vtkIntArray * clindex = vtkIntArray::New();
  clindex -> SetName("ElemIndex_1");
  clindex -> SetNumberOfComponents(1);
  for(int ii=0; ii<numcels; ++ii)
    clindex -> InsertNextValue( ele_index_1[ii] );
  grid_w -> GetCellData() -> AddArray( clindex );
  clindex -> Delete();

  vtkIntArray * clindex2 = vtkIntArray::New();
  clindex2 -> SetName("ElemIndex_2");
  clindex2 -> SetNumberOfComponents(1);
  for(int ii=0; ii<numcels; ++ii)
    clindex2 -> InsertNextValue( ele_index_2[ii] );
  grid_w -> GetCellData() -> AddArray( clindex2 );
  clindex2 -> Delete();

  // write vtk
  vtkXMLPolyDataWriter * writer = vtkXMLPolyDataWriter::New();
  std::string name_to_write(filename);
  name_to_write.append(".vtp");
  writer -> SetFileName( name_to_write.c_str() );

  writer->SetInputData(grid_w);
  writer->Write();
  writer->Delete();
  grid_w->Delete();
}


double TET_T::get_aspect_ratio( const std::vector<double> &pt )
{
  double edge[6];
  // e01
  edge[0] = MATH_T::norm2(pt[3]-pt[0], pt[4]-pt[1], pt[5]-pt[2]) ;
  // e12
  edge[1] = MATH_T::norm2(pt[6]-pt[3], pt[7]-pt[4], pt[8]-pt[5]) ;
  // e02
  edge[2] = MATH_T::norm2(pt[6]-pt[0], pt[7]-pt[1], pt[8]-pt[2]) ;
  // e03
  edge[3] = MATH_T::norm2(pt[9]-pt[0], pt[10]-pt[1], pt[11]-pt[2]) ;
  // e13
  edge[4] = MATH_T::norm2(pt[9]-pt[3], pt[10]-pt[4], pt[11]-pt[5]) ;
  // e23
  edge[5] = MATH_T::norm2(pt[9]-pt[6], pt[10]-pt[7], pt[11]-pt[8]) ;

  const double emax = *std::max_element(edge, edge+6);
  const double emin = *std::min_element(edge, edge+6);

  return emax / emin;
}


void TET_T::get_out_normal( const std::string &vtpfile,
    const std::vector<double> &vol_ctrlPts,
    const IIEN * const &vol_ien,
    std::vector<double> &outVec )
{
  int numpts, numcels;
  std::vector<double> pts;
  std::vector<int> ien, gnode, gelem;

  TET_T::read_vtp_grid( vtpfile, numpts, numcels, pts, ien, gnode, gelem );

  std::vector<int> trn; trn.resize(3); 
  trn[0] = gnode[ ien[0] ]; // triangle nodes' global indices
  trn[1] = gnode[ ien[1] ];
  trn[2] = gnode[ ien[2] ];

  int tete0 = gelem[0]; // triangle's associated tet element indices

  SYS_T::print_fatal_if(tete0 == -1, "Error: TET_T::get_out_normal requires the element indices for the vtp file.\n");

  std::vector<int> ten; ten.resize(4);
  ten[0] = vol_ien->get_IEN(tete0, 0);
  ten[1] = vol_ien->get_IEN(tete0, 1);
  ten[2] = vol_ien->get_IEN(tete0, 2);
  ten[3] = vol_ien->get_IEN(tete0, 3);

  bool gotnode[4];

  int inside_node = 0, node_check = 0;
  for(int ii=0; ii<4; ++ii)
  {
    gotnode[ii] = VEC_T::is_invec( trn, ten[ii] );

    if(!gotnode[ii]) inside_node = ten[ii];
    else node_check += 1;
  }

  SYS_T::print_fatal_if(node_check!=3, "Error: TET_T::get_out_normal, the associated tet element is incompatible with the triangle element. \n");

  // make cross line-0-1 and line-0-2
  const double l01x = vol_ctrlPts[3*trn[1]] - vol_ctrlPts[3*trn[0]];
  const double l01y = vol_ctrlPts[3*trn[1]+1] - vol_ctrlPts[3*trn[0]+1];
  const double l01z = vol_ctrlPts[3*trn[1]+2] - vol_ctrlPts[3*trn[0]+2];

  const double l02x = vol_ctrlPts[3*trn[2]] - vol_ctrlPts[3*trn[0]];
  const double l02y = vol_ctrlPts[3*trn[2]+1] - vol_ctrlPts[3*trn[0]+1];
  const double l02z = vol_ctrlPts[3*trn[2]+2] - vol_ctrlPts[3*trn[0]+2];

  double oux, ouy, ouz;

  MATH_T::cross3d( l01x, l01y, l01z, l02x, l02y, l02z, oux, ouy, ouz );

  MATH_T::normalize3d( oux, ouy, ouz );

  // obtain the line from tri node 0 to the out-of-surface node
  const double inwx = vol_ctrlPts[inside_node*3] - vol_ctrlPts[3*trn[0]];
  const double inwy = vol_ctrlPts[inside_node*3+1] - vol_ctrlPts[3*trn[0]+1];
  const double inwz = vol_ctrlPts[inside_node*3+2] - vol_ctrlPts[3*trn[0]+2];

  // inner product outward with inward
  const double out_dot_in = MATH_T::dot3d(oux, ouy, ouz, inwx, inwy, inwz);

  if(out_dot_in > 0)
  {
    oux *= -1.0;
    ouy *= -1.0;
    ouz *= -1.0;
  }

  outVec.clear(); outVec.resize(3);
  outVec[0] = oux; outVec[1] = ouy; outVec[2] = ouz;
}


void TET_T::tetgenio2vtu( const tetgenio &meshout, const std::string &fName )
{
  const int index_offset = meshout.firstnumber;
  const int numpts = meshout.numberofpoints;
  const int numcel = meshout.numberoftetrahedra;

  std::vector<double> pt;
  pt.resize( numpts * 3 );
  std::vector<int> ien;
  ien.resize( numcel * 4 );

  std::cout<<"Preparing for writing .vtu grid file ... \n";
  std::cout<<"  -- "<<numpts<<" mesh points detected. \n";
  std::cout<<"  -- "<<numcel<<" mesh tetrahedral cell detected. \n";

  for(int ii=0; ii<numpts*3; ++ii) pt[ii] = meshout.pointlist[ii];

  for(int ii=0; ii<numcel*4; ++ii) ien[ii] = meshout.tetrahedronlist[ii] - index_offset;

  write_tet_grid(fName, numpts, numcel, pt, ien);

  std::cout<<"Volumetric grid file "<<fName<<".vtu has been written on disk. \n";
}


void TET_T::tetgenio2vtu_windex( const tetgenio &meshout, const std::string &fName )
{
  const int index_offset = meshout.firstnumber;
  const int numpts = meshout.numberofpoints;
  const int numcel = meshout.numberoftetrahedra;

  std::vector<double> pt;
  pt.resize( numpts * 3 );
  std::vector<int> ien;
  ien.resize( numcel * 4 );

  std::cout<<"Preparing for writing .vtu grid file ... \n";
  std::cout<<"  -- "<<numpts<<" mesh points detected. \n";
  std::cout<<"  -- "<<numcel<<" mesh tetrahedral cell detected. \n";

  for(int ii=0; ii<numpts*3; ++ii) pt[ii] = meshout.pointlist[ii];

  for(int ii=0; ii<numcel*4; ++ii) ien[ii] = meshout.tetrahedronlist[ii] - index_offset;

  write_tet_grid_node_elem_index(fName, numpts, numcel, pt, ien);

  std::cout<<"Volumetric grid file "<<fName<<".vtu has been written on disk. \n";
}


void TET_T::tetgenio2vtp( const tetgenio &meshout, const std::string &fName,
    const int &bcmarker )
{
  const int index_offset = meshout.firstnumber;
  const int numpts = meshout.numberofpoints;
  const int numcel = meshout.numberoftetrahedra;

  // global point coordinates
  std::vector<double> pt;
  pt.resize( numpts * 3 );
  for(int ii=0; ii<numpts*3; ++ii) pt[ii] = meshout.pointlist[ii];

  // Get global volumetric IEN array
  std::vector<int> ien;
  ien.resize( numcel * 4 );
  for(int ii=0; ii<numcel*4; ++ii) ien[ii] = meshout.tetrahedronlist[ii] - index_offset;

  // Store the point on the boundary in bcpt
  std::cout<<"Preparing for writing boundary "<<bcmarker<<" mesh file \n";
  std::vector<int> bcpt;
  const int numface = meshout.numberoftrifaces;
  bcpt.clear();
  std::cout<<"  -- "<<numface<<" faces detected. \n";
  int bcnumcl = 0;
  std::vector<int> trien_global;
  trien_global.clear();
  for(int ee=0; ee<numface; ++ee)
  {
    if( meshout.trifacemarkerlist[ee] == bcmarker )
    {
      bcnumcl += 1;
      bcpt.push_back( meshout.trifacelist[ee*3] - index_offset );
      bcpt.push_back( meshout.trifacelist[ee*3+1] - index_offset );
      bcpt.push_back( meshout.trifacelist[ee*3+2] - index_offset );
    }
  }
  trien_global = bcpt;

  std::cout<<"  -- "<<bcnumcl<<" faces on BC "<<bcmarker<<" detected. \n";

  // list the bc node in unique ascending order
  VEC_T::sort_unique_resize(bcpt);
  const int bcnumpt = static_cast<int>( bcpt.size() );

  std::cout<<"  -- "<<bcnumpt<<" points on BC "<<bcmarker<<" detected. \n";

  // tript stores the coordinates of the boundary points 
  std::vector<double> tript; tript.clear(); tript.resize(3*bcpt.size());
  for( int ii=0; ii<bcnumpt; ++ii )
  {
    tript[ii*3]   = pt[bcpt[ii]*3] ;
    tript[ii*3+1] = pt[bcpt[ii]*3+1] ;
    tript[ii*3+2] = pt[bcpt[ii]*3+2] ;
  }

  // generate a mapper that maps the bc node to 1, other node to 0
  bool * bcmap = new bool [numpts];
  for(int ii=0; ii<numpts; ++ii) bcmap[ii] = 0;
  for(int ii=0; ii<bcnumpt; ++ii) bcmap[bcpt[ii]] = 1;

  std::vector<int> gelem; gelem.clear();
  for( int ee=0; ee<numcel; ++ee )
  {
    int total = 0;
    total += bcmap[ ien[4*ee] ];
    total += bcmap[ ien[4*ee+1] ];
    total += bcmap[ ien[4*ee+2] ];
    total += bcmap[ ien[4*ee+3] ];

    // get a list of volume element that has at least 3 poitns on
    // the boundary surface 
    if(total >= 3) gelem.push_back(ee);
  }
  delete [] bcmap;
  std::cout<<"  -- "<<gelem.size()<<" tets have face on this BC surface. \n";

  // generate the local triangle IEN array
  std::vector<int> trien; trien.clear(); 
  int node0, node1, node2;
  std::vector<int>::iterator it;
  for( int ee=0; ee<bcnumcl; ++ee)
  {
    node0 = trien_global[3*ee];
    node1 = trien_global[3*ee+1];
    node2 = trien_global[3*ee+2];

    // trien_global and bcpt stores the same items, just ordering is
    // changed. So, it should always find the node in bcpt list.
    it = find(bcpt.begin(), bcpt.end(), node0);
    trien.push_back( it - bcpt.begin() );

    it = find(bcpt.begin(), bcpt.end(), node1);
    trien.push_back( it - bcpt.begin() );

    it = find(bcpt.begin(), bcpt.end(), node2);
    trien.push_back( it - bcpt.begin() );
  }
  std::cout<<"  -- triangle IEN generated. \n"; 

  // determine the face-to-element mapping
  std::vector<int> face2elem; face2elem.resize(bcnumcl);
  int vol_elem; 
  int vnode[4];
  bool got0, got1, got2, gotit;
  for(int ff=0; ff<bcnumcl; ++ff)
  {
    node0 = trien_global[3*ff];
    node1 = trien_global[3*ff+1];
    node2 = trien_global[3*ff+2];
    gotit = false;
    int ee = -1;
    while( !gotit )
    {
      ee += 1;
      vol_elem = gelem[ee];
      vnode[0] = ien[4*vol_elem];
      vnode[1] = ien[4*vol_elem+1];
      vnode[2] = ien[4*vol_elem+2];
      vnode[3] = ien[4*vol_elem+3];
      std::sort(vnode, vnode+4);
      got0 = ( std::find(vnode, vnode+4, node0) != vnode+4 );
      got1 = ( std::find(vnode, vnode+4, node1) != vnode+4 );
      got2 = ( std::find(vnode, vnode+4, node2) != vnode+4 );

      gotit = got0 && got1 && got2;
    }
    face2elem[ff] = gelem[ee];
  }

  // preparing for mesh file (.vtp) name
  std::string outname(fName);
  outname.append(".");
  std::string strbcindex = SYS_T::to_string(bcmarker);
  outname.append(strbcindex);

  write_triangle_grid( outname, bcnumpt, bcnumcl, tript, 
      trien, bcpt, face2elem );
  std::cout<<"Surface file "<<outname<<".vtp has been written on disk. \n";  
}


namespace TET_T
{
  Tet4::Tet4()
  {
    pts[0] = 0.0; pts[1] = 0.0;  pts[2] = 0.0;
    pts[3] = 1.0; pts[4] = 0.0;  pts[5] = 0.0;
    pts[6] = 0.0; pts[7] = 1.0;  pts[8] = 0.0;
    pts[9] = 0.0; pts[10] = 0.0; pts[11] = 1.0;

    gindex[0] = 0; gindex[1] = 1;
    gindex[2] = 2; gindex[3] = 3;
  }


  Tet4::Tet4( const std::vector<double> &in_nodes )
  {
    SYS_T::print_exit_if( in_nodes.size() != 12, "Error: input nodal list shall have 4 nodes with xyz-coordinates. \n");

    for(int ii=0; ii<12; ++ii) pts[ii] = in_nodes[ii];

    gindex[0] = 0; gindex[1] = 1;
    gindex[2] = 2; gindex[3] = 3;
  }


  Tet4::Tet4( const std::vector<double> &ctrlPts,
      const int &ien0, const int &ien1,
      const int &ien2, const int &ien3 )
  {
    gindex[0] = ien0;
    gindex[1] = ien1;
    gindex[2] = ien2;
    gindex[3] = ien3;

    for(int ii=0; ii<4; ++ii)
    {
      pts[ii*3]   = ctrlPts[gindex[ii]*3];
      pts[ii*3+1] = ctrlPts[gindex[ii]*3+1];
      pts[ii*3+2] = ctrlPts[gindex[ii]*3+2];
    }
  }


  Tet4::~Tet4()
  {
  }


  void Tet4::reset( const std::vector<double> &in_nodes )
  {
    SYS_T::print_exit_if( in_nodes.size() != 12, "Error: input nodal list shall have 4 nodes with xyz-coordinates. \n");
    for(int ii=0; ii<12; ++ii) pts[ii] = in_nodes[ii];

    gindex[0] = 0; gindex[1] = 1;
    gindex[2] = 2; gindex[3] = 3;
  }


  void Tet4::reset( const std::vector<double> &ctrlPts,
      const int &ien0, const int &ien1,
      const int &ien2, const int &ien3 )
  {
    gindex[0] = ien0;
    gindex[1] = ien1;
    gindex[2] = ien2;
    gindex[3] = ien3;

    for(int ii=0; ii<4; ++ii)
    {
      pts[ii*3]   = ctrlPts[gindex[ii]*3];
      pts[ii*3+1] = ctrlPts[gindex[ii]*3+1];
      pts[ii*3+2] = ctrlPts[gindex[ii]*3+2];
    }
  }


  void Tet4::reset( const int &ien0, const int &ien1,
      const int &ien2, const int &ien3 )
  {
    gindex[0] = ien0;
    gindex[1] = ien1;
    gindex[2] = ien2;
    gindex[3] = ien3;
  }


  void Tet4::reset( const std::vector<double> &ctrlPts,
      const IIEN * const &ien_ptr, const int &ee )
  {
    gindex[0] = ien_ptr->get_IEN(ee, 0);
    gindex[1] = ien_ptr->get_IEN(ee, 1);
    gindex[2] = ien_ptr->get_IEN(ee, 2);
    gindex[3] = ien_ptr->get_IEN(ee, 3);

    for(int ii=0; ii<4; ++ii)
    {
      pts[ii*3]   = ctrlPts[gindex[ii]*3];
      pts[ii*3+1] = ctrlPts[gindex[ii]*3+1];
      pts[ii*3+2] = ctrlPts[gindex[ii]*3+2];
    }
  }


  void Tet4::print_info() const
  {
    std::cout<<"Tet4 object : \n";
    std::cout<<" -- global indices: "<<gindex[0]<<'\t';
    std::cout<<gindex[1]<<'\t'<<gindex[2]<<'\t'<<gindex[3]<<'\n';
    std::cout<<" -- pt0 : "<<pts[0]<<'\t'<<pts[1]<<'\t'<<pts[2]<<'\n';
    std::cout<<" -- pt1 : "<<pts[3]<<'\t'<<pts[4]<<'\t'<<pts[5]<<'\n';
    std::cout<<" -- pt2 : "<<pts[6]<<'\t'<<pts[7]<<'\t'<<pts[8]<<'\n';
    std::cout<<" -- pt3 : "<<pts[9]<<'\t'<<pts[10]<<'\t'<<pts[11]<<'\n';
  }


  double Tet4::get_aspect_ratio() const
  {
    double edge[6];
    // e01
    edge[0] = MATH_T::norm2(pts[3]-pts[0], pts[4]-pts[1], pts[5]-pts[2]) ;
    // e12
    edge[1] = MATH_T::norm2(pts[6]-pts[3], pts[7]-pts[4], pts[8]-pts[5]) ;
    // e02
    edge[2] = MATH_T::norm2(pts[6]-pts[0], pts[7]-pts[1], pts[8]-pts[2]) ;
    // e03
    edge[3] = MATH_T::norm2(pts[9]-pts[0], pts[10]-pts[1], pts[11]-pts[2]) ;
    // e13
    edge[4] = MATH_T::norm2(pts[9]-pts[3], pts[10]-pts[4], pts[11]-pts[5]) ;
    // e23
    edge[5] = MATH_T::norm2(pts[9]-pts[6], pts[10]-pts[7], pts[11]-pts[8]) ;

    const double emax = *std::max_element(edge, edge+6);
    const double emin = *std::min_element(edge, edge+6);

    return emax / emin;
  }


  int Tet4::get_face_id( const int &n0, const int &n1, const int &n2) const
  {
    int temp[] = { n0, n1, n2 };

    int * it0, * it1, * it2, * it3;

    it0 = std::find(temp, temp+3, gindex[0]);
    it1 = std::find(temp, temp+3, gindex[1]);
    it2 = std::find(temp, temp+3, gindex[2]);
    it3 = std::find(temp, temp+3, gindex[3]);

    // flgx is true if the node is found in temp
    bool flg[4];
    flg[0] = (it0 != temp+3);
    flg[1] = (it1 != temp+3);
    flg[2] = (it2 != temp+3);
    flg[3] = (it3 != temp+3);

    int sum = flg[0] + flg[1] + flg[2] + flg[3];

    SYS_T::print_exit_if( sum != 3, "Error: Tet4::find_face_id input is not a proper face of this element. \n");

    int loc = 0;
    while( flg[loc] ) loc+=1;

    return loc;
  }


  double Tet4::get_diameter() const
  {
    double x, y, z, r;
    MATH_T::get_tet_sphere_info(
        pts[0], pts[3], pts[6], pts[9],
        pts[1], pts[4], pts[7], pts[10],
        pts[2], pts[5], pts[8], pts[11],
        x, y, z, r );

    return 2.0 * r;
  }


  double Tet4::get_volume() const
  {
    const double x10 = pts[3] - pts[0];
    const double y01 = pts[1] - pts[4];
    const double z01 = pts[2] - pts[5];

    const double x21 = pts[6] - pts[3];
    const double y12 = pts[4] - pts[7];
    const double z12 = pts[5] - pts[8];

    const double x32 = pts[9] - pts[6];
    const double z23 = pts[8] - pts[11];
    const double y23 = pts[7] - pts[10];

    const double sum = x10 * (y12*z23 - y23*z12)
      + x21*(y23*z01 - y01*z23) + x32*(y01*z12 - y12*z01);

    return sum / 6.0;
  }


  void Tet4::write_vtu( const std::string &fileName ) const
  {
    // Setup the VTK objects
    vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

    // 1. four nodal points
    vtkPoints * ppt = vtkPoints::New(); 
    double coor[3];
    for(int ii=0; ii<4; ++ii)
    {
      coor[0] = pts[3*ii]; coor[1] = pts[3*ii+1]; coor[2] = pts[3*ii+2];
      ppt -> InsertPoint(ii, coor);
    }

    grid_w -> SetPoints(ppt);
    ppt -> Delete();

    vtkIntArray * ptindex = vtkIntArray::New();
    ptindex -> SetNumberOfComponents(1);
    ptindex -> SetName("NodalIndex");
    for(int ii=0; ii<4; ++ii) ptindex->InsertComponent(ii, 0, gindex[ii]);

    grid_w -> GetPointData() -> AddArray( ptindex );
    ptindex -> Delete();

    // 2. Cell
    vtkDoubleArray * edge_aspect_ratio = vtkDoubleArray::New();
    edge_aspect_ratio -> SetName("Aspect_ratio");
    edge_aspect_ratio -> SetNumberOfComponents(1);

    vtkCell * cl = vtkTetra::New();

    cl->GetPointIds()->SetId( 0, 0 );
    cl->GetPointIds()->SetId( 1, 1 );
    cl->GetPointIds()->SetId( 2, 2 );
    cl->GetPointIds()->SetId( 3, 3 );

    grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );

    cl -> Delete();

    edge_aspect_ratio -> InsertNextValue( get_aspect_ratio() );

    grid_w -> GetCellData() -> AddArray( edge_aspect_ratio );
    edge_aspect_ratio -> Delete();

    // write vtk
    vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
    std::string name_to_write(fileName);
    name_to_write.append(".vtu");
    writer -> SetFileName( name_to_write.c_str() );

    writer->SetInputData(grid_w);
    writer->Write();
    writer->Delete();
    grid_w->Delete();
  }
}


void TET_T::tetmesh_check(const std::vector<double> &cpts,
    const IIEN * const &ienptr, const int &nelem,
    const double &crit_aspect_ratio )
{
  std::cout<<"\n======= Tet4 mesh quality =======\n";
  TET_T::Tet4 * teton = new TET_T::Tet4();
  teton -> reset( cpts, ienptr, 0 );
  double teton_max_vol = teton -> get_volume();
  double teton_min_vol = teton_max_vol;
  double teton_max_h = teton -> get_diameter();
  double teton_min_h = teton_max_h;
  int num_dist_elem = 0; // number of element that has negative volume
  int num_aspt_elem = 0; // number of element that has aspect ratio larger
  // than the given critical aspect ratio value
  for(int ee = 0; ee<nelem; ++ee)
  {
    teton->reset( cpts, ienptr, ee );

    if( teton->get_aspect_ratio() > crit_aspect_ratio )
      num_aspt_elem += 1;

    double teton_ee_vol = teton -> get_volume();
    if( teton_ee_vol < 0.0 ) num_dist_elem += 1;

    if( teton_max_vol < teton_ee_vol) teton_max_vol = teton_ee_vol;

    if( teton_min_vol > teton_ee_vol) teton_min_vol = teton_ee_vol;

    double teton_he = teton -> get_diameter();

    if( teton_max_h < teton_he ) teton_max_h = teton_he;

    if( teton_min_h > teton_he ) teton_min_h = teton_he;
  }
  delete teton;
  cout<<"- maximum tetrahedron volume : "<<teton_max_vol<<endl;
  cout<<"- minimum tetrahedron volume : "<<teton_min_vol<<endl;
  cout<<"- maximum tetrahedron diameter : "<<teton_max_h<<endl;
  cout<<"- minimum tetrahedron diameter : "<<teton_min_h<<endl;
  cout<<"- number of distorted element : "<<num_dist_elem<<endl;
  cout<<"- number of element with aspect ratio larger than "<<crit_aspect_ratio<<" : "<<num_aspt_elem<<endl;
  std::cout<<"==================================\n";
}

// EOF
