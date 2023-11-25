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
    else if( cell->GetCellType() == 12 ) 
    {
      // cell type 12 is eight-node hex
      ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(3) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(4) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(5) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(6) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(7) ) );
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
    else if( cell-> GetCellType() == 28 )
    {
      // cell type 28 is nine-node quad
      ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(3) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(4) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(5) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(6) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(7) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(8) ) );
    }
    else if( cell-> GetCellType() == 29 )
    {
      // cell type 29 is 27-node hex
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
      ien_array.push_back( static_cast<int>( cell->GetPointId(10) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(11) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(12) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(13) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(14) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(15) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(16) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(17) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(18) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(19) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(20) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(21) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(22) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(23) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(24) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(25) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(26) ) );
    }
    else SYS_T::print_fatal("Error: VTK_T::read_vtu_grid read a mesh with VTK cell type %d is not supported.\n", cell-> GetCellType() ); 
  }

  reader->Delete();
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
      // cell type 5 is three-node triangle
      ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
    }
    else if(cell->GetCellType() == 9)
    {
      // cell type 9 is four-node quad
      ien_array.push_back( static_cast<int>( cell->GetPointId(0) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(1) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(2) ) );
      ien_array.push_back( static_cast<int>( cell->GetPointId(3) ) );
    }
    else SYS_T::print_fatal("Error: read_vtp_grid read a mesh with VTK cell type 5 or 9. \n");
  }

  reader->Delete();
}

int VTK_T::read_grid( const std::string &filename,
    int &numpts, int &numcels,
    std::vector<double> &pt, std::vector<int> &ien_array )
{
  int file_type = 0;

  vtkXMLGenericDataObjectReader * reader = vtkXMLGenericDataObjectReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();

  // Obtain the filename extension
  std::string fend;
  fend.assign( filename.end()-4 , filename.end() );

  if(dynamic_cast<vtkPolyData*>(reader->GetOutput()))
  {
    SYS_T::print_fatal_if(fend.compare(".vtp") !=0, "Error: VTK::read_grid, the filename %s does not end with vtp. \n", filename.c_str());
    
    file_type = 1; 

    read_vtp_grid(filename, numpts, numcels, pt, ien_array);
  }
  else if(dynamic_cast<vtkUnstructuredGrid*>(reader->GetOutput()))
  {
    SYS_T::print_fatal_if(fend.compare(".vtu") !=0, "Error: VTK::read_grid, the filename %s does not end with vtu. \n", filename.c_str());
    
    file_type = 2;

    read_vtu_grid(filename, numpts, numcels, pt, ien_array);
  }
  else
    SYS_T::print_fatal("VTK_T::read_grid unknown vtk object type.\n");

  reader -> Delete();

  return file_type;
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

std::vector<Vector_3> VTK_T::read_Vector_3_PointData( const std::string &filename,
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
    SYS_T::print_fatal("VTK_T::read_Vector_3_PointData unknown vtk object type.\n");

  vtkDataArray * pd = pointdata->GetScalars( dataname.c_str() );

  std::vector<Vector_3> data( numpts );
  for(int ii=0; ii<numpts; ++ii)
  {
    const double data_x = static_cast<double>( pd->GetComponent(ii, 0) );
    const double data_y = static_cast<double>( pd->GetComponent(ii, 1) );
    const double data_z = static_cast<double>( pd->GetComponent(ii, 2) );
    data[ii] = Vector_3( data_x, data_y, data_z );
  }

  reader -> Delete();

  return data;
}

int VTK_T::read_num_pt( const std::string &filename )
{
  vtkXMLGenericDataObjectReader * reader = vtkXMLGenericDataObjectReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  
  int numpts = -1;

  // Downcasting will return null if fails
  if(dynamic_cast<vtkPolyData*>(reader->GetOutput()))
  {
    vtkPolyData * vtkgrid = reader -> GetPolyDataOutput (); 
    numpts = static_cast<int>( vtkgrid -> GetNumberOfPoints() );
  }
  else if(dynamic_cast<vtkUnstructuredGrid*>(reader->GetOutput()))
  {
    vtkUnstructuredGrid * vtkgrid = reader -> GetUnstructuredGridOutput();
    numpts = static_cast<int>( vtkgrid -> GetNumberOfPoints() );
  }
  else
    SYS_T::print_fatal("VTK_T::read_num_pt unknown vtk object type.\n");

  reader -> Delete();

  return numpts;
}

int VTK_T::read_num_cl( const std::string &filename )
{
  vtkXMLGenericDataObjectReader * reader = vtkXMLGenericDataObjectReader::New();
  reader -> SetFileName( filename.c_str() );
  reader -> Update();
  
  int numcels = -1;

  // Downcasting will return null if fails
  if(dynamic_cast<vtkPolyData*>(reader->GetOutput()))
  {
    vtkPolyData * vtkgrid = reader -> GetPolyDataOutput (); 
    numcels = static_cast<int>( vtkgrid -> GetNumberOfCells() );
  }
  else if(dynamic_cast<vtkUnstructuredGrid*>(reader->GetOutput()))
  {
    vtkUnstructuredGrid * vtkgrid = reader -> GetUnstructuredGridOutput();
    numcels = static_cast<int>( vtkgrid -> GetNumberOfCells() );
  }
  else
    SYS_T::print_fatal("VTK_T::read_num_cl unknown vtk object type.\n");

  reader -> Delete();

  return numcels;
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
