#include "Vis_Tools.hpp"

void VIS_T::writeVisFile( vtkUnstructuredGrid * gridData,
    const std::string &outputBName, const std::string &pre_pvtuname,
    const int &rank, const int &size, 
    const double &sol_time, const bool &isXML )
{
  // write vtk/vtu files
  SYS_T::commPrint("-- Write vtu files.\n");
  writeGridObject(gridData, pre_pvtuname, rank, isXML);

  // write pvtu file
  writepvtuFile(gridData, pre_pvtuname, rank, size, isXML);
  MPI_Barrier(PETSC_COMM_WORLD);

  // write pvd file
  SYS_T::commPrint("-- Update pvd file.\n");
  writepvdFile(outputBName, pre_pvtuname, rank, size, sol_time, isXML);

  MPI_Barrier(PETSC_COMM_WORLD);
}

void VIS_T::writeGridObject( vtkUnstructuredGrid * gridData,
    const std::string &baseName, const int &rank, const bool &isXML )
{
  std::string fName( baseName );
  fName.append("_p");
  // append processor rank
  if( rank / 10 == 0 )
    fName.append("000");
  else if( rank / 100 == 0 )
    fName.append("00");
  else if( rank/ 1000 == 0)
    fName.append("0");
  std::stringstream sstrm;
  sstrm << rank;
  fName.append(sstrm.str());

  // write vtu or vtk file
  if( isXML )
  {
    fName.append(".vtu");
    vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
    writer -> SetFileName( fName.c_str() );
#if VTK_MAJOR_VERSION <= 5    
    writer -> SetInput( gridData );
#else    
    writer -> SetInputData( gridData );
#endif    
    if( 1 != writer->Write() )
    {
      PetscPrintf(PETSC_COMM_WORLD, "ERROR: Failed to write vtu file %s. \n", fName.c_str());
      writer->Delete();
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    writer->Delete();
  }
  else
  {
    fName.append(".vtk");
    vtkUnstructuredGridWriter * writer = vtkUnstructuredGridWriter::New();
    writer -> SetFileName( fName.c_str() );
#if VTK_MAJOR_VERSION <= 5    
    writer -> SetInput(gridData);
#else    
    writer -> SetInputData(gridData);
#endif
    if( 1 != writer->Write() )
    {
      PetscPrintf(PETSC_COMM_WORLD, "ERROR: Failed to write vtk file %s. \n", fName.c_str());
      writer->Delete();
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
    writer->Delete();
  }
}

void VIS_T::writepvtuFile( vtkUnstructuredGrid * gridData,
    const std::string &baseName, const int &rank, const int &size,
    const bool &isXML)
{
  if( isXML && rank == 0 && size > 1 )
  {
    // setup name
    std::string pName(baseName);
    pName.append(".pvtu");
    std::ofstream pvtuFile;
    pvtuFile.open(pName.c_str());
    if( !pvtuFile.is_open() )
    {
      std::cerr<<"ERROR: failed to open pvtu file "<<pName<<endl;
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    pvtuFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"";
    // Test for big or little endian system
    long one = 1;
    if( (*((char *)(&one))) )
      pvtuFile << "LittleEndian";
    else
      pvtuFile << "BigEndian";
    pvtuFile << "\">" << std::endl;
    pvtuFile << "<PUnstructuredGrid GhostLevel=\"0\">" << std::endl;

    // Cell Data: processor ID associated with the parition
    pvtuFile << "<PCellData>" << std::endl;
    
    vtkCellData * cellData = gridData->GetCellData();
    const int numCellData = cellData->GetNumberOfArrays();
    for( int ii = 0; ii < numCellData; ++ii )
    {
      vtkDataArray * cellArray = cellData -> GetArray( ii );

      if( NULL == cellArray ) continue;

      pvtuFile << "<PDataArray type=\"";
      
      if( VTK_INT == cellArray -> GetDataType() )
        pvtuFile << "Int32";
      else
        SYS_T::print_fatal("Error: unknown VTK_INT type.\n");

      pvtuFile<<"\" Name=\""<< cellArray->GetName()<<"\""
        <<" NumberOfComponents=\""<<cellArray->GetNumberOfComponents()
        <<"\"/>"<<std::endl;
    }
    pvtuFile<<"</PCellData>" << std::endl;

    // Point Data: associated with the gridData object
    vtkPointData* pointData = gridData->GetPointData();
    pvtuFile << "<PPointData>" << std::endl;
    const int numArrays = pointData->GetNumberOfArrays();

    for( int ii = 0; ii < numArrays; ii++ )
    {
      vtkDataArray * currArray = pointData->GetArray( ii );
      if( NULL == currArray )
        continue;
      pvtuFile << "<PDataArray type=\"";
      if( VTK_DOUBLE == currArray->GetDataType()  )
        pvtuFile << "Float64";
      else
        SYS_T::print_fatal("Error: unknown VTK_DOUBLE type.\n");

      pvtuFile << "\" Name=\"" << currArray->GetName() << "\""
        << " NumberOfComponents=\"" << currArray->GetNumberOfComponents() << "\"/>"
        << std::endl;
    }
    pvtuFile << "</PPointData>" << std::endl;

    // Points
    pvtuFile << "<PPoints>" << std::endl
      << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>" << std::endl
      << "</PPoints>" << std::endl;

    for( int ii=0; ii<size; ++ii )
    {
      std::string pieceName(baseName);
      pieceName.append("_p");
      if( ii / 10 == 0 )
        pieceName.append("000");
      else if( ii / 100 == 0 )
        pieceName.append("00");
      else if( ii / 1000 == 0)
        pieceName.append("0");
      std::stringstream sstrm;
      sstrm << ii;
      pieceName.append(sstrm.str());
      pieceName.append(".vtu");
      pvtuFile << "<Piece Source=\"" << pieceName << "\"/>" << std::endl; 
    }

    pvtuFile << "</PUnstructuredGrid>" << std::endl;
    pvtuFile << "</VTKFile>" << std::endl;
    pvtuFile.close();
    std::cout << "-- Results output writen to file: " << pName << std::endl;
  }
}

void VIS_T::writepvdFile_Init( const std::string &pvdFName )
{
  std::ofstream pvdFile;

  pvdFile.open( pvdFName.c_str(), std::ios::out | std::ios::trunc );

  if( !pvdFile.is_open() )
  {
    std::cerr<<"ERROR: UNABLE TO OPEN THE PVD FILE : "<<pvdFName<<endl;
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  // writes the xml top part of the pvd file
  pvdFile << "<?xml version=\"1.0\"?>"  << std::endl
    << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
    << "<Collection>" << std::endl
    << "</Collection>" << std::endl
    << "</VTKFile>" << std::endl; 

  pvdFile.close();
}

void VIS_T::writepvdFile( const std::string &baseName,
    const std::string &pre_pvtuname, const int &rank, const int &size,
    const double &sol_time, const bool &isXML )
{
  if( rank == 0 && isXML == true )
  {
    std::string pvdname( baseName );
    pvdname.append(".pvd");

    // if baseName.pvd does not exists, write the file with top initialized part
    std::fstream fileTest;
    fileTest.open(pvdname.c_str(), std::fstream::in);
    if( !fileTest.is_open() )
      writepvdFile_Init(pvdname);

    fileTest.close();

    // now we are sure the pvd file exists, we just need to append the pvtu and
    // time information to it.
    std::fstream pvdFile;
    pvdFile.open(pvdname.c_str());
    if( !pvdFile.is_open() )
    {
      std::cerr<<"ERROR: FAIL TO OPEN THE PVD FILE "<<pvdname<<endl;
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    std::string pvtuname( pre_pvtuname );
    pvtuname.append(".pvtu");

    // if runs in serial, use vtu name
    std::string serialname( pre_pvtuname );
    serialname.append("_p0000.vtu");

    pvdFile.seekp(-25 , ios::end);
    if( size > 1 )
    {
      pvdFile << "<DataSet timestep=\"" << std::setprecision(12) << sol_time
        << "\" groups=\"\" part=\"0\" "
        << "file=\"" << pvtuname << "\"/>" << std::endl;
    }
    else 
    {
      pvdFile << "<DataSet timestep=\"" << sol_time << "\" groups=\"\" part=\"0\" "
        << "file=\"" << serialname << "\"/>" << std::endl;
    }

    // finish up and close the file with end part
    pvdFile << "</Collection>" << std::endl
      << "</VTKFile>" << std::endl;

    pvdFile.close();
  }
}

void VIS_T::setHexelem( const int &ptoffset, vtkUnstructuredGrid * gridData )
{
  vtkCell * cell = vtkHexahedron::New();

  cell->GetPointIds()->SetId( 0, ptoffset + 0 );
  cell->GetPointIds()->SetId( 1, ptoffset + 1 );
  cell->GetPointIds()->SetId( 2, ptoffset + 2 );
  cell->GetPointIds()->SetId( 3, ptoffset + 3 );
  cell->GetPointIds()->SetId( 4, ptoffset + 4 );
  cell->GetPointIds()->SetId( 5, ptoffset + 5 );
  cell->GetPointIds()->SetId( 6, ptoffset + 6 );
  cell->GetPointIds()->SetId( 7, ptoffset + 7 );

  gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
  cell->Delete();
  
}

void VIS_T::setHexelem( const int &ptid0, const int &ptid1,
      const int &ptid2, const int &ptid3, const int &ptid4,
      const int &ptid5, const int &ptid6, const int &ptid7,
      vtkUnstructuredGrid * gridData )
{
  vtkCell * cell = vtkHexahedron::New();
  
  cell->GetPointIds()->SetId( 0, ptid0 );
  cell->GetPointIds()->SetId( 1, ptid1 );
  cell->GetPointIds()->SetId( 2, ptid2 );
  cell->GetPointIds()->SetId( 3, ptid3 );
  cell->GetPointIds()->SetId( 4, ptid4 );
  cell->GetPointIds()->SetId( 5, ptid5 );
  cell->GetPointIds()->SetId( 6, ptid6 );
  cell->GetPointIds()->SetId( 7, ptid7 );

  gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
  cell->Delete();
}

void VIS_T::setTriQuadHexelem( const int &ptoffset, vtkUnstructuredGrid * gridData )
{
  vtkCell * cell = vtkTriQuadraticHexahedron::New();
  
  cell->GetPointIds()->SetId( 0, ptoffset + 0 );
  cell->GetPointIds()->SetId( 1, ptoffset + 1 );
  cell->GetPointIds()->SetId( 2, ptoffset + 2 );
  cell->GetPointIds()->SetId( 3, ptoffset + 3 );
  cell->GetPointIds()->SetId( 4, ptoffset + 4 );
  cell->GetPointIds()->SetId( 5, ptoffset + 5 );
  cell->GetPointIds()->SetId( 6, ptoffset + 6 );
  cell->GetPointIds()->SetId( 7, ptoffset + 7 );
  cell->GetPointIds()->SetId( 8, ptoffset + 8 );
  cell->GetPointIds()->SetId( 9, ptoffset + 9 );  
  cell->GetPointIds()->SetId( 10, ptoffset + 10 );
  cell->GetPointIds()->SetId( 11, ptoffset + 11 );
  cell->GetPointIds()->SetId( 12, ptoffset + 12 );
  cell->GetPointIds()->SetId( 13, ptoffset + 13 );
  cell->GetPointIds()->SetId( 14, ptoffset + 14 );
  cell->GetPointIds()->SetId( 15, ptoffset + 15 );
  cell->GetPointIds()->SetId( 16, ptoffset + 16 );
  cell->GetPointIds()->SetId( 17, ptoffset + 17 );
  cell->GetPointIds()->SetId( 18, ptoffset + 18 );
  cell->GetPointIds()->SetId( 19, ptoffset + 19 );
  cell->GetPointIds()->SetId( 20, ptoffset + 20 );
  cell->GetPointIds()->SetId( 21, ptoffset + 21 );
  cell->GetPointIds()->SetId( 22, ptoffset + 22 );
  cell->GetPointIds()->SetId( 23, ptoffset + 23 );
  cell->GetPointIds()->SetId( 24, ptoffset + 24 );
  cell->GetPointIds()->SetId( 25, ptoffset + 25 );
  cell->GetPointIds()->SetId( 26, ptoffset + 26 );

  gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
  cell->Delete();
}

void VIS_T::setTriQuadHexelem( const int &ptid0, const int &ptid1,
      const int &ptid2, const int &ptid3, const int &ptid4,
      const int &ptid5, const int &ptid6, const int &ptid7,
      const int &ptid8, const int &ptid9, const int &ptid10,
      const int &ptid11, const int &ptid12, const int &ptid13,
      const int &ptid14, const int &ptid15, const int &ptid16,
      const int &ptid17, const int &ptid18, const int &ptid19,
      const int &ptid20, const int &ptid21, const int &ptid22,
      const int &ptid23, const int &ptid24, const int &ptid25,
      const int &ptid26, vtkUnstructuredGrid * gridData )
{
  vtkCell * cell = vtkTriQuadraticHexahedron::New();
  
  cell->GetPointIds()->SetId( 0, ptid0 );
  cell->GetPointIds()->SetId( 1, ptid1 );
  cell->GetPointIds()->SetId( 2, ptid2 );
  cell->GetPointIds()->SetId( 3, ptid3 );
  cell->GetPointIds()->SetId( 4, ptid4 );
  cell->GetPointIds()->SetId( 5, ptid5 );
  cell->GetPointIds()->SetId( 6, ptid6 );
  cell->GetPointIds()->SetId( 7, ptid7 );
  cell->GetPointIds()->SetId( 8, ptid8 );
  cell->GetPointIds()->SetId( 9, ptid9 );  
  cell->GetPointIds()->SetId( 10, ptid10 );
  cell->GetPointIds()->SetId( 11, ptid11 );
  cell->GetPointIds()->SetId( 12, ptid12 );
  cell->GetPointIds()->SetId( 13, ptid13 );
  cell->GetPointIds()->SetId( 14, ptid14 );
  cell->GetPointIds()->SetId( 15, ptid15 );
  cell->GetPointIds()->SetId( 16, ptid16 );
  cell->GetPointIds()->SetId( 17, ptid17 );
  cell->GetPointIds()->SetId( 18, ptid18 );
  cell->GetPointIds()->SetId( 19, ptid19 );
  cell->GetPointIds()->SetId( 20, ptid20 );
  cell->GetPointIds()->SetId( 21, ptid21 );
  cell->GetPointIds()->SetId( 22, ptid22 );
  cell->GetPointIds()->SetId( 23, ptid23 );
  cell->GetPointIds()->SetId( 24, ptid24 );
  cell->GetPointIds()->SetId( 25, ptid25 );
  cell->GetPointIds()->SetId( 26, ptid26 );

  gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
  cell->Delete();
}

void VIS_T::setTetraelem( const int &ptoffset, vtkUnstructuredGrid * gridData )
{
  vtkCell * cell = vtkTetra::New();

  cell->GetPointIds()->SetId( 0, ptoffset + 0 );
  cell->GetPointIds()->SetId( 1, ptoffset + 1 );
  cell->GetPointIds()->SetId( 2, ptoffset + 2 );
  cell->GetPointIds()->SetId( 3, ptoffset + 3 );

  gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
  cell->Delete();
}

void VIS_T::setTetraelem( const int &ptid0, const int &ptid1,
    const int &ptid2, const int &ptid3, vtkUnstructuredGrid * gridData )
{
  vtkCell * cell = vtkTetra::New();

  cell->GetPointIds()->SetId( 0, ptid0 );
  cell->GetPointIds()->SetId( 1, ptid1 );
  cell->GetPointIds()->SetId( 2, ptid2 );
  cell->GetPointIds()->SetId( 3, ptid3 );

  gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
  cell->Delete();
}

void VIS_T::setQuadTetraelem( const int &ptoffset, vtkUnstructuredGrid * gridData )
{
  vtkCell * cell = vtkQuadraticTetra::New();

  cell->GetPointIds()->SetId( 0, ptoffset + 0 );
  cell->GetPointIds()->SetId( 1, ptoffset + 1 );
  cell->GetPointIds()->SetId( 2, ptoffset + 2 );
  cell->GetPointIds()->SetId( 3, ptoffset + 3 );
  cell->GetPointIds()->SetId( 4, ptoffset + 4 );
  cell->GetPointIds()->SetId( 5, ptoffset + 5 );
  cell->GetPointIds()->SetId( 6, ptoffset + 6 );
  cell->GetPointIds()->SetId( 7, ptoffset + 7 );
  cell->GetPointIds()->SetId( 8, ptoffset + 8 );
  cell->GetPointIds()->SetId( 9, ptoffset + 9 );

  gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
  cell->Delete();
}

void VIS_T::setQuadTetraelem( const int &ptid0, const int &ptid1,
    const int &ptid2, const int &ptid3, const int &ptid4, const int &ptid5, 
    const int &ptid6, const int &ptid7, const int &ptid8, const int &ptid9, 
    vtkUnstructuredGrid * gridData )
{
  vtkCell * cell = vtkQuadraticTetra::New();
  
  cell->GetPointIds()->SetId( 0, ptid0 );
  cell->GetPointIds()->SetId( 1, ptid1 );
  cell->GetPointIds()->SetId( 2, ptid2 );
  cell->GetPointIds()->SetId( 3, ptid3 );
  cell->GetPointIds()->SetId( 4, ptid4 );
  cell->GetPointIds()->SetId( 5, ptid5 );
  cell->GetPointIds()->SetId( 6, ptid6 );
  cell->GetPointIds()->SetId( 7, ptid7 );
  cell->GetPointIds()->SetId( 8, ptid8 );
  cell->GetPointIds()->SetId( 9, ptid9 );

  gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
  cell->Delete();
}

std::vector<int> VIS_T::read_epart( const std::string &epart_file, int esize )
{
  std::string fname(epart_file);
  fname.erase( fname.end()-3, fname.end() );
  fname.append(".h5");
  
  const auto elem_part = HDF5_T::read_intVector( fname.c_str(), "/", "part" );

  SYS_T::print_fatal_if( int(elem_part.size()) != esize, "Error: the epart file's part length does not match given size. \n" );

  return elem_part;
}

std::vector<int> VIS_T::readNodeMapping( const std::string &node_mapping_file,
    const char * const &mapping_type )
{
  return HDF5_T::read_intVector( node_mapping_file.c_str(), "/", mapping_type );
}

std::vector<double> VIS_T::readPETSc_vec(const std::string &solution_file_name)
{
  Vec sol_temp;
  VecCreate(PETSC_COMM_SELF, &sol_temp);
  VecSetType(sol_temp, VECSEQ);

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, solution_file_name.c_str(), FILE_MODE_READ, &viewer);
  VecLoad(sol_temp, viewer);
  PetscViewerDestroy(&viewer);

  // Check the sol_temp has correct size
  PetscInt vec_size;
  VecGetSize(sol_temp, &vec_size);

  // read in array
  double * array_temp;
  VecGetArray(sol_temp, &array_temp);

  std::vector<double> sol_vec(vec_size, 0.0);
  for(int ii=0; ii<vec_size; ++ii)
    sol_vec[ii] = array_temp[ii];

  VecRestoreArray(sol_temp, &array_temp);
  VecDestroy(&sol_temp);

  return sol_vec;
}

std::vector<double> VIS_T::readPETSc_vec( const std::string &solution_file_name,
    const std::vector<int> &nodemap, int in_dof )
{
  Vec sol_temp;
  VecCreate(PETSC_COMM_SELF, &sol_temp);
  VecSetType(sol_temp, VECSEQ);

  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_SELF, solution_file_name.c_str(),
      FILE_MODE_READ, &viewer);
  VecLoad(sol_temp, viewer);
  PetscViewerDestroy(&viewer);

  // Check the solution length
  PetscInt vec_size;
  VecGetSize(sol_temp, &vec_size);

  std::vector<double> veccopy(vec_size, 0.0);
  double * array_temp;
  VecGetArray(sol_temp, &array_temp);

  for(int ii=0; ii<vec_size; ++ii)
    veccopy[ii] = array_temp[ii];

  VecRestoreArray(sol_temp, &array_temp);
  VecDestroy(&sol_temp);

  // copy the solution varibles to the correct location
  std::vector<double> sol_vec(vec_size, 0.0);

  // check the nodemap size
  if( (int)nodemap.size() * in_dof != vec_size ) SYS_T::print_fatal("Error: node map size is incompatible with the solution length. \n");

  for(unsigned int ii=0; ii<nodemap.size(); ++ii)
  {
    const int index = nodemap[ii];
    for(int jj=0; jj<in_dof; ++jj)
      sol_vec[in_dof*index+jj] = veccopy[in_dof*ii+jj];
  }

  return sol_vec;
}

// EOF
