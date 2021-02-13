#include "VTK_Writer.hpp"

VTK_Writer::VTK_Writer( const IAGlobal_Mesh_Info * const &gmesh_ptr )
{
  nLocBas = gmesh_ptr->get_nLocBas();
  
  R = new double [nLocBas];
}

VTK_Writer::~VTK_Writer()
{
  delete [] R; R = nullptr;
}

void VTK_Writer::interpolateNURBS( const double * const &inputVal,
    const FEAElement * const &elem, double * const &output ) const
{
  const int nqp = elem->get_numQuapts();

  for(int ii=0; ii<nqp; ++ii)
  {
    output[ii] = 0.0;
    elem->get_R(ii, R);
    for(int jj=0; jj<nLocBas; ++jj)
      output[ii] += inputVal[jj] * R[jj];
  }
}

void VTK_Writer::interpolateNURBS(const double * const &inputVal_1,
    const double * const &inputVal_2, const FEAElement * const &elem,
    double * const &output_1, double * const &output_2 ) const
{
  const int nqp = elem->get_numQuapts();

  for(int ii=0; ii<nqp; ++ii)
  {
    output_1[ii] = 0.0;
    output_2[ii] = 0.0;
    elem->get_R(ii, R);
    for(int jj=0; jj<nLocBas; ++jj)
    {
      output_1[ii] += inputVal_1[jj] * R[jj];
      output_2[ii] += inputVal_2[jj] * R[jj];
    }
  }
}

void VTK_Writer::interpolateNURBS(const double * const &inputVal_1,
    const double * const &inputVal_2, const double * const &inputVal_3,
    const FEAElement * const &elem, double * const &output_1, 
    double * const &output_2, double * const &output_3 ) const
{
  const int nqp = elem->get_numQuapts();

  for(int ii=0; ii<nqp; ++ii)
  {
    output_1[ii] = 0.0;
    output_2[ii] = 0.0;
    output_3[ii] = 0.0;
    elem->get_R(ii, R);
    for(int jj=0; jj<nLocBas; ++jj)
    {
      output_1[ii] += inputVal_1[jj] * R[jj];
      output_2[ii] += inputVal_2[jj] * R[jj];
      output_3[ii] += inputVal_3[jj] * R[jj];
    }
  }
}

void VTK_Writer::interpolatePts( 
    const int &ptoffset, const double * const &ctrlPts_x,
    const double * const &ctrlPts_y, const double * const &ctrlPts_z,
    const FEAElement * const &elem, vtkPoints * const &vtkpts ) const
{
  const int nqp = elem->get_numQuapts();
  double * out_x = new double [nqp];
  double * out_y = new double [nqp];
  double * out_z = new double [nqp];

  interpolateNURBS(ctrlPts_x, ctrlPts_y, ctrlPts_z, elem, out_x, out_y, out_z);
  
  double xyz[3];

  for(int ii=0; ii<nqp; ++ii)
  {
    xyz[0] = out_x[ii];
    xyz[1] = out_y[ii];
    xyz[2] = out_z[ii];

    vtkpts->InsertPoint(ptoffset+ii, xyz);
  }
  delete [] out_x; delete [] out_y; delete [] out_z;
}

void VTK_Writer::interpolatePts( 
    const int &ptoffset, const double * const &ctrlPts_x,
    const double * const &ctrlPts_y, const FEAElement * const &elem,
    vtkPoints * const &vtkpts ) const
{
  const int nqp = elem->get_numQuapts();
  double * out_x = new double [nqp];
  double * out_y = new double [nqp];

  interpolateNURBS(ctrlPts_x, ctrlPts_y, elem, out_x, out_y);
  
  double xyz[2];

  for(int ii=0; ii<nqp; ++ii)
  {
    xyz[0] = out_x[ii];
    xyz[1] = out_y[ii];

    vtkpts->InsertPoint(ptoffset+ii, xyz);
  }
  delete [] out_x; delete [] out_y;
}

void VTK_Writer::interpolateData( const int &size, const int &ptoffset,
    const double * const &inputData, const FEAElement * const &elem,
    vtkDoubleArray * const &vtkData ) const
{
  const int nqp = elem->get_numQuapts();

  double * compData = new double [nLocBas];
  double * outData  = new double [nqp];

  for(int comp=0; comp<size; ++comp)
  {
    for(int ii=0; ii<nLocBas; ++ii)
      compData[ii] = inputData[ii*size + comp]; // extract the comp-th component
    interpolateNURBS(compData, elem, outData);  // interpolate it over the element
    
    // Insert the comp-th component into vtkData
    for(int ii=0; ii<nqp; ++ii)
      vtkData->InsertComponent(ptoffset+ii, comp, outData[ii]);
  }

  delete [] compData; delete [] outData;
}

void VTK_Writer::setHexelem( const int &segs, const int &segt, 
    const int &segu, const int &ptoffset, 
    vtkUnstructuredGrid * gridData ) const
{
  for(int ii=0; ii<segs-1; ++ii)
  {
    for(int jj=0; jj<segt-1; ++jj)
    {
      for(int kk=0; kk<segu-1; ++kk)
      {
        vtkCell * cell = vtkHexahedron::New();
        
        cell->GetPointIds()->SetId( 0, ptoffset + kk*segs*segt + jj*segs + ii);
        cell->GetPointIds()->SetId( 1, ptoffset + kk*segs*segt + jj*segs + ii + 1);
        cell->GetPointIds()->SetId( 2, ptoffset + kk*segs*segt + (jj+1)*segs + ii + 1);
        cell->GetPointIds()->SetId( 3, ptoffset + kk*segs*segt + (jj+1)*segs + ii);

        cell->GetPointIds()->SetId( 4, ptoffset + (kk+1)*segs*segt + jj*segs + ii);
        cell->GetPointIds()->SetId( 5, ptoffset + (kk+1)*segs*segt + jj*segs + ii + 1);
        cell->GetPointIds()->SetId( 6, ptoffset + (kk+1)*segs*segt + (jj+1)*segs + ii + 1);
        cell->GetPointIds()->SetId( 7, ptoffset + (kk+1)*segs*segt + (jj+1)*segs + ii);
       
        gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );

        cell->Delete();
      }
    }
  }
}

void VTK_Writer::setQuadelem( const int &segs, const int &segt, 
    const int &ptoffset, vtkUnstructuredGrid * gridData ) const
{
  for(int ii=0; ii<segs-1; ++ii)
  {
    for(int jj=0; jj<segt-1; ++jj)
    {
      vtkCell * cell = vtkQuad::New();
      cell->GetPointIds()->SetId(0, ptoffset + jj*segs + ii);
      cell->GetPointIds()->SetId(1, ptoffset + jj*segs + ii + 1);
      cell->GetPointIds()->SetId(2, ptoffset + (jj+1)*segs + ii + 1);
      cell->GetPointIds()->SetId(3, ptoffset + (jj+1)*segs + ii);

      gridData->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
      cell->Delete();
    }
  }
}

void VTK_Writer::read_epart( const std::string &epart_file, const int &esize,
    int * const &elem_part ) const
{
  hid_t file_id = H5Fopen(epart_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  hid_t data_id = H5Dopen(file_id, "part", H5P_DEFAULT);
  hid_t data_space = H5Dget_space( data_id );
  hid_t data_rank = H5Sget_simple_extent_ndims( data_space );

  if( data_rank != 1 )
  {
    PetscPrintf(PETSC_COMM_WORLD, 
        "Error: the epart file's part is not one dimensional. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  hsize_t * data_dims = new hsize_t [1];
  H5Sget_simple_extent_dims( data_space, data_dims, NULL );

  hid_t mem_space = H5Screate_simple(data_rank, data_dims, NULL);

  if( int(data_dims[0]) != esize )
  {
    PetscPrintf(PETSC_COMM_WORLD, 
        "Error: the epart file's part length does not match given size. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  H5Dread( data_id, H5T_NATIVE_INT, mem_space, data_space,
      H5P_DEFAULT, elem_part );

  delete [] data_dims;
  H5Sclose( mem_space );
  H5Sclose( data_space );
  H5Dclose( data_id );
  H5Fclose( file_id );
}

void VTK_Writer::build3doutput(vtkUnstructuredGrid * const &gridData,
    const std::string &epart_file, const int nElem,
    const FEANode * const &fNode_ptr,
    const IALocal_meshSize * const &lmsize_ptr,
    const IAExtractor * const &ext_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    const BernsteinBasis_Array * const &bern_x,
    const BernsteinBasis_Array * const &bern_y,
    const BernsteinBasis_Array * const &bern_z,
    const double * const * const &pointArrays ) const
{
  // point struct initialized
  vtkPoints * points = vtkPoints::New();

  // vtkDoubleArray & vtkIntArray allocation
  const int numDArrays = vdata_ptr->get_arrayCompSize();

  vtkDoubleArray ** dataVecs = new vtkDoubleArray * [numDArrays];
  for(int ii=0; ii<numDArrays; ++ii)
  {
    dataVecs[ii] = vtkDoubleArray::New();
    dataVecs[ii] -> SetNumberOfComponents( vdata_ptr->get_arraySizes(ii) );
    std::string temp_name = vdata_ptr->get_arrayNames(ii);
    dataVecs[ii] -> SetName( temp_name.c_str() );
  }

  vtkIntArray * anaprocId = vtkIntArray::New();
  anaprocId -> SetName("Analysis_Partition");
  anaprocId -> SetNumberOfComponents(1);

  // read in epart file to get analysis element partition
  int * epart_map = new int [nElem];
  read_epart(epart_file, nElem, epart_map);

  // allocate memory for some arrays needed in for loop
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];

  int ptOffset = 0;
  // loop over element to interpolate data and coordinates
  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    FEAElement * elem_ptr = new FEAElement_NURBS_3D_der0_v3(ee, lmsize_ptr,
        bern_x, bern_y, bern_z, fNode_ptr, ext_ptr, lien_ptr );

    if(elem_ptr->is_sizeNonzero())
    {
      // number of sampling points in each direction:
      int sSegs   = bern_x->get_nQuapts();
      int tSegs   = bern_y->get_nQuapts();
      int uSegs   = bern_z->get_nQuapts();

      lien_ptr->get_LIEN_e(ee, IEN_e);

      fNode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      // interpolate for the xyz coordinates of sampling points
      interpolatePts(ptOffset, ectrl_x, ectrl_y, ectrl_z, elem_ptr, points );

      // interpolate for data to assign dataVecs
      for(int ii=0; ii<numDArrays; ++ii)
      {
        std::vector<double> inputInfo;
        inputInfo.clear();
        int asize = vdata_ptr->get_arraySizes(ii);
        for(int jj=0; jj<nLocBas; ++jj)
        {
          int pt_index = IEN_e[jj];

          for(int kk=0; kk<asize; ++kk)
            inputInfo.push_back( pointArrays[ii][pt_index * asize + kk ] );
        }
        interpolateData(asize, ptOffset, &inputInfo[0], elem_ptr, dataVecs[ii]);
      }
      // build connectivity for sampling points' sub-elements
      setHexelem(sSegs, tSegs, uSegs, ptOffset, gridData);

      // record this element's analysis partition index
      int e_global = lelem_ptr->get_elem_loc(ee);
      for(int ii=0; ii<(sSegs-1)*(tSegs-1)*(uSegs-1); ++ii)
        anaprocId->InsertNextValue(epart_map[e_global]);

      // update the sampling point offset for the next element
      ptOffset += sSegs * tSegs * uSegs;
    }

    delete elem_ptr;
  }

  // pass build info into gridData
  gridData->SetPoints(points);
  points->Delete();

  for(int ii=0; ii<numDArrays; ++ii)
  {
    gridData->GetPointData()->AddArray( dataVecs[ii] );
    dataVecs[ii]->Delete();
  }
  gridData->GetCellData()->AddArray(anaprocId);

  // free memory
  delete [] dataVecs;
  anaprocId->Delete();
  delete [] IEN_e; delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z;
  delete [] epart_map;
}

void VTK_Writer::build2doutput(vtkUnstructuredGrid * const &gridData,
    const std::string &epart_file, const int nElem,
    const FEANode * const &fNode_ptr,
    const IALocal_meshSize * const &lmsize_ptr,
    const IAExtractor * const &ext_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    const BernsteinBasis_Array * const &bern_x,
    const BernsteinBasis_Array * const &bern_y,
    const double * const * const &pointArrays ) const
{
  // point struct initialized
  vtkPoints * points = vtkPoints::New();

  // vtkDoubleArray & vtkIntArray allocation
  const int numDArrays = vdata_ptr->get_arrayCompSize();

  vtkDoubleArray ** dataVecs = new vtkDoubleArray * [numDArrays];
  for(int ii=0; ii<numDArrays; ++ii)
  {
    dataVecs[ii] = vtkDoubleArray::New();
    dataVecs[ii] -> SetNumberOfComponents( vdata_ptr->get_arraySizes(ii) );
    std::string temp_name = vdata_ptr->get_arrayNames(ii);
    dataVecs[ii] -> SetName( temp_name.c_str() );
  }

  vtkIntArray * anaprocId = vtkIntArray::New();
  anaprocId -> SetName("Analysis_Partition");
  anaprocId -> SetNumberOfComponents(1);

  // read in epart file to get analysis element partition
  int * epart_map = new int [nElem];
  read_epart(epart_file, nElem, epart_map);

  // allocate memory for some arrays needed in for loop
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];

  int ptOffset = 0;
  // loop over element to interpolate data and coordinates
  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    FEAElement * elem_ptr = new FEAElement_NURBS_2D_der0(ee, lmsize_ptr,
        bern_x, bern_y, fNode_ptr, ext_ptr, lien_ptr );

    if(elem_ptr->is_sizeNonzero())
    {
      // number of sampling points in each direction:
      int sSegs   = bern_x->get_nQuapts();
      int tSegs   = bern_y->get_nQuapts();

      lien_ptr->get_LIEN_e(ee, IEN_e);

      fNode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      // interpolate for the xyz coordinates of sampling points
      interpolatePts(ptOffset, ectrl_x, ectrl_y, elem_ptr, points );

      // interpolate for data to assign dataVecs
      for(int ii=0; ii<numDArrays; ++ii)
      {
        std::vector<double> inputInfo;
        inputInfo.clear();
        int asize = vdata_ptr->get_arraySizes(ii);
        for(int jj=0; jj<nLocBas; ++jj)
        {
          int pt_index = IEN_e[jj];

          for(int kk=0; kk<asize; ++kk)
            inputInfo.push_back( pointArrays[ii][pt_index * asize + kk ] );

          interpolateData(asize, ptOffset, &inputInfo[0], elem_ptr, dataVecs[ii]);
        }
      }
      // build connectivity for sampling points' sub-elements
      setQuadelem(sSegs, tSegs, ptOffset, gridData);

      // record this element's analysis partition index
      int e_global = lelem_ptr->get_elem_loc(ee);
      for(int ii=0; ii<(sSegs-1)*(tSegs-1); ++ii)
        anaprocId->InsertNextValue(epart_map[e_global]);

      // update the sampling point offset for the next element
      ptOffset += sSegs * tSegs;
    }

    delete elem_ptr;
  }

  // pass build info into gridData
  gridData->SetPoints(points);
  points->Delete();

  for(int ii=0; ii<numDArrays; ++ii)
  {
    gridData->GetPointData()->AddArray( dataVecs[ii] );
    dataVecs[ii]->Delete();
  }
  gridData->GetCellData()->AddArray(anaprocId);

  // free memory
  delete [] dataVecs;
  anaprocId->Delete();
  delete [] IEN_e; delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z;
  delete [] epart_map;
}

void VTK_Writer::writeOutput( const std::string &epart_file,
    const IAGlobal_Mesh_Info * const &gmesh_ptr,
    const FEANode * const &fnode_ptr,
    const IALocal_meshSize * const &lmsize_ptr,
    const IAExtractor * const &ext_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    const BernsteinBasis_Array * const &bern_x,
    const BernsteinBasis_Array * const &bern_y,
    const BernsteinBasis_Array * const &bern_z,
    const double * const * const &pointArrays,
    const int &rank, const int &size,
    const double &sol_time,
    const std::string &basename,
    const std::string &outputBName,
    const std::string &outputName,
    const bool &isXML ) const
{
  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();

  // build the vtk gridData for writer
  build3doutput(gridData, epart_file, gmesh_ptr->get_nElem(),
      fnode_ptr, lmsize_ptr, ext_ptr, lien_ptr, lelem_ptr,
      vdata_ptr, bern_x, bern_y, bern_z, pointArrays );

  // add visualization id if run in parallel
  if(size > 1)
  {
    int numCells = gridData->GetNumberOfCells();
    vtkIntArray * procId = vtkIntArray::New();
    procId->SetName("PostProcess_ID");
    procId->SetNumberOfComponents(1);
    for(int ii=0; ii<numCells; ++ii)
      procId->InsertComponent(ii, 0 , rank);
    gridData->GetCellData()->AddArray(procId);
    procId->Delete();
  }

  // write the vtk/vtu files
  PetscPrintf(PETSC_COMM_WORLD, "-- Write vtu files ... \n");
  writeGridObject(gridData, outputName, rank, size, isXML);

  // write the pvtu file
  writepvtuFile(gridData, outputName, rank, size, isXML);

  MPI_Barrier(PETSC_COMM_WORLD);

  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object ... \n");
  gridData->Delete();

  // write the pvd file
  PetscPrintf(PETSC_COMM_WORLD, "-- Update pvd file. \n");
  writepvdFile(outputBName, outputName, rank, size, sol_time, isXML);

  MPI_Barrier(PETSC_COMM_WORLD);
}

void VTK_Writer::writeOutput( const std::string &epart_file,
    const IAGlobal_Mesh_Info * const &gmesh_ptr,
    const FEANode * const &fnode_ptr,
    const IALocal_meshSize * const &lmsize_ptr,
    const IAExtractor * const &ext_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    const BernsteinBasis_Array * const &bern_x,
    const BernsteinBasis_Array * const &bern_y,
    const double * const * const &pointArrays,
    const int &rank, const int &size,
    const double &sol_time,
    const std::string &basename,
    const std::string &outputBName,
    const std::string &outputName,
    const bool &isXML ) const
{
  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();

  // build the vtk gridData for writer
  build2doutput(gridData, epart_file, gmesh_ptr->get_nElem(),
      fnode_ptr, lmsize_ptr, ext_ptr, lien_ptr, lelem_ptr,
      vdata_ptr, bern_x, bern_y, pointArrays );

  // add visualization id if run in parallel
  if(size > 1)
  {
    int numCells = gridData->GetNumberOfCells();
    vtkIntArray * procId = vtkIntArray::New();
    procId->SetName("PostProcess_ID");
    procId->SetNumberOfComponents(1);
    for(int ii=0; ii<numCells; ++ii)
      procId->InsertComponent(ii, 0 , rank);
    gridData->GetCellData()->AddArray(procId);
    procId->Delete();
  }

  // write the vtk/vtu files
  PetscPrintf(PETSC_COMM_WORLD, "-- Write vtu files ... \n");
  writeGridObject(gridData, outputName, rank, size, isXML);

  // write the pvtu file
  writepvtuFile(gridData, outputName, rank, size, isXML);

  MPI_Barrier(PETSC_COMM_WORLD);

  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object ... \n");
  gridData->Delete();

  // write the pvd file
  PetscPrintf(PETSC_COMM_WORLD, "-- Update pvd file. \n");
  writepvdFile(outputBName, outputName, rank, size, sol_time, isXML);

  MPI_Barrier(PETSC_COMM_WORLD);
}

void VTK_Writer::writeGridObject( vtkUnstructuredGrid * gridData,
    const std::string &baseName, const int &rank, const int &size,
    const bool &isXML ) const
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

void VTK_Writer::writepvtuFile( vtkUnstructuredGrid * gridData,
    const std::string &baseName, const int &rank, const int &size,
    const bool &isXML) const
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
      std::cerr<<"ERROR: failed to open pvtu file "<<pName<<std::endl;
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
    pvtuFile << "<PCellData>" << std::endl
      << "<PDataArray type=\"Int32\" Name=\"Analysis_Partition\""
      << " NumberOfComponents=\"1\"/>" << std::endl
      << "<PDataArray type=\"Int32\" Name=\"PostProcess_ID\""
      << " NumberOfComponents=\"1\"/>" << std::endl
      << "</PCellData>" << std::endl;


    // Point Data: associated with the gridData object
    vtkPointData* pointData = gridData->GetPointData();
    pvtuFile << "<PPointData>" << std::endl;
    int ii;
    int numArrays = pointData->GetNumberOfArrays();

    for( ii = 0; ii < numArrays; ii++ )
    {
      vtkDataArray* currArray = pointData->GetArray( ii );
      if( NULL == currArray )
        continue;
      pvtuFile << "<PDataArray type=\"";
      if( VTK_DOUBLE == currArray->GetDataType()  )
        pvtuFile << "Float64";
      else
        assert( 0 );

      pvtuFile << "\" Name=\"" << currArray->GetName() << "\""
        << " NumberOfComponents=\"" << currArray->GetNumberOfComponents() << "\"/>"
        << std::endl;
    }

    pvtuFile << "</PPointData>" << std::endl;

    pvtuFile << "<PPoints>" << std::endl
      << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>" << std::endl
      << "</PPoints>" << std::endl;

    for(ii=0; ii<size; ++ii)
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

void VTK_Writer::writepvdFile_Init( const std::string &pvdFName ) const
{  
  std::ofstream pvdFile;

  pvdFile.open( pvdFName.c_str(), std::ios::out | std::ios::trunc );

  if( !pvdFile.is_open() )
  {
    std::cerr<<"ERROR: UNABLE TO OPEN THE PVD FILE : "<<pvdFName<<std::endl;
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

void VTK_Writer::writepvdFile( const std::string &baseName,  
    const std::string &pre_pvtuname, const int &rank, const int &size,
    const double &sol_time, const bool &isXML ) const
{
  if( rank == 0 && isXML == true )
  {
    std::string pvdname( baseName );
    pvdname.append(".pvd");

    // if baseName.pvd does not exists, write the file with top initialized part
    fstream fileTest;
    fileTest.open(pvdname.c_str(), fstream::in);
    if( !fileTest.is_open() )
      writepvdFile_Init(pvdname);

    fileTest.close();

    // now we are sure the pvd file exists, we just need to append the pvtu and
    // time information to it.
    fstream pvdFile;
    pvdFile.open(pvdname.c_str());
    if( !pvdFile.is_open() )
    {
      std::cerr<<"ERROR: FAIL TO OPEN THE PVD FILE "<<pvdname<<std::endl;
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

// EOF
