#include "VTK_Writer_Fluids_ALE_Tet4.hpp"

VTK_Writer_Fluids_ALE_Tet4::VTK_Writer_Fluids_ALE_Tet4( 
    const int &in_nelem,
    const std::string &epart_file )
: nLocBas(4), nElem(in_nelem)
{
  VIS_T::read_epart( epart_file, nElem, epart_map );
}

VTK_Writer_Fluids_ALE_Tet4::~VTK_Writer_Fluids_ALE_Tet4()
{
  VEC_T::clean(epart_map);
}

void VTK_Writer_Fluids_ALE_Tet4::writeOutput(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    FEAElement * const &elemptr,
    const IQuadPts * const &quad,
    const double * const * const &pointArrays,
    const int &rank, const int &size,
    const double &sol_time,
    const std::string &basename,
    const std::string &outputBName,
    const std::string &outputName,
    const bool &isXML )
{
  Interpolater intep( nLocBas );

  // Allocate gridData
  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();

  // Now generate and store data to gridData
  // Point object initialized
  vtkPoints * points = vtkPoints::New();

  // DoubleArray & IntArray allocation
  const int numDArrays = vdata_ptr->get_arrayCompSize();

  // Check if the number of visualization quantities equal a given number
  // I will have to manually interpolate these quantities in the for-loop
  if(numDArrays != 3) SYS_T::print_fatal("Error: vdata size numDArrays != 3. \n");

  // dataVecs is the holder of the interpolated data
  vtkDoubleArray ** dataVecs = new vtkDoubleArray * [numDArrays];
  for(int ii=0; ii<numDArrays; ++ii)
  {
    dataVecs[ii] = vtkDoubleArray::New();
    dataVecs[ii] -> SetNumberOfComponents( vdata_ptr->get_arraySizes(ii) );
    std::string temp_name = vdata_ptr->get_arrayNames(ii);
    dataVecs[ii] -> SetName( temp_name.c_str() );
  }

  // An additional array for the analysis partition info
  vtkIntArray * anaprocId = vtkIntArray::New();
  anaprocId -> SetName("Analysis_Partition");
  anaprocId -> SetNumberOfComponents(1);

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    const std::vector<int> IEN_e = lien_ptr -> get_LIEN(ee);

    double ectrl_x[4]; double ectrl_y[4]; double ectrl_z[4];

    fnode_ptr -> get_ctrlPts_xyz(nLocBas, &IEN_e[0], ectrl_x, ectrl_y, ectrl_z);

    elemptr -> buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

    // Interpolate data and assign to dataVecs
    int visCompIdx = 0;
    std::vector<double> inputInfo; inputInfo.clear();
    int asize = vdata_ptr -> get_arraySizes( visCompIdx );

    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[0][pt_index * asize + kk] );
    }

    // displacement interpolation
    intep.interpolateVTKData( asize, &IEN_e[0], &inputInfo[0],
        elemptr, dataVecs[visCompIdx] );

    // use displacement to update points
    intep.interpolateVTKPts( &IEN_e[0], ectrl_x, ectrl_y, ectrl_z,
        &inputInfo[0], elemptr, points ); 

    // Interpolate the pressure scalar
    visCompIdx += 1;
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes( visCompIdx );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, &IEN_e[0], &inputInfo[0],
        elemptr, dataVecs[visCompIdx] );

    // Interpolate the velocity vector
    inputInfo.clear();
    visCompIdx += 1;
    asize = vdata_ptr->get_arraySizes( visCompIdx );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, &IEN_e[0], &inputInfo[0],
        elemptr, dataVecs[visCompIdx] );

    // Set mesh connectivity
    VIS_T::setTetraelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3], gridData );

    // Analysis mesh partition 
    int e_global = lelem_ptr->get_elem_loc(ee);
    anaprocId->InsertNextValue( epart_map[e_global] );
  }

  gridData -> SetPoints( points );
  points -> Delete();

  for(int ii=0; ii<numDArrays; ++ii)
  {
    gridData->GetPointData()->AddArray( dataVecs[ii] );
    dataVecs[ii]->Delete();
  }

  gridData->GetCellData()->AddArray(anaprocId);

  delete [] dataVecs;
  anaprocId->Delete();

  // If postprocess is parallel, record its partition
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

  // Write gridData
  VIS_T::writeVisFile( gridData, outputBName, outputName,
      rank, size, sol_time, isXML );

  // Clean gridData
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object.\n");
  gridData->Delete();
}

// EOF
