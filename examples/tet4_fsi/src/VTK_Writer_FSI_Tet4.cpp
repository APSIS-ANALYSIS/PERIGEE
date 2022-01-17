#include "VTK_Writer_FSI_Tet4.hpp"

VTK_Writer_FSI_Tet4::VTK_Writer_FSI_Tet4( const int &in_nelem,
    const std::string &epart_file ) : nLocBas(4), nElem(in_nelem)
{
  VIS_T::read_epart( epart_file, nElem, epart_map );
}


VTK_Writer_FSI_Tet4::~VTK_Writer_FSI_Tet4()
{
  VEC_T::clean(epart_map);
}

void VTK_Writer_FSI_Tet4::interpolateJ( const int * const &ptid,
    const std::vector<double> &inputData,
    const FEAElement * const &elem,
    vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  std::vector<double> u (nLocBas, 0.0), v (nLocBas, 0.0), w (nLocBas, 0.0);

  for(int ii=0; ii<nLocBas; ++ii)
  {
    u[ii] = inputData[ii*3];
    v[ii] = inputData[ii*3+1];
    w[ii] = inputData[ii*3+2];
  }

  Interpolater intep( nLocBas );

  std::vector<double> ux, uy, uz, vx, vy, vz, wx, wy, wz;

  intep.interpolateFE_Grad(u, elem, ux, uy, uz);
  intep.interpolateFE_Grad(v, elem, vx, vy, vz);
  intep.interpolateFE_Grad(w, elem, wx, wy, wz);

  for(int ii=0; ii<nqp; ++ii)
  {
    Matrix_3x3 F( ux[ii] + 1.0, uy[ii],       uz[ii],
                  vx[ii],       vy[ii] + 1.0, vz[ii],
                  wx[ii],       wy[ii],       wz[ii] + 1.0 );

    vtkData->InsertComponent( ptid[ii], 0, F.det() );
  }
}

void VTK_Writer_FSI_Tet4::writeOutput(
        const FEANode * const &fnode_ptr,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const ALocal_Elem * const &lelem_ptr,
        const IVisDataPrep * const &vdata_ptr,
        FEAElement * const &elemptr,
        const IQuadPts * const &quad,
        const double * const * const &pointArrays,
        const int &rank, const int &size,
        const int &num_of_nodes,
        const double &sol_time,
        const std::string &basename,
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML )
{
  // This routine requires nqp = 4
  SYS_T::print_fatal_if(quad->get_num_quadPts() != 4, "Error: VTK_Writer_Tet4 requires 4 quadrature points.\n");

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
  if(numDArrays != 4) SYS_T::print_fatal("Error: vdata size numDArrays != 4.\n");

  // dataVecs is the holder of the interpolated data
  vtkDoubleArray ** dataVecs = new vtkDoubleArray * [numDArrays];
  for(int ii=0; ii<numDArrays; ++ii)
  {
    dataVecs[ii] = vtkDoubleArray::New();
    dataVecs[ii] -> SetNumberOfComponents( vdata_ptr->get_arraySizes(ii) );
    std::string temp_name = vdata_ptr->get_arrayNames(ii);
    dataVecs[ii] -> SetName( temp_name.c_str() );
    dataVecs[ii] -> SetNumberOfTuples( num_of_nodes );
  }

  // An additional array for the analysis partition info
  vtkIntArray * anaprocId = vtkIntArray::New();
  anaprocId -> SetName("Analysis_Partition");
  anaprocId -> SetNumberOfComponents(1);

  // An additional array for the FSI physical sub-domain info
  vtkIntArray * phyDomain = vtkIntArray::New();
  phyDomain -> SetName("Physical_subdomain");
  phyDomain -> SetNumberOfComponents(1);

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    const std::vector<int> IEN_v = lien_v -> get_LIEN( ee );
    const std::vector<int> IEN_p = lien_p -> get_LIEN( ee );

    double ectrl_x[4], ectrl_y[4], ectrl_z[4];

    fnode_ptr -> get_ctrlPts_xyz(nLocBas, &IEN_v[0], ectrl_x, ectrl_y, ectrl_z);

    elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

    // Interpolate data and assign to dataVecs
    std::vector<double> inputInfo; inputInfo.clear();
    int asize = vdata_ptr -> get_arraySizes(0);

    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_v[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[0][pt_index * asize + kk] );
    }

    // displacement interpolation
    intep.interpolateVTKData( asize, &IEN_p[0], inputInfo, elemptr, dataVecs[0] );

    // use displacement to update points
    intep.interpolateVTKPts( &IEN_p[0], ectrl_x, ectrl_y, ectrl_z, inputInfo, elemptr, points);

    // Interpolate detF
    interpolateJ( &IEN_p[0], inputInfo, elemptr, dataVecs[1] );

    // Interpolate the pressure scalar
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes(2);
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_p[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, &IEN_p[0], inputInfo, elemptr, dataVecs[2] );

    // Interpolate the velocity vector
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes(3);
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_v[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, &IEN_p[0], inputInfo, elemptr, dataVecs[3] );

    // Set mesh connectivity
    VIS_T::setTetraelem( IEN_p[0], IEN_p[1], IEN_p[2], IEN_p[3], gridData ); 

    // Analysis mesh partition
    const int e_global = lelem_ptr->get_elem_loc(ee);
    anaprocId->InsertNextValue( epart_map[e_global] );

    // Physical domain info
    phyDomain->InsertNextValue( lelem_ptr->get_elem_tag(ee) );
  }

  gridData -> SetPoints( points );
  points -> Delete();

  for(int ii=0; ii<numDArrays; ++ii)
  {
    gridData->GetPointData()->AddArray( dataVecs[ii] );
    dataVecs[ii]->Delete();
  }

  gridData->GetCellData()->AddArray(anaprocId);
  gridData->GetCellData()->AddArray(phyDomain);

  delete [] dataVecs;
  anaprocId->Delete();
  phyDomain->Delete();

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
  VIS_T::writeVisFile( gridData, outputBName, outputName, rank, size, sol_time, isXML );

  // Clean gridData
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object.\n");
  gridData->Delete();
}

// EOF
