#include "VTK_Writer_Solids_Tet4.hpp"

VTK_Writer_Solids_Tet4::VTK_Writer_Solids_Tet4( const int &in_nelem,
    const std::string &epart_file, const int &in_nLocBas )
: nLocBas( in_nLocBas ), nElem(in_nelem), intep(nLocBas, true)
{
  VIS_T::read_epart( epart_file, nElem, epart_map );

  IEN_e   = new int [nLocBas];
  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];
}


VTK_Writer_Solids_Tet4::~VTK_Writer_Solids_Tet4()
{
  VEC_T::clean(epart_map);

  delete IEN_e; delete ectrl_x; delete ectrl_y; delete ectrl_z;
}


void VTK_Writer_Solids_Tet4::writeOutput_cur_isotropic(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    IMaterialModel * const &model,
    FEAElement * const &elemptr,
    const IQuadPts * const &quad,
    const double * const * const &pointArrays,
    const int &rank, const int &size,
    const double &sol_time,
    const std::string &outputBName,
    const std::string &outputName,
    const bool &isXML )
{
  // Allocate gridData
  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();

  // Now generate and store data to gridData
  // Point object initialized
  vtkPoints * points = vtkPoints::New();

  // DoubleArray & IntArray allocation
  const int numDArrays = vdata_ptr->get_arrayCompSize();

  // Check if the number of visualization quantities equal a given number
  // I will have to manually interpolate these quantities in the for-loop
  if(numDArrays != 5) SYS_T::print_fatal("Error: vdata size numDArrays != 5. \n");

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

  int ptOffset = 0;
  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    lien_ptr -> get_LIEN_e(ee, IEN_e);

    fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

    // Interpolate data and assign to dataVecs
    std::vector<double> inputInfo; inputInfo.clear();
    int asize = vdata_ptr -> get_arraySizes( 0 );

    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[0][pt_index * asize + kk] );
    }

    // displacement interpolation
    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[0] );

    // use displacement to update points
    intep.interpolateVTKPts(ptOffset, ectrl_x, ectrl_y, ectrl_z,
        &inputInfo[0], elemptr, points ); 

    // Interpolate detF
    interpolateJ( ptOffset, &inputInfo[0], elemptr, dataVecs[1] );
    
    // Interpolate von Mises
    interpolateVonMise( ptOffset, &inputInfo[0], elemptr, model, dataVecs[2] );

    // Interpolate the pressure scalar
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes( 3 );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
    }

    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[3] );

    // Interpolate the velocity vector
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes( 4 );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[4] );

    // Set mesh connectivity
    if( nLocBas == 4)
      VIS_T::setTetraelem( ptOffset, gridData );
    else if(nLocBas == 10)
      VIS_T::setQuadTetraelem( ptOffset, gridData );

    // Analysis mesh partition 
    int e_global = lelem_ptr->get_elem_loc(ee);
    anaprocId->InsertNextValue( epart_map[e_global] );

    // update offset
    ptOffset += nLocBas;
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
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object ... \n");
  gridData->Delete();
}


void VTK_Writer_Solids_Tet4::writeOutput_cur_anisotropic(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    IMaterialModel * const &model,
    FEAElement * const &elemptr,
    const IQuadPts * const &quad,
    const double * const * const &pointArrays,
    const int &rank, const int &size,
    const double &sol_time,
    const std::string &outputBName,
    const std::string &outputName,
    const bool &isXML )
{
  // Allocate gridData
  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();

  // Now generate and store data to gridData
  // Point object initialized
  vtkPoints * points = vtkPoints::New();

  // DoubleArray & IntArray allocation
  const int numDArrays = vdata_ptr->get_arrayCompSize();

  // Check if the number of visualization quantities equal a given number
  // I will have to manually interpolate these quantities in the for-loop
  if(numDArrays != 7) SYS_T::print_fatal("Error: vdata size numDArrays != 7. \n");

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

  int ptOffset = 0;
  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    lien_ptr -> get_LIEN_e(ee, IEN_e);

    fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

    // Interpolate data and assign to dataVecs
    std::vector<double> inputInfo; inputInfo.clear();
    int asize = vdata_ptr -> get_arraySizes( 0 );

    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[0][pt_index * asize + kk] );
    }

    std::vector<double> inputInfo_p; inputInfo_p.clear();
    for(int jj=0; jj<nLocBas; ++jj)
      inputInfo_p.push_back( pointArrays[1][ IEN_e[jj] ] );
    
    // displacement interpolation
    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[0] );

    // use displacement to update points
    intep.interpolateVTKPts(ptOffset, ectrl_x, ectrl_y, ectrl_z,
        &inputInfo[0], elemptr, points ); 

    // Interpolate detF
    interpolateJ( ptOffset, &inputInfo[0], elemptr, dataVecs[1] );

    interpolateJele( ptOffset, &inputInfo[0], elemptr, dataVecs[2] );

    interpolateCauchyZpPres( ptOffset, &inputInfo[0], &inputInfo_p[0], 
        elemptr, model, dataVecs[3] );
    
    //interpolateCauchyZ( ptOffset, &inputInfo[0], elemptr, model, dataVecs[3] );
    
    interpolateOrientation( ptOffset, &inputInfo[0], elemptr, model, dataVecs[4] );

    // Interpolate the pressure scalar
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes( 5 );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
    }

    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[5] );

    // Interpolate the velocity vector
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes( 6 );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[6] );

    // Set mesh connectivity
    VIS_T::setTetraelem( ptOffset, gridData );

    // Analysis mesh partition 
    int e_global = lelem_ptr->get_elem_loc(ee);
    anaprocId->InsertNextValue( epart_map[e_global] );

    // update offset
    ptOffset += 4;
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
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object ... \n");
  gridData->Delete();
}


void VTK_Writer_Solids_Tet4::writeOutput_ref(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    IMaterialModel * const &model,
    FEAElement * const &elemptr,
    const IQuadPts * const &quad,
    const double * const * const &pointArrays,
    const int &rank, const int &size,
    const double &sol_time,
    const std::string &outputBName,
    const std::string &outputName,
    const bool &isXML )
{ 
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

  int ptOffset = 0;
  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    lien_ptr -> get_LIEN_e(ee, IEN_e);

    fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

    // Interpolate data and assign to dataVecs
    std::vector<double> inputInfo; inputInfo.clear();
    int asize = vdata_ptr -> get_arraySizes( 0 );

    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[0][pt_index * asize + kk] );
    }

    // displacement interpolation
    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[0] );

    // interpolate points
    intep.interpolateVTKPts(ptOffset, ectrl_x, ectrl_y, ectrl_z,
        elemptr, points ); 

    // Interpolate the pressure scalar
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes( 1 );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
    }

    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[1] );

    // Interpolate the velocity vector
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes( 2 );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[2] );

    // Set mesh connectivity
    VIS_T::setTetraelem( ptOffset, gridData );

    // Analysis mesh partition 
    int e_global = lelem_ptr->get_elem_loc(ee);
    anaprocId->InsertNextValue( epart_map[e_global] );

    // update offset
    ptOffset += 4;
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
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object ... \n");
  gridData->Delete();
}


void VTK_Writer_Solids_Tet4::interpolateJ(
    const int &ptoffset, const double * const &inputData,
    const FEAElement * const &elem, vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  double * u = new double [nLocBas];
  double * v = new double [nLocBas];
  double * w = new double [nLocBas];

  std::vector<double> ux, uy, uz, vx, vy, vz, wx, wy, wz;

  for(int ii=0; ii<nLocBas; ++ii)
  {
    u[ii] = inputData[ii*3];
    v[ii] = inputData[ii*3+1];
    w[ii] = inputData[ii*3+2];
  }

  intep.interpolateFE_Grad(u, elem, ux, uy, uz);
  intep.interpolateFE_Grad(v, elem, vx, vy, vz);
  intep.interpolateFE_Grad(w, elem, wx, wy, wz);

  Matrix_3x3 F;
  double detF;

  for(int ii=0; ii<nqp; ++ii)
  {
    F(0) = ux[ii] + 1.0; F(1) = uy[ii];       F(2) = uz[ii];
    F(3) = vx[ii];       F(4) = vy[ii] + 1.0; F(5) = vz[ii];
    F(6) = wx[ii];       F(7) = wy[ii];       F(8) = wz[ii] + 1.0;

    detF = F.det();

    vtkData->InsertComponent(ptoffset+ii, 0, detF);
  }

  delete [] u; delete [] v; delete [] w;
}


void VTK_Writer_Solids_Tet4::interpolateJele( const int &ptoffset,
    const double * const &inputData,
    FEAElement * const &elem, vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  double * u = new double [nLocBas];
  double * v = new double [nLocBas];
  double * w = new double [nLocBas];

  std::vector<double> ux, uy, uz, vx, vy, vz, wx, wy, wz;

  for(int ii=0; ii<nLocBas; ++ii)
  {
    u[ii] = inputData[ii*3];
    v[ii] = inputData[ii*3+1];
    w[ii] = inputData[ii*3+2];
  }

  intep.interpolateFE_Grad(u, elem, ux, uy, uz);
  intep.interpolateFE_Grad(v, elem, vx, vy, vz);
  intep.interpolateFE_Grad(w, elem, wx, wy, wz);

  Matrix_3x3 F;
  std::vector<double> detF;

  for(int ii=0; ii<nqp; ++ii)
  {
    F(0) = ux[ii] + 1.0; F(1) = uy[ii];       F(2) = uz[ii];
    F(3) = vx[ii];       F(4) = vy[ii] + 1.0; F(5) = vz[ii];
    F(6) = wx[ii];       F(7) = wy[ii];       F(8) = wz[ii] + 1.0;

    detF.push_back( F.det() );
  }

  double Jave = 0.0;
  for(int ii=0; ii<nqp; ++ii) Jave += detF[ii];

  Jave = Jave / double(nqp);

  for(int ii=0; ii<nqp; ++ii)
    vtkData->InsertComponent(ptoffset+ii, 0, Jave);

  delete [] u; delete [] v; delete [] w;
}


void VTK_Writer_Solids_Tet4::interpolateCauchyZ( const int &ptoffset,
    const double * const &inputData, FEAElement * const &elem,
    IMaterialModel * const &model, vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  double * u = new double [nLocBas];
  double * v = new double [nLocBas];
  double * w = new double [nLocBas];

  std::vector<double> ux, uy, uz, vx, vy, vz, wx, wy, wz;

  for(int ii=0; ii<nLocBas; ++ii)
  {
    u[ii] = inputData[ii*3];
    v[ii] = inputData[ii*3+1];
    w[ii] = inputData[ii*3+2];
  }

  intep.interpolateFE_Grad(u, elem, ux, uy, uz);
  intep.interpolateFE_Grad(v, elem, vx, vy, vz);
  intep.interpolateFE_Grad(w, elem, wx, wy, wz);

  Matrix_3x3 F, sigma;

  for(int ii=0; ii<nqp; ++ii)
  {
    F(0) = ux[ii] + 1.0; F(1) = uy[ii];       F(2) = uz[ii];
    F(3) = vx[ii];       F(4) = vy[ii] + 1.0; F(5) = vz[ii];
    F(6) = wx[ii];       F(7) = wy[ii];       F(8) = wz[ii] + 1.0;

    model->get_Cauchy_stress(F, sigma);

    vtkData->InsertComponent(ptoffset+ii, 0, sigma(0,2));
    vtkData->InsertComponent(ptoffset+ii, 1, sigma(1,2));
    vtkData->InsertComponent(ptoffset+ii, 2, sigma(2,2));
  }

  delete [] u; delete [] v; delete [] w;
}


void VTK_Writer_Solids_Tet4::interpolateCauchyZpPres( const int &ptoffset,
    const double * const &inputData_u,
    const double * const &inputData_p,
    FEAElement * const &elem,
    IMaterialModel * const &model,
    vtkDoubleArray * const &vtkData )
{
  const int nqp = elem -> get_numQuapts();

  double * u = new double [nLocBas];
  double * v = new double [nLocBas];
  double * w = new double [nLocBas];

  std::vector<double> ux, uy, uz, vx, vy, vz, wx, wy, wz;

  for(int ii=0; ii<nLocBas; ++ii)
  {
    u[ii] = inputData_u[ii*3];
    v[ii] = inputData_u[ii*3+1];
    w[ii] = inputData_u[ii*3+2];
  }

  intep.interpolateFE_Grad(u, elem, ux, uy, uz);
  intep.interpolateFE_Grad(v, elem, vx, vy, vz);
  intep.interpolateFE_Grad(w, elem, wx, wy, wz);

  double * p = new double [nLocBas];
  
  for(int ii=0; ii<nLocBas; ++ii) p[ii] = inputData_p[ii];

  std::vector<double> pp;
  intep.interpolateFE(p, elem, pp);

  Matrix_3x3 F, sigma;

  for(int ii=0; ii<nqp; ++ii)
  {
    F(0) = ux[ii] + 1.0; F(1) = uy[ii];       F(2) = uz[ii];
    F(3) = vx[ii];       F(4) = vy[ii] + 1.0; F(5) = vz[ii];
    F(6) = wx[ii];       F(7) = wy[ii];       F(8) = wz[ii] + 1.0;

    model->get_Cauchy_stress(F, sigma);

    vtkData->InsertComponent(ptoffset+ii, 0, sigma(0,2));
    vtkData->InsertComponent(ptoffset+ii, 1, sigma(1,2));
    
    // At the last component slot, I put the scaled Cauchy stress. 
    double val = sigma(0,0) - pp[ii];

    if(val > 1.5e7) val = 1.5e7;
    else if(val < 5.0e6) val = 5.0e6;
    
    vtkData->InsertComponent(ptoffset+ii, 2, val);
  }

  delete [] u; delete [] v; delete [] w; delete [] p;
}



void VTK_Writer_Solids_Tet4::interpolateOrientation( const int &ptoffset,
    const double * const &inputData, FEAElement * const &elem,
    IMaterialModel * const &model,
    vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  double * u = new double [nLocBas];
  double * v = new double [nLocBas];
  double * w = new double [nLocBas];

  std::vector<double> ux, uy, uz, vx, vy, vz, wx, wy, wz;

  Vector_3 a, b;

  model -> get_fibre_dir(0, a(0), a(1), a(2));
  model -> get_fibre_dir(1, b(0), b(1), b(2));

  for(int ii=0; ii<nLocBas; ++ii)
  {
    u[ii] = inputData[ii*3];
    v[ii] = inputData[ii*3+1];
    w[ii] = inputData[ii*3+2];
  }

  intep.interpolateFE_Grad(u, elem, ux, uy, uz);
  intep.interpolateFE_Grad(v, elem, vx, vy, vz);
  intep.interpolateFE_Grad(w, elem, wx, wy, wz);

  Matrix_3x3 F, C;

  for(int ii=0; ii<nqp; ++ii)
  {
    F(0) = ux[ii] + 1.0; F(1) = uy[ii];       F(2) = uz[ii];
    F(3) = vx[ii];       F(4) = vy[ii] + 1.0; F(5) = vz[ii];
    F(6) = wx[ii];       F(7) = wy[ii];       F(8) = wz[ii] + 1.0;

    C.MatMultTransposeLeft(F);

    double ori_measure = C.VecMatVec(a,b);

    double len_a = C.VecMatVec(a,a);
    double len_b = C.VecMatVec(b,b);

    vtkData->InsertComponent(ptoffset+ii, 0, ori_measure/(len_a * len_b));
  }

  delete [] u; delete [] v; delete [] w;
}


void VTK_Writer_Solids_Tet4::interpolateVonMise( const int &ptoffset,
    const double * const &inputData, FEAElement * const &elem,
    IMaterialModel * const &model,
    vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  double * u = new double [nLocBas];
  double * v = new double [nLocBas];
  double * w = new double [nLocBas];

  std::vector<double> ux, uy, uz, vx, vy, vz, wx, wy, wz;

  for(int ii=0; ii<nLocBas; ++ii)
  {
    u[ii] = inputData[ii*3];
    v[ii] = inputData[ii*3+1];
    w[ii] = inputData[ii*3+2];
  }

  intep.interpolateFE_Grad(u, elem, ux, uy, uz);
  intep.interpolateFE_Grad(v, elem, vx, vy, vz);
  intep.interpolateFE_Grad(w, elem, wx, wy, wz);

  Matrix_3x3 F, sigma;

  double vonMises, a, b, c, d;

  for(int ii=0; ii<nqp; ++ii)
  {
    F(0) = ux[ii] + 1.0; F(1) = uy[ii];       F(2) = uz[ii];
    F(3) = vx[ii];       F(4) = vy[ii] + 1.0; F(5) = vz[ii];
    F(6) = wx[ii];       F(7) = wy[ii];       F(8) = wz[ii] + 1.0;

    model->get_Cauchy_stress(F, sigma);

    a = (sigma(0) - sigma(4)) * (sigma(0) - sigma(4));
    b = (sigma(4) - sigma(8)) * (sigma(4) - sigma(8));
    c = (sigma(0) - sigma(8)) * (sigma(0) - sigma(8));
    d = 6.0 * ( sigma(1) * sigma(1) + sigma(2) * sigma(2) + sigma(5) * sigma(5) );
    vonMises = sqrt( 0.5*(a+b+c+d) );

    vtkData->InsertComponent(ptoffset+ii, 0, vonMises);
  }

  delete [] u; delete [] v; delete [] w;
}


// EOF
