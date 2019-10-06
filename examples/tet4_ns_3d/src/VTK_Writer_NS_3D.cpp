#include "VTK_Writer_NS_3D.hpp"

VTK_Writer_NS_3D::VTK_Writer_NS_3D( const int &in_nelem,
    const std::string &epart_file )
: nLocBas(4), nElem(in_nelem), intep(4, true)
{
  VIS_T::read_epart( epart_file, nElem, epart_map );
}


VTK_Writer_NS_3D::~VTK_Writer_NS_3D()
{
  VEC_T::clean( epart_map );
}


void VTK_Writer_NS_3D::writeOutput(
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
  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();
  vtkPoints * points = vtkPoints::New();
  const int numDArrays = vdata_ptr->get_arrayCompSize();
  if(numDArrays != 2) SYS_T::print_fatal("Error: vdata size numDArrays != 2. \n");

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

  int ptOffset = 0;
  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    lien_ptr -> get_LIEN_e(ee, IEN_e);
    fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);
    elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

    intep.interpolateVTKPts(ptOffset, ectrl_x, ectrl_y, ectrl_z,
        elemptr, points );

    // Interpolate the pressure scalar
    std::vector<double> inputInfo; inputInfo.clear();
    int asize = vdata_ptr->get_arraySizes( 0 );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[0][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[0] );

    // Interpolate velo 
    asize = vdata_ptr -> get_arraySizes( 1 );
    inputInfo.clear();
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[1][pt_index * asize + kk] );
    }
    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
        elemptr, dataVecs[1] );

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


void VTK_Writer_NS_3D::writeOutput_compact(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_ptr,
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
  // Make sure nqp == nLocBas == 4
  SYS_T::print_fatal_if(quad->get_num_quadPts() != 4, "Error: VTK_Writer_Tet4 requires 4 quadrature points.\n");

  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();
  vtkPoints * points = vtkPoints::New();
  const int numDArrays = vdata_ptr->get_arrayCompSize();
  if(numDArrays != 2) SYS_T::print_fatal("Error: vdata size numDArrays != 2.\n");

  vtkDoubleArray ** dataVecs = new vtkDoubleArray * [numDArrays];
  for(int ii=0; ii<numDArrays; ++ii)
  {
    dataVecs[ii] = vtkDoubleArray::New();
    dataVecs[ii] -> SetNumberOfComponents( vdata_ptr->get_arraySizes(ii) );
    std::string temp_name = vdata_ptr->get_arrayNames(ii);
    dataVecs[ii] -> SetName( temp_name.c_str() );
    dataVecs[ii] -> SetNumberOfTuples( num_of_nodes );
  }

  vtkIntArray * anaprocId = vtkIntArray::New();
  anaprocId -> SetName("Analysis_Partition");
  anaprocId -> SetNumberOfComponents(1);

  // allocate holder for grad velo
  std::vector<double> ux, uy, uz, vx, vy, vz, wx, wy, wz, num_adj_cell;
  ux.resize( num_of_nodes ); uy.resize( num_of_nodes ); uz.resize( num_of_nodes );
  vx.resize( num_of_nodes ); vy.resize( num_of_nodes ); vz.resize( num_of_nodes );
  wx.resize( num_of_nodes ); wy.resize( num_of_nodes ); wz.resize( num_of_nodes );
  num_adj_cell.resize( num_of_nodes );
  
  for(int ii=0; ii<num_of_nodes; ++ii)
  {
    ux[ii] = 0.0; uy[ii] = 0.0; uz[ii] = 0.0;
    vx[ii] = 0.0; vy[ii] = 0.0; vz[ii] = 0.0;
    wx[ii] = 0.0; wy[ii] = 0.0; wz[ii] = 0.0;
    num_adj_cell[ii] = 0.0;
  }

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    lien_ptr -> get_LIEN_e(ee, IEN_e);
    fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);
    elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

    intep.interpolateVTKPts(IEN_e, ectrl_x, ectrl_y, ectrl_z, elemptr, points );

    // Interpolate the pressure scalar
    std::vector<double> inputInfo; inputInfo.clear();
    int asize = vdata_ptr->get_arraySizes( 0 ); // 1
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[0][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, IEN_e, inputInfo, elemptr, dataVecs[0] );

    // Interpolate velo 
    asize = vdata_ptr -> get_arraySizes( 1 ); // 3
    inputInfo.clear();
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[1][pt_index * asize + kk] );
    }
    intep.interpolateVTKData( asize, IEN_e, inputInfo, elemptr, dataVecs[1] );

    // Interpolate velocity gradient as algebraic average of
    // values in adjacent cells for each node.
    std::vector<double> input_u, input_v, input_w;
    input_u.clear(); input_v.clear(); input_w.clear();
    for(int jj=0; jj<nLocBas; ++jj)
    {
      input_u.push_back( pointArrays[1][ IEN_e[jj]*3 + 0 ] );
      input_v.push_back( pointArrays[1][ IEN_e[jj]*3 + 1 ] );
      input_w.push_back( pointArrays[1][ IEN_e[jj]*3 + 2 ] );
    }
    std::vector<double> dudx, dudy, dudz;
    intep.interpolateFE_Grad(input_u, elemptr, dudx, dudy, dudz);

    std::vector<double> dvdx, dvdy, dvdz;
    intep.interpolateFE_Grad(input_v, elemptr, dvdx, dvdy, dvdz);

    std::vector<double> dwdx, dwdy, dwdz;
    intep.interpolateFE_Grad(input_w, elemptr, dwdx, dwdy, dwdz);
    
    for(int ii=0; ii<4; ++ii)
    {
      ux[ IEN_e[ii] ] += dudx[ii];
      uy[ IEN_e[ii] ] += dudy[ii];
      uz[ IEN_e[ii] ] += dudz[ii];

      vx[ IEN_e[ii] ] += dvdx[ii];
      vy[ IEN_e[ii] ] += dvdy[ii];
      vz[ IEN_e[ii] ] += dvdz[ii];

      wx[ IEN_e[ii] ] += dwdx[ii];
      wy[ IEN_e[ii] ] += dwdy[ii];
      wz[ IEN_e[ii] ] += dwdz[ii];

      num_adj_cell[ IEN_e[ii] ] += 1.0;
    }

    // Set mesh connectivity
    VIS_T::setTetraelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3], gridData );

    // Analysis mesh partition 
    const int e_global = lelem_ptr->get_elem_loc(ee);
    anaprocId->InsertNextValue( epart_map[e_global] );
  }

  gridData -> SetPoints( points );
  points -> Delete();

  for(int ii=0; ii<numDArrays; ++ii)
  {
    gridData->GetPointData()->AddArray( dataVecs[ii] );
    dataVecs[ii]->Delete();
  }
  
  delete [] dataVecs;

  // Do notal average and insert grad v data
  vtkDoubleArray * vorticity = vtkDoubleArray::New();

  vorticity -> SetName("Velocity Gradient");
  vorticity -> SetNumberOfComponents( 9 );
  vorticity -> SetNumberOfTuples( num_of_nodes );

  for(int ii=0; ii<num_of_nodes; ++ii)
  {
    vorticity -> InsertComponent(ii, 0, ux[ii] / num_adj_cell[ii] );
    vorticity -> InsertComponent(ii, 1, uy[ii] / num_adj_cell[ii] );
    vorticity -> InsertComponent(ii, 2, uz[ii] / num_adj_cell[ii] );

    vorticity -> InsertComponent(ii, 3, vx[ii] / num_adj_cell[ii] );
    vorticity -> InsertComponent(ii, 4, vy[ii] / num_adj_cell[ii] );
    vorticity -> InsertComponent(ii, 5, vz[ii] / num_adj_cell[ii] );

    vorticity -> InsertComponent(ii, 6, wx[ii] / num_adj_cell[ii] );
    vorticity -> InsertComponent(ii, 7, wy[ii] / num_adj_cell[ii] );
    vorticity -> InsertComponent(ii, 8, wz[ii] / num_adj_cell[ii] );
  }

  gridData->GetPointData()->AddArray( vorticity );
  vorticity -> Delete();

  // Add cell data
  gridData->GetCellData()->AddArray(anaprocId);
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

// EOF
