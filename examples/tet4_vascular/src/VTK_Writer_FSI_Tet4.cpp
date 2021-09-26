#include "VTK_Writer_FSI_Tet4.hpp"

VTK_Writer_FSI_Tet4::VTK_Writer_FSI_Tet4( const int &in_nelem,
    const std::string &epart_file )
: nLocBas(4), nElem(in_nelem)
{
  VIS_T::read_epart( epart_file, nElem, epart_map );
}


VTK_Writer_FSI_Tet4::~VTK_Writer_FSI_Tet4()
{
  VEC_T::clean(epart_map);
}


void VTK_Writer_FSI_Tet4::writeOutput(
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
    lien_ptr -> get_LIEN(ee, IEN_e);
    fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);
    elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

    // Interpolate data and assign to dataVecs
    std::vector<double> inputInfo; inputInfo.clear();
    int asize = vdata_ptr -> get_arraySizes(0);

    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[0][pt_index * asize + kk] );
    }

    // displacement interpolation
    intep.interpolateVTKData( asize, IEN_e, inputInfo, elemptr, dataVecs[0] );

    // use displacement to update points
    intep.interpolateVTKPts(IEN_e, ectrl_x, ectrl_y, ectrl_z,
        inputInfo, elemptr, points); 

    // Interpolate detF
    interpolateJ( IEN_e, inputInfo, elemptr, dataVecs[1] );

    // Interpolate the pressure scalar
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes(2);
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, IEN_e, inputInfo, elemptr, dataVecs[2] );

    // Interpolate the velocity vector
    inputInfo.clear();
    asize = vdata_ptr->get_arraySizes(3);
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, IEN_e, inputInfo, elemptr, dataVecs[3] );

    // Set mesh connectivity
    VIS_T::setTetraelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3], gridData );

    // Analysis mesh partition 
    int e_global = lelem_ptr->get_elem_loc(ee);
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
  VIS_T::writeVisFile( gridData, outputBName, outputName,
      rank, size, sol_time, isXML );

  // Clean gridData
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object.\n");
  gridData->Delete();
}


void VTK_Writer_FSI_Tet4::writeOutput_fluid(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_ptr,
    const std::vector<int> &fien,
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
  if(numDArrays != 3) SYS_T::print_fatal("Error: vdata size numDArrays != 3. \n");

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

  // Holder for grad v and number of adjacent elements for node points
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

  int IEN_f[4];

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    if( lelem_ptr->get_elem_tag(ee) == 0 )
    {
      lien_ptr -> get_LIEN(ee, IEN_e);

      fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

      IEN_f[0] = fien[ee*4+0]; IEN_f[1] = fien[ee*4+1]; 
      IEN_f[2] = fien[ee*4+2]; IEN_f[3] = fien[ee*4+3]; 

      // Interpolate data and assign to dataVecs
      std::vector<double> inputInfo; inputInfo.clear();
      int asize = vdata_ptr -> get_arraySizes(0);

      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_e[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[0][pt_index * asize + kk] );
      }

      // displacement interpolation
      intep.interpolateVTKData( asize, IEN_f, inputInfo, elemptr, dataVecs[0] );

      // use displacement to update points
      intep.interpolateVTKPts( IEN_f, ectrl_x, ectrl_y, ectrl_z,
          inputInfo, elemptr, points ); 

      // Interpolate the pressure scalar
      inputInfo.clear();
      asize = vdata_ptr->get_arraySizes(1);
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_e[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, IEN_f, inputInfo, elemptr, dataVecs[1] );

      // Interpolate the velocity vector
      inputInfo.clear();
      asize = vdata_ptr->get_arraySizes( 2 );
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_e[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, IEN_f, inputInfo, elemptr, dataVecs[2] );

      // Interpolate for velocity gradient
      std::vector<double> input_u, input_v, input_w;
      input_u.clear(); input_v.clear(); input_w.clear();
      for(int jj=0; jj<nLocBas; ++jj)
      {
        input_u.push_back( pointArrays[2][ IEN_e[jj]*3 + 0 ] );
        input_v.push_back( pointArrays[2][ IEN_e[jj]*3 + 1 ] );
        input_w.push_back( pointArrays[2][ IEN_e[jj]*3 + 2 ] );
      }
      
      std::vector<double> dudx, dudy, dudz;
      intep.interpolateFE_Grad(input_u, elemptr, dudx, dudy, dudz);

      std::vector<double> dvdx, dvdy, dvdz;
      intep.interpolateFE_Grad(input_v, elemptr, dvdx, dvdy, dvdz);

      std::vector<double> dwdx, dwdy, dwdz;
      intep.interpolateFE_Grad(input_w, elemptr, dwdx, dwdy, dwdz);

      for(int ii=0; ii<nLocBas; ++ii)
      {
        ux[ IEN_f[ii] ] += dudx[ii]; 
        uy[ IEN_f[ii] ] += dudy[ii];
        uz[ IEN_f[ii] ] += dudz[ii];
        
        vx[ IEN_f[ii] ] += dvdx[ii]; 
        vy[ IEN_f[ii] ] += dvdy[ii];
        vz[ IEN_f[ii] ] += dvdz[ii];
        
        wx[ IEN_f[ii] ] += dwdx[ii]; 
        wy[ IEN_f[ii] ] += dwdy[ii];
        wz[ IEN_f[ii] ] += dwdz[ii];
        
        num_adj_cell[ IEN_f[ii] ] += 1.0;
      }

      // Set mesh connectivity
      VIS_T::setTetraelem( IEN_f[0], IEN_f[1], IEN_f[2], IEN_f[3], gridData );

      // Analysis mesh partition 
      int e_global = lelem_ptr->get_elem_loc(ee);
      anaprocId->InsertNextValue( epart_map[e_global] );
    }
  }

  gridData -> SetPoints( points ); points -> Delete();

  for(int ii=0; ii<numDArrays; ++ii)
  {
    gridData->GetPointData()->AddArray( dataVecs[ii] );
    dataVecs[ii]->Delete();
  }
  delete [] dataVecs;

  // Insert Grad v data
  vtkDoubleArray * vorticity = vtkDoubleArray::New();

  vorticity -> SetNumberOfComponents( 9 );
  vorticity -> SetName("Velocity Gradient");
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
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object.\n");
  gridData->Delete();
}


void VTK_Writer_FSI_Tet4::writeOutput_solid(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_ptr,
    const std::vector<int> &fien,
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

  int IEN_f[4];

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    if( lelem_ptr->get_elem_tag(ee) == 1 )
    {
      lien_ptr -> get_LIEN(ee, IEN_e);

      fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

      IEN_f[0] = fien[ee*4+0]; IEN_f[1] = fien[ee*4+1];
      IEN_f[2] = fien[ee*4+2]; IEN_f[3] = fien[ee*4+3];

      // Interpolate data and assign to dataVecs
      std::vector<double> inputInfo; inputInfo.clear();
      int asize = vdata_ptr -> get_arraySizes(0); // 3

      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_e[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[0][pt_index * asize + kk] );
      }

      // displacement interpolation
      intep.interpolateVTKData( asize, IEN_f, inputInfo, elemptr, dataVecs[0] );

      // use displacement to update points
      intep.interpolateVTKPts( IEN_f, ectrl_x, ectrl_y, ectrl_z,
          inputInfo, elemptr, points ); 

      // Interpolate detF
      interpolateJ( IEN_f, inputInfo, elemptr, dataVecs[1] );

      // Interpolate the pressure scalar
      inputInfo.clear();
      asize = vdata_ptr->get_arraySizes(2); // 1
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_e[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, IEN_f, inputInfo, elemptr, dataVecs[2] );

      // Interpolate the velocity vector
      inputInfo.clear();
      asize = vdata_ptr->get_arraySizes(3); // 3
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_e[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, IEN_f, inputInfo, elemptr, dataVecs[3] );

      // Set mesh connectivity
      VIS_T::setTetraelem( IEN_f[0], IEN_f[1], IEN_f[2], IEN_f[3], gridData );

      // Analysis mesh partition 
      int e_global = lelem_ptr->get_elem_loc(ee);
      anaprocId->InsertNextValue( epart_map[e_global] );
    }
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


void VTK_Writer_FSI_Tet4::writeOutput_solid_ref(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_ptr,
    const std::vector<int> &fien,
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
  if(numDArrays != 3) SYS_T::print_fatal("Error: vdata size numDArrays != 3. \n");

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

  int IEN_f[4];

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    if( lelem_ptr->get_elem_tag(ee) == 1 )
    {
      lien_ptr -> get_LIEN(ee, IEN_e);

      fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

      IEN_f[0] = fien[ee*4+0]; IEN_f[1] = fien[ee*4+1];
      IEN_f[2] = fien[ee*4+2]; IEN_f[3] = fien[ee*4+3];

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
      intep.interpolateVTKData( asize, IEN_f, &inputInfo[0],
          elemptr, dataVecs[0] );

      // interpolate points
      intep.interpolateVTKPts(IEN_f, ectrl_x, ectrl_y, ectrl_z,
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

      intep.interpolateVTKData( asize, IEN_f, &inputInfo[0],
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
      intep.interpolateVTKData( asize, IEN_f, &inputInfo[0],
          elemptr, dataVecs[2] );

      // Set mesh connectivity
      VIS_T::setTetraelem( IEN_f[0], IEN_f[1], IEN_f[2], IEN_f[3], gridData );

      // Analysis mesh partition 
      int e_global = lelem_ptr->get_elem_loc(ee);
      anaprocId->InsertNextValue( epart_map[e_global] );
    }
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



void VTK_Writer_FSI_Tet4::interpolateJ(
    const int * const &ptid, const std::vector<double> &inputData,
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
  
  Interpolater intep( nLocBas );

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

    vtkData->InsertComponent(ptid[ii], 0, detF);
  }

  delete [] u; delete [] v; delete [] w;
}

// EOF
