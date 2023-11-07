#include "VTK_Writer_FSI.hpp"

VTK_Writer_FSI::VTK_Writer_FSI( const int &in_nelem,
    const int &in_nlocbas, const std::string &epart_file ) : nLocBas(in_nlocbas), nElem(in_nelem)
{
  VIS_T::read_epart( epart_file, nElem, epart_map );
}


VTK_Writer_FSI::~VTK_Writer_FSI()
{
  VEC_T::clean(epart_map);
}

void VTK_Writer_FSI::interpolateJ( const int * const &ptid,
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
    Tensor2_3D F( ux[ii] + 1.0, uy[ii],       uz[ii],
                  vx[ii],       vy[ii] + 1.0, vz[ii],
                  wx[ii],       wy[ii],       wz[ii] + 1.0 );

    vtkData->InsertComponent( ptid[ii], 0, F.det() );
  }
}

void VTK_Writer_FSI::writeOutput(
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
        const std::string &outputBName,
        const std::string &outputName,
        const bool &isXML )
{
  if(nLocBas == 4)  // elemType 501
  {
    // This routine requires nqp = 4
    SYS_T::print_fatal_if(quad->get_num_quadPts() != 4, "Error: VTK_Writer requires 4 quadrature points for Tet4.\n");
  }
  else if(nLocBas == 8)  // elemType 601
  {
    // This routine requires nqp = 8
    SYS_T::print_fatal_if(quad->get_num_quadPts() != 8, "Error: VTK_Writer requires 4 quadrature points for Hex8.\n");    
  }
  else SYS_T::print_fatal( "Error: VTK_Writer_FSI::writeOutput function: unsupported element type \n" );

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

    double ectrl_x[nLocBas], ectrl_y[nLocBas], ectrl_z[nLocBas];

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
    if( elemptr->get_Type() == 501 )
      VIS_T::setTetraelem( IEN_p[0], IEN_p[1], IEN_p[2], IEN_p[3], gridData );
    else if( elemptr->get_Type() == 601 )
      VIS_T::setHexelem( IEN_p[0], IEN_p[1], IEN_p[2], IEN_p[3], 
        IEN_p[4], IEN_p[5], IEN_p[6], IEN_p[7], gridData );
    else SYS_T::print_fatal("Error: unknown element type.\n");

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
  SYS_T::commPrint("-- Clean gridData object.\n");
  gridData->Delete();
}

void VTK_Writer_FSI::writeOutput_fluid(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const std::vector<int> &fien,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    FEAElement * const &elemptr,
    const IQuadPts * const &quad,
    const double * const * const &pointArrays,
    const int &rank, const int &size,
    const int &num_of_nodes,
    const double &sol_time,
    const std::string &outputBName,
    const std::string &outputName,
    const bool &isXML )
{
  if(nLocBas == 4)  // elemType 501
  {
    // This routine requires nqp = 4
    SYS_T::print_fatal_if(quad->get_num_quadPts() != 4, "Error: VTK_Writer requires 4 quadrature points for Tet4.\n");
  }
  else if(nLocBas == 8)  // elemType 601
  {
    // This routine requires nqp = 8
    SYS_T::print_fatal_if(quad->get_num_quadPts() != 8, "Error: VTK_Writer requires 4 quadrature points for Hex8.\n");    
  }
  else SYS_T::print_fatal( "Error: VTK_Writer_FSI::writeOutput_fluid function: unsupported element type \n" );

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
    dataVecs[ii] -> SetName( vdata_ptr->get_arrayNames(ii).c_str() );
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

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    if( lelem_ptr->get_elem_tag(ee) == 0 )
    {
      const std::vector<int> IEN_v = lien_v -> get_LIEN(ee);
      const std::vector<int> IEN_p = lien_p -> get_LIEN(ee);

      double ectrl_x[nLocBas], ectrl_y[nLocBas], ectrl_z[nLocBas];

      fnode_ptr -> get_ctrlPts_xyz(nLocBas, &IEN_v[0], ectrl_x, ectrl_y, ectrl_z);

      elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

      int IEN_f[nLocBas]{};

      for(int ii=0; ii<nLocBas; ++ii) IEN_f[ii] = fien[ee * nLocBas + ii]; 

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
      intep.interpolateVTKData( asize, IEN_f, inputInfo, elemptr, dataVecs[0] );

      // use displacement to update points
      intep.interpolateVTKPts( IEN_f, ectrl_x, ectrl_y, ectrl_z, inputInfo, elemptr, points );

      // Interpolate the pressure scalar
      inputInfo.clear();
      asize = vdata_ptr->get_arraySizes(1);
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_p[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, IEN_f, inputInfo, elemptr, dataVecs[1] );

      // Interpolate the velocity vector
      inputInfo.clear();
      asize = vdata_ptr->get_arraySizes( 2 );
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_v[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, IEN_f, inputInfo, elemptr, dataVecs[2] );
      
      // Interpolate for velocity gradient
      std::vector<double> input_u, input_v, input_w;
      input_u.clear(); input_v.clear(); input_w.clear();
      for(int jj=0; jj<nLocBas; ++jj)
      {
        input_u.push_back( pointArrays[2][ IEN_v[jj]*3 + 0 ] );
        input_v.push_back( pointArrays[2][ IEN_v[jj]*3 + 1 ] );
        input_w.push_back( pointArrays[2][ IEN_v[jj]*3 + 2 ] );
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
      if( elemptr->get_Type() == 501 )
        VIS_T::setTetraelem( IEN_f[0], IEN_f[1], IEN_f[2], IEN_f[3], gridData );
      else if( elemptr->get_Type() == 601 )
        VIS_T::setHexelem( IEN_f[0], IEN_f[1], IEN_f[2], IEN_f[3], 
          IEN_f[4], IEN_f[5], IEN_f[6], IEN_f[7], gridData );
      else SYS_T::print_fatal("Error: unknown element type.\n");

      // Analysis mesh partition
      const int e_global = lelem_ptr->get_elem_loc(ee);
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
  VIS_T::writeVisFile( gridData, outputBName, outputName, rank, size, sol_time, isXML );

  // Clean gridData
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object.\n");
  gridData->Delete();
}


void VTK_Writer_FSI::writeOutput_solid_cur(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const std::vector<int> &sien,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    FEAElement * const &elemptr,
    const IQuadPts * const &quad,
    const double * const * const &pointArrays,
    const int &rank, const int &size,
    const int &num_of_nodes,
    const double &sol_time,
    const std::string &outputBName,
    const std::string &outputName,
    const bool &isXML )
{
  // This routine requires nqp = 4
  SYS_T::print_fatal_if(quad->get_num_quadPts() != 4, "Error: VTK_Writer requires 4 quadrature points for Tet4.\n");

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
  if(numDArrays != 4) SYS_T::print_fatal("Error: vdata size numDArrays != 4. \n");

  // dataVecs is the holder of the interpolated data
  vtkDoubleArray ** dataVecs = new vtkDoubleArray * [numDArrays];
  for(int ii=0; ii<numDArrays; ++ii)
  {
    dataVecs[ii] = vtkDoubleArray::New();
    dataVecs[ii] -> SetNumberOfComponents( vdata_ptr->get_arraySizes(ii) );
    dataVecs[ii] -> SetName( vdata_ptr->get_arrayNames(ii).c_str() );
    dataVecs[ii] -> SetNumberOfTuples( num_of_nodes );
  }

  // An additional array for the analysis partition info
  vtkIntArray * anaprocId = vtkIntArray::New();
  anaprocId -> SetName("Analysis_Partition");
  anaprocId -> SetNumberOfComponents(1);

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    if( lelem_ptr->get_elem_tag(ee) == 1 )
    {
      const std::vector<int> IEN_v = lien_v -> get_LIEN(ee);
      const std::vector<int> IEN_p = lien_p -> get_LIEN(ee);

      double ectrl_x[nLocBas], ectrl_y[nLocBas], ectrl_z[nLocBas];

      fnode_ptr -> get_ctrlPts_xyz(nLocBas, &IEN_v[0], ectrl_x, ectrl_y, ectrl_z);

      elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

      const int IEN_s[4] = { sien[ee*4+0], sien[ee*4+1], sien[ee*4+2], sien[ee*4+3] };

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
      intep.interpolateVTKData( asize, IEN_s, inputInfo, elemptr, dataVecs[0] );

      // use displacement to update points
      intep.interpolateVTKPts( IEN_s, ectrl_x, ectrl_y, ectrl_z, inputInfo, elemptr, points );

      // Interpolate detF
      interpolateJ( IEN_s, inputInfo, elemptr, dataVecs[1] );

      // Interpolate the pressure scalar
      inputInfo.clear();
      asize = vdata_ptr->get_arraySizes(2);  
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_p[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, IEN_s, inputInfo, elemptr, dataVecs[2] ); 
      
      // Interpolate the velocity vector
      inputInfo.clear();
      asize = vdata_ptr->get_arraySizes( 3 );  
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_v[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, IEN_s, inputInfo, elemptr, dataVecs[3] );
      
      // Set mesh connectivity
      VIS_T::setTetraelem( IEN_s[0], IEN_s[1], IEN_s[2], IEN_s[3], gridData );

      // Analysis mesh partition
      const int e_global = lelem_ptr->get_elem_loc(ee);
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
  VIS_T::writeVisFile( gridData, outputBName, outputName, rank, size, sol_time, isXML );

  // Clean gridData
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object.\n");
  gridData->Delete();
}

void VTK_Writer_FSI::writeOutput_solid_ref(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_v,
    const ALocal_IEN * const &lien_p,
    const std::vector<int> &sien,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    FEAElement * const &elemptr,
    const IQuadPts * const &quad,
    const double * const * const &pointArrays,
    const int &rank, const int &size,
    const int &num_of_nodes,
    const double &sol_time,
    const std::string &outputBName,
    const std::string &outputName,
    const bool &isXML )
{
  // This routine requires nqp = 4
  SYS_T::print_fatal_if(quad->get_num_quadPts() != 4, "Error: VTK_Writer requires 4 quadrature points for Tet4.\n");

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
    dataVecs[ii] -> SetName( vdata_ptr->get_arrayNames(ii).c_str() );
    dataVecs[ii] -> SetNumberOfTuples( num_of_nodes );
  }

  // An additional array for the analysis partition info
  vtkIntArray * anaprocId = vtkIntArray::New();
  anaprocId -> SetName("Analysis_Partition");
  anaprocId -> SetNumberOfComponents(1);

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    if( lelem_ptr->get_elem_tag(ee) == 1 )
    {
      const std::vector<int> IEN_v = lien_v -> get_LIEN(ee);
      const std::vector<int> IEN_p = lien_p -> get_LIEN(ee);

      double ectrl_x[nLocBas], ectrl_y[nLocBas], ectrl_z[nLocBas];

      fnode_ptr -> get_ctrlPts_xyz(nLocBas, &IEN_v[0], ectrl_x, ectrl_y, ectrl_z);

      elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

      const int IEN_s[4] = { sien[ee*4+0], sien[ee*4+1], sien[ee*4+2], sien[ee*4+3] };

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
      intep.interpolateVTKData( asize, IEN_s, inputInfo, elemptr, dataVecs[0] );

      // interpolate the coordinates of the points
      intep.interpolateVTKPts( IEN_s, ectrl_x, ectrl_y, ectrl_z, elemptr, points );

      // Interpolate the pressure scalar
      inputInfo.clear();
      asize = vdata_ptr->get_arraySizes(1);  
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_p[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, IEN_s, inputInfo, elemptr, dataVecs[1] ); 
     
      // Interpolate the velocity vector
      inputInfo.clear();
      asize = vdata_ptr->get_arraySizes( 2 );  
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_v[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, IEN_s, inputInfo, elemptr, dataVecs[2] );
      
      // Set mesh connectivity
      VIS_T::setTetraelem( IEN_s[0], IEN_s[1], IEN_s[2], IEN_s[3], gridData );

      // Analysis mesh partition
      const int e_global = lelem_ptr->get_elem_loc(ee);
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
  VIS_T::writeVisFile( gridData, outputBName, outputName, rank, size, sol_time, isXML );

  // Clean gridData
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object.\n");
  gridData->Delete();
}

// EOF
