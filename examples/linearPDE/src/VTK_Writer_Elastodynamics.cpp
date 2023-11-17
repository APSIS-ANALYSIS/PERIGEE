#include "VTK_Writer_Elastodynamics.hpp"

VTK_Writer_Elastodynamics::VTK_Writer_Elastodynamics( const int &in_nelem,
    const int &in_nlocbas, const std::string &epart_file,
    const double &in_module_E, const double &in_nu )
: nLocBas( in_nlocbas ), nElem( in_nelem ),
  lambda( in_nu * in_module_E / ((1.0 + in_nu) * (1.0 - 2.0 * in_nu)) ),
  mu( 0.5 * in_module_E / (1.0 + in_nu) )
{
  VIS_T::read_epart( epart_file, nElem, epart_map );
}

void VTK_Writer_Elastodynamics::writeOutput(
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

  // Allocate VTK gridData object
  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();

  // Allocate vtk point
  vtkPoints * points = vtkPoints::New();

  // Number of arrays visualized
  const int numDArrays = vdata_ptr->get_arrayCompSize();

  vtkDoubleArray ** dataVecs = new vtkDoubleArray * [numDArrays];
  for(int ii=0; ii<numDArrays; ++ii)
  {
    dataVecs[ii] = vtkDoubleArray::New();
    dataVecs[ii] -> SetNumberOfComponents( vdata_ptr->get_arraySizes(ii) );
    std::string temp_name = vdata_ptr->get_arrayNames(ii);
    dataVecs[ii] -> SetName( temp_name.c_str() );
  }

  // An integer array for the analysis mesh partition
  vtkIntArray * anaprocId = vtkIntArray::New();
  anaprocId -> SetName("Analysis_Partition");
  anaprocId -> SetNumberOfComponents(1);

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    const std::vector<int> IEN_e = lien_ptr -> get_LIEN( ee );

    std::vector<double> ectrl_x ( nLocBas, 0.0 );
    std::vector<double> ectrl_y ( nLocBas, 0.0 );
    std::vector<double> ectrl_z ( nLocBas, 0.0 );

    fnode_ptr -> get_ctrlPts_xyz(nLocBas, &IEN_e[0], &ectrl_x[0], &ectrl_y[0], &ectrl_z[0]);

    elemptr->buildBasis( quad, &ectrl_x[0], &ectrl_y[0], &ectrl_z[0] );
    
    // Interpolate nodal coordinates
    intep.interpolateVTKPts(&IEN_e[0], &ectrl_x[0], &ectrl_y[0], &ectrl_z[0],
        elemptr, points );
  
    std::vector<double> inputInfo; inputInfo.clear();
    
    // Interpolate velocity vector
    int asize = vdata_ptr->get_arraySizes(0);
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[0][pt_index * asize + kk ] );
    }
    intep.interpolateVTKData( asize, &IEN_e[0], &inputInfo[0],
        elemptr, dataVecs[0] );

    interpolateCauchy( &IEN_e[0], &inputInfo[0], elemptr, dataVecs[1] );

    // Set mesh connectivity
    if( elemptr->get_Type() == 501 )
      VIS_T::setTetraelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3], gridData );
    else if( elemptr->get_Type() == 502 )
      VIS_T::setQuadTetraelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3], 
          IEN_e[4], IEN_e[5], IEN_e[6], IEN_e[7], IEN_e[8], IEN_e[9],
          gridData );
    else if( elemptr->get_Type() == 601 )
      VIS_T::setHexelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3], 
        IEN_e[4], IEN_e[5], IEN_e[6], IEN_e[7], gridData );
    else if( elemptr->get_Type() == 602 )
      VIS_T::setTriQuadHexelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3], 
        IEN_e[4], IEN_e[5], IEN_e[6], IEN_e[7], IEN_e[8], IEN_e[9],
        IEN_e[10], IEN_e[11], IEN_e[12], IEN_e[13], IEN_e[14], IEN_e[15],
        IEN_e[16], IEN_e[17], IEN_e[18], IEN_e[19], IEN_e[20], IEN_e[21],
        IEN_e[22], IEN_e[23], IEN_e[24], IEN_e[25], IEN_e[26], gridData );
    else SYS_T::print_fatal("Error: unknown element type.\n");

    // Mesh partition info
    anaprocId->InsertNextValue( epart_map[ lelem_ptr->get_elem_loc(ee) ] );
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

  // Write the grid Data
  VIS_T::writeVisFile( gridData, outputBName, outputName,
      rank, size, sol_time, isXML );

  // Clean gridData
  SYS_T::commPrint("-- Clean gridData object.\n");
  gridData->Delete();
}

void VTK_Writer_Elastodynamics::interpolateCauchy( const int * const &ptid,
    const double * const &inputData, FEAElement * const &elem,
    vtkDoubleArray * const &vtkData)
{
  Interpolater intep( nLocBas );

  const double l2mu = lambda + 2.0 * mu;

  const int nqp = elem->get_numQuapts();

  double * ux = new double [nLocBas];
  double * uy = new double [nLocBas];
  double * uz = new double [nLocBas];

  std::vector<double> ux_x, ux_y, ux_z, uy_x, uy_y, uy_z, uz_x, uz_y, uz_z;

  for(int ii=0; ii<nLocBas; ++ii)
  {
    ux[ii] = inputData[ii*3];
    uy[ii] = inputData[ii*3+1];
    uz[ii] = inputData[ii*3+2];
  }

  intep.interpolateFE_Grad(ux, elem, ux_x, ux_y, ux_z);
  intep.interpolateFE_Grad(uy, elem, uy_x, uy_y, uy_z);
  intep.interpolateFE_Grad(uz, elem, uz_x, uz_y, uz_z);

  double sigma_xx = 0.0, sigma_yy = 0.0, sigma_zz = 0.0;
  double sigma_xy = 0.0, sigma_xz = 0.0, sigma_yz = 0.0;

  for(int ii=0; ii<nqp; ++ii)
  {
    sigma_xx = l2mu * ux_x[ii] + mu * (uy_y[ii] + uz_z[ii]);
    sigma_yy = l2mu * uy_y[ii] + mu * (ux_x[ii] + uz_z[ii]);
    sigma_zz = l2mu * uz_z[ii] + mu * (uy_y[ii] + ux_x[ii]);
    sigma_xy = mu * ( ux_y[ii] + uy_x[ii] );
    sigma_xz = mu * ( ux_z[ii] + uz_x[ii] );
    sigma_yz = mu * ( uy_z[ii] + uz_y[ii] );
    vtkData->InsertComponent(ptid[ii], 0, sigma_xx);
    vtkData->InsertComponent(ptid[ii], 1, sigma_yy);
    vtkData->InsertComponent(ptid[ii], 2, sigma_zz);
    vtkData->InsertComponent(ptid[ii], 3, sigma_xy);
    vtkData->InsertComponent(ptid[ii], 4, sigma_yz);
    vtkData->InsertComponent(ptid[ii], 5, sigma_xz);
  }

  delete [] ux; delete [] uy; delete [] uz;
}

// EOF