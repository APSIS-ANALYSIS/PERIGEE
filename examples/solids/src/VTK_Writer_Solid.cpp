#include "VTK_Writer_Solid.hpp"
#include "Sys_Tools.hpp"

VTK_Writer_Solid::VTK_Writer_Solid( const int &in_nelem,
    const int &in_nlocbas, const std::string &epart_file,
    std::unique_ptr<MaterialModel_Mixed_Elasticity> in_matmodel )
: nLocBas( in_nlocbas ), nElem( in_nelem ),
  matmodel(std::move(in_matmodel))
{
  epart_map = VIS_T::read_epart( epart_file, nElem );
}

void VTK_Writer_Solid::interpolateJ_Cauchy( const int * const &ptid,
    const double * const &dispData,
    const double * const &presData,
    const FEAElement * const &elem,
    vtkDoubleArray * const &detFData,
    vtkDoubleArray * const &cauchyData ) const
{
  const int nqp = elem->get_numQuapts();

  std::vector<double> ux(nLocBas, 0.0), uy(nLocBas, 0.0), uz(nLocBas, 0.0);
  std::vector<double> pp(nLocBas, 0.0);

  for(int ii=0; ii<nLocBas; ++ii)
  {
    ux[ii] = dispData[ii*3  ];
    uy[ii] = dispData[ii*3+1];
    uz[ii] = dispData[ii*3+2];
    pp[ii] = presData[ii];
  }

  std::vector<double> ux_x, ux_y, ux_z;
  std::vector<double> uy_x, uy_y, uy_z;
  std::vector<double> uz_x, uz_y, uz_z;
  std::vector<double> p_qp;

  Interp::FE_Grad( &ux[0], elem, ux_x, ux_y, ux_z );
  Interp::FE_Grad( &uy[0], elem, uy_x, uy_y, uy_z );
  Interp::FE_Grad( &uz[0], elem, uz_x, uz_y, uz_z );
  Interp::FE( &pp[0], elem, p_qp );

  for(int ii=0; ii<nqp; ++ii)
  {
    const Tensor2_3D F( ux_x[ii] + 1.0, ux_y[ii],       ux_z[ii],
                        uy_x[ii],       uy_y[ii] + 1.0, uy_z[ii],
                        uz_x[ii],       uz_y[ii],       uz_z[ii] + 1.0 );

    const double detF = F.det();
    SymmTensor2_3D sigma = matmodel->get_Cauchy_stress( F );
    sigma.xx() -= p_qp[ii];
    sigma.yy() -= p_qp[ii];
    sigma.zz() -= p_qp[ii];

    detFData->InsertComponent( ptid[ii], 0, detF );
    cauchyData->InsertComponent( ptid[ii], 0, sigma.xx() );
    cauchyData->InsertComponent( ptid[ii], 1, sigma.yy() );
    cauchyData->InsertComponent( ptid[ii], 2, sigma.zz() );
    cauchyData->InsertComponent( ptid[ii], 3, sigma.xy() );
    cauchyData->InsertComponent( ptid[ii], 4, sigma.yz() );
    cauchyData->InsertComponent( ptid[ii], 5, sigma.xz() );
  }
}

void VTK_Writer_Solid::writeOutput(
    const FEANode * const &fnode_ptr,
    const ALocal_IEN * const &lien_ptr,
    const ALocal_Elem * const &lelem_ptr,
    const IVisDataPrep * const &vdata_ptr,
    FEAElement * const &elemptr,
    const IQuadPts * const &quad,
    const double * const * const &pointArrays,
    const int &rank, const int &size,
    const double &sol_time,
    const std::string &outputBName,
    const std::string &outputName,
    const bool &isXML,
    const bool &is_ref )
{
  // Allocate VTK gridData object
  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();

  // Allocate vtk point
  vtkPoints * points = vtkPoints::New();

  const int numDArrays = vdata_ptr->get_arrayCompSize();
  if(numDArrays != 5) SYS_T::print_fatal("Error: vdata size numDArrays != 5. \n");

  vtkDoubleArray ** dataVecs = new vtkDoubleArray * [numDArrays];
  for(int ii=0; ii<numDArrays; ++ii)
  {
    dataVecs[ii] = vtkDoubleArray::New();
    dataVecs[ii] -> SetNumberOfComponents( vdata_ptr->get_arraySizes(ii) );
    dataVecs[ii] -> SetName( vdata_ptr->get_arrayNames(ii).c_str() );
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

    // Displacement
    std::vector<double> dispInfo; dispInfo.clear();
    int asize = vdata_ptr->get_arraySizes(0);
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        dispInfo.push_back( pointArrays[0][pt_index * asize + kk ] );
    }
    Interp::VTKData( asize, &IEN_e[0], &dispInfo[0], elemptr, dataVecs[0] );

    if( is_ref )
      Interp::VTKPts( &IEN_e[0], &ectrl_x[0], &ectrl_y[0], &ectrl_z[0], elemptr, points );
    else
      Interp::VTKPts( &IEN_e[0], &ectrl_x[0], &ectrl_y[0], &ectrl_z[0],
          &dispInfo[0], elemptr, points );

    // Pressure
    std::vector<double> presInfo; presInfo.clear();
    asize = vdata_ptr->get_arraySizes(2);
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        presInfo.push_back( pointArrays[1][pt_index * asize + kk ] );
    }
    Interp::VTKData( asize, &IEN_e[0], &presInfo[0], elemptr, dataVecs[2] );

    // Velocity
    std::vector<double> veloInfo; veloInfo.clear();
    asize = vdata_ptr->get_arraySizes(3);
    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        veloInfo.push_back( pointArrays[2][pt_index * asize + kk ] );
    }
    Interp::VTKData( asize, &IEN_e[0], &veloInfo[0], elemptr, dataVecs[3] );

    // detF and Cauchy stress
    interpolateJ_Cauchy( &IEN_e[0], &dispInfo[0], &presInfo[0],
        elemptr, dataVecs[1], dataVecs[4] );

    // Set mesh connectivity
    if( elemptr->get_Type() == FEType::Tet4 )
      VIS_T::setTetraelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3], gridData );
    else if( elemptr->get_Type() == FEType::Tet10 )
      VIS_T::setQuadTetraelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3],
          IEN_e[4], IEN_e[5], IEN_e[6], IEN_e[7], IEN_e[8], IEN_e[9], gridData );
    else if( elemptr->get_Type() == FEType::Hex8 )
      VIS_T::setHexelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3],
          IEN_e[4], IEN_e[5], IEN_e[6], IEN_e[7], gridData );
    else if( elemptr->get_Type() == FEType::Hex27 )
      VIS_T::setTriQuadHexelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3],
          IEN_e[4], IEN_e[5], IEN_e[6], IEN_e[7], IEN_e[8], IEN_e[9],
          IEN_e[10], IEN_e[11], IEN_e[12], IEN_e[13], IEN_e[14], IEN_e[15],
          IEN_e[16], IEN_e[17], IEN_e[18], IEN_e[19], IEN_e[20], IEN_e[21],
          IEN_e[22], IEN_e[23], IEN_e[24], IEN_e[25], IEN_e[26], gridData );
    else SYS_T::print_fatal("Error: unknown element type.\n");

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

  VIS_T::writeVisFile( gridData, outputBName, outputName,
      rank, size, sol_time, isXML );

  SYS_T::commPrint("-- Clean gridData object.\n");
  gridData->Delete();
}

// EOF
