#include "VTK_Writer_Lap.hpp"

VTK_Writer_Lap::VTK_Writer_Lap( const int &in_nelem,
    const int &in_nlocbas, const std::string &epart_file )
: nLocBas( in_nlocbas ), nElem( in_nelem )
{
  VIS_T::read_epart( epart_file, nElem, epart_map );
}


VTK_Writer_Lap::~VTK_Writer_Lap()
{
  VEC_T::clean( epart_map );
}

void VTK_Writer_Lap::writeOutput(
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

  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();

  vtkPoints * points = vtkPoints::New();

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

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    const std::vector<int> IEN_e = lien_ptr -> get_LIEN( ee );

    std::vector<double> ectrl_x ( nLocBas, 0.0 );
    std::vector<double> ectrl_y ( nLocBas, 0.0 );
    std::vector<double> ectrl_z ( nLocBas, 0.0 );

    fnode_ptr -> get_ctrlPts_xyz(nLocBas, &IEN_e[0], &ectrl_x[0], &ectrl_y[0], &ectrl_z[0]);

    elemptr->buildBasis( quad, &ectrl_x[0], &ectrl_y[0], &ectrl_z[0] );

    intep.interpolateVTKPts(&IEN_e[0], &ectrl_x[0], &ectrl_y[0], &ectrl_z[0],
        elemptr, points );

    std::vector<double> inputInfo; inputInfo.clear();

    int asize = vdata_ptr -> get_arraySizes(0);

    for(int jj=0; jj<nLocBas; ++jj)
    {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk)
        inputInfo.push_back( pointArrays[0][pt_index * asize + kk] );
    }

    intep.interpolateVTKData( asize, &IEN_e[0], &inputInfo[0],
        elemptr, dataVecs[0] );

    if( elemptr->get_Type() == 501 )
      VIS_T::setTetraelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3], gridData );
    else if( elemptr->get_Type() == 502 )
      VIS_T::setQuadTetraelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3],
          IEN_e[4], IEN_e[5], IEN_e[6], IEN_e[7], IEN_e[8], IEN_e[9],
          gridData );
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

  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object.\n");
  gridData->Delete();
}

// EOF
