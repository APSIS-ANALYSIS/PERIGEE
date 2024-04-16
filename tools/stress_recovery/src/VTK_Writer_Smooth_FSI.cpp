#include "VTK_Writer_Smooth_FSI.hpp"

VTK_Writer_Smooth_FSI::VTK_Writer_Smooth_FSI( const int &in_nelem,
    const int &in_nlocbas, const std::string &epart_file )
: nLocBas( in_nlocbas ), nElem( in_nelem )
{
  VIS_T::read_epart( epart_file, nElem, epart_map );
}

VTK_Writer_Smooth_FSI::~VTK_Writer_Smooth_FSI()
{
  VEC_T::clean(epart_map);
}

void VTK_Writer_Smooth_FSI::interpolateF( const int * const &ptid,
    const std::vector<double> &inputGradDisp,
    const FEAElement * const &elem,
    vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  std::vector<double> ux (nLocBas, 0.0), uy (nLocBas, 0.0), uz (nLocBas, 0.0),
                      vx (nLocBas, 0.0), vy (nLocBas, 0.0), vz (nLocBas, 0.0),
                      wx (nLocBas, 0.0), wy (nLocBas, 0.0), wz (nLocBas, 0.0);

  for(int ii=0; ii<nLocBas; ++ii)
  {
    ux[ii] = inputGradDisp[ii*9];
    uy[ii] = inputGradDisp[ii*9+1];
    uz[ii] = inputGradDisp[ii*9+2];
    vx[ii] = inputGradDisp[ii*9+3];
    vy[ii] = inputGradDisp[ii*9+4];
    vz[ii] = inputGradDisp[ii*9+5];
    wx[ii] = inputGradDisp[ii*9+6];
    wy[ii] = inputGradDisp[ii*9+7];
    wz[ii] = inputGradDisp[ii*9+8];
  }

  // interpolate reference points
  Interpolater intep( nLocBas );

  for(int ii=0; ii<nqp; ++ii)
  {
    // F
    const Tensor2_3D FF( ux[ii] + 1.0, uy[ii],       uz[ii],
                         vx[ii],       vy[ii] + 1.0, vz[ii],
                         wx[ii],       wy[ii],       wz[ii] + 1.0 );

    vtkData->InsertComponent( ptid[ii], 0, FF.xx() );
    vtkData->InsertComponent( ptid[ii], 1, FF.xy() );
    vtkData->InsertComponent( ptid[ii], 2, FF.xz() );
    vtkData->InsertComponent( ptid[ii], 3, FF.yx() );
    vtkData->InsertComponent( ptid[ii], 4, FF.yy() );
    vtkData->InsertComponent( ptid[ii], 5, FF.yz() );
    vtkData->InsertComponent( ptid[ii], 6, FF.zx() );
    vtkData->InsertComponent( ptid[ii], 7, FF.zy() );
    vtkData->InsertComponent( ptid[ii], 8, FF.zz() );
  }
}

void VTK_Writer_Smooth_FSI::interpolateStrain( const int * const &ptid,
    const std::vector<double> &inputGradDisp,
    const FEAElement * const &elem,
    IMaterialModel * const &model,
    vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  std::vector<double> ux (nLocBas, 0.0), uy (nLocBas, 0.0), uz (nLocBas, 0.0),
                      vx (nLocBas, 0.0), vy (nLocBas, 0.0), vz (nLocBas, 0.0),
                      wx (nLocBas, 0.0), wy (nLocBas, 0.0), wz (nLocBas, 0.0);

  for(int ii=0; ii<nLocBas; ++ii)
  {
    ux[ii] = inputGradDisp[ii*9];
    uy[ii] = inputGradDisp[ii*9+1];
    uz[ii] = inputGradDisp[ii*9+2];
    vx[ii] = inputGradDisp[ii*9+3];
    vy[ii] = inputGradDisp[ii*9+4];
    vz[ii] = inputGradDisp[ii*9+5];
    wx[ii] = inputGradDisp[ii*9+6];
    wy[ii] = inputGradDisp[ii*9+7];
    wz[ii] = inputGradDisp[ii*9+8];
  }

  // interpolate reference points
  Interpolater intep( nLocBas );

  for(int ii=0; ii<nqp; ++ii)
  {
    // F
    const Tensor2_3D FF( ux[ii] + 1.0, uy[ii],       uz[ii],
                         vx[ii],       vy[ii] + 1.0, vz[ii],
                         wx[ii],       wy[ii],       wz[ii] + 1.0 );

    Tensor2_3D CC; CC.MatMultTransposeLeft(FF);

    const Tensor2_3D I( 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    const Tensor2_3D EE = 0.5 * ( CC - I );

    vtkData->InsertComponent( ptid[ii], 0, EE.xx() );
    vtkData->InsertComponent( ptid[ii], 1, EE.xy() );
    vtkData->InsertComponent( ptid[ii], 2, EE.xz() );
    vtkData->InsertComponent( ptid[ii], 3, EE.yx() );
    vtkData->InsertComponent( ptid[ii], 4, EE.yy() );
    vtkData->InsertComponent( ptid[ii], 5, EE.yz() );
    vtkData->InsertComponent( ptid[ii], 6, EE.zx() );
    vtkData->InsertComponent( ptid[ii], 7, EE.zy() );
    vtkData->InsertComponent( ptid[ii], 8, EE.zz() );
  }
}

void VTK_Writer_Smooth_FSI::interpolateCauchy( const int * const &ptid,
    const std::vector<Vector_3> &eleBasis_r,
    const std::vector<Vector_3> &eleBasis_c,
    const std::vector<Vector_3> &eleBasis_l,
    const std::vector<double> &inputGradDisp,
    const FEAElement * const &elem,
    IMaterialModel * const &model,
    vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  std::vector<double> ux (nLocBas, 0.0), uy (nLocBas, 0.0), uz (nLocBas, 0.0),
                      vx (nLocBas, 0.0), vy (nLocBas, 0.0), vz (nLocBas, 0.0),
                      wx (nLocBas, 0.0), wy (nLocBas, 0.0), wz (nLocBas, 0.0);                   

  Vector_3 basis_r(0.0, 0.0, 0.0);
  Vector_3 basis_c(0.0, 0.0, 0.0);
  Vector_3 basis_l(0.0, 0.0, 0.0);

  for(int ii=0; ii<nLocBas; ++ii)
  {
    ux[ii] = inputGradDisp[ii*9];
    uy[ii] = inputGradDisp[ii*9+1];
    uz[ii] = inputGradDisp[ii*9+2];
    vx[ii] = inputGradDisp[ii*9+3];
    vy[ii] = inputGradDisp[ii*9+4];
    vz[ii] = inputGradDisp[ii*9+5];
    wx[ii] = inputGradDisp[ii*9+6];
    wy[ii] = inputGradDisp[ii*9+7];
    wz[ii] = inputGradDisp[ii*9+8];

    basis_r = eleBasis_r[ii];
    basis_c = eleBasis_c[ii];
    basis_l = eleBasis_l[ii];
  }

  Interpolater intep( nLocBas );

  for(int ii=0; ii<nqp; ++ii)
  {
    // update fibre direction by the basis vector on an arbitrary node
    //model->update_fibre_dir(basis_r, basis_c, basis_l);
    model->update_fibre_dir(eleBasis_r[ii], eleBasis_c[ii], eleBasis_l[ii]);

    // F
    const Tensor2_3D FF( ux[ii] + 1.0, uy[ii],       uz[ii],
                         vx[ii],       vy[ii] + 1.0, vz[ii],
                         wx[ii],       wy[ii],       wz[ii] + 1.0 );

    Tensor2_3D sigma = model -> get_Cauchy_stress( FF );

    vtkData->InsertComponent( ptid[ii], 0, sigma.xx() );
    vtkData->InsertComponent( ptid[ii], 1, sigma.xy() );
    vtkData->InsertComponent( ptid[ii], 2, sigma.xz() );
    vtkData->InsertComponent( ptid[ii], 3, sigma.yx() );
    vtkData->InsertComponent( ptid[ii], 4, sigma.yy() );
    vtkData->InsertComponent( ptid[ii], 5, sigma.yz() );
    vtkData->InsertComponent( ptid[ii], 6, sigma.zx() );
    vtkData->InsertComponent( ptid[ii], 7, sigma.zy() );
    vtkData->InsertComponent( ptid[ii], 8, sigma.zz() );
  }
}

void VTK_Writer_Smooth_FSI::interpolateVonMises( const int * const &ptid,
    const std::vector<Vector_3> &eleBasis_r,
    const std::vector<Vector_3> &eleBasis_c,
    const std::vector<Vector_3> &eleBasis_l,
    const std::vector<double> &inputGradDisp,
    const std::vector<double> &inputPres,
    const FEAElement * const &elem,
    IMaterialModel * const &model,
    vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  std::vector<double> ux (nLocBas, 0.0), uy (nLocBas, 0.0), uz (nLocBas, 0.0),
                      vx (nLocBas, 0.0), vy (nLocBas, 0.0), vz (nLocBas, 0.0),
                      wx (nLocBas, 0.0), wy (nLocBas, 0.0), wz (nLocBas, 0.0),
                      pp (nLocBas, 0.0);                    

  Vector_3 basis_r(0.0, 0.0, 0.0);
  Vector_3 basis_c(0.0, 0.0, 0.0);
  Vector_3 basis_l(0.0, 0.0, 0.0);

  for(int ii=0; ii<nLocBas; ++ii)
  {
    ux[ii] = inputGradDisp[ii*9];
    uy[ii] = inputGradDisp[ii*9+1];
    uz[ii] = inputGradDisp[ii*9+2];
    vx[ii] = inputGradDisp[ii*9+3];
    vy[ii] = inputGradDisp[ii*9+4];
    vz[ii] = inputGradDisp[ii*9+5];
    wx[ii] = inputGradDisp[ii*9+6];
    wy[ii] = inputGradDisp[ii*9+7];
    wz[ii] = inputGradDisp[ii*9+8];

    pp[ii] = inputPres[ii];

    basis_r = eleBasis_r[ii];
    basis_c = eleBasis_c[ii];
    basis_l = eleBasis_l[ii];
  }

  Interpolater intep( nLocBas );

  for(int ii=0; ii<nqp; ++ii)
  {
    // update fibre direction by the basis vector on an arbitrary node
    //model->update_fibre_dir(basis_r, basis_c, basis_l);
    model->update_fibre_dir(eleBasis_r[ii], eleBasis_c[ii], eleBasis_l[ii]);

    // F
    const Tensor2_3D FF( ux[ii] + 1.0, uy[ii],       uz[ii],
                         vx[ii],       vy[ii] + 1.0, vz[ii],
                         wx[ii],       wy[ii],       wz[ii] + 1.0 );

    Tensor2_3D sigma = model -> get_Cauchy_stress( FF );

    Tensor2_3D pres( -pp[ii], 0.0, 0.0, 0.0, -pp[ii], 0.0, 0.0, 0.0, -pp[ii]);

    Tensor2_3D sigma_p = sigma + pres;

    // Principal stress
    double eta1 = 0.0; 
    double eta2 = 0.0;
    double eta3 = 0.0;
    Vector_3 vec1(0.0, 0.0, 0.0); 
    Vector_3 vec2(0.0, 0.0, 0.0); 
    Vector_3 vec3(0.0, 0.0, 0.0); 
    Vector_3 temp_vec(0.0, 0.0, 0.0);
    
    const int num_dif_eigen = sigma_p.eigen_decomp(eta1, eta2, eta3, vec1, vec2, vec3);

    double sigma_v = 0.5 * std::sqrt(2.0) * std::sqrt( (eta1 - eta2) * (eta1 - eta2)
                    + (eta2 - eta3) * (eta2 - eta3) + (eta3 - eta1) * (eta3 - eta1) );

    vtkData->InsertComponent( ptid[ii], 0, sigma_v );

  }
}

void VTK_Writer_Smooth_FSI::interpolateVonMises_nop( const int * const &ptid,
    const std::vector<Vector_3> &eleBasis_r,
    const std::vector<Vector_3> &eleBasis_c,
    const std::vector<Vector_3> &eleBasis_l,
    const std::vector<double> &inputGradDisp,
    const FEAElement * const &elem,
    IMaterialModel * const &model,
    vtkDoubleArray * const &vtkData )
{
  const int nqp = elem->get_numQuapts();

  std::vector<double> ux (nLocBas, 0.0), uy (nLocBas, 0.0), uz (nLocBas, 0.0),
                      vx (nLocBas, 0.0), vy (nLocBas, 0.0), vz (nLocBas, 0.0),
                      wx (nLocBas, 0.0), wy (nLocBas, 0.0), wz (nLocBas, 0.0),
                      pp (nLocBas, 0.0);                    

  Vector_3 basis_r(0.0, 0.0, 0.0);
  Vector_3 basis_c(0.0, 0.0, 0.0);
  Vector_3 basis_l(0.0, 0.0, 0.0);

  for(int ii=0; ii<nLocBas; ++ii)
  {
    ux[ii] = inputGradDisp[ii*9];
    uy[ii] = inputGradDisp[ii*9+1];
    uz[ii] = inputGradDisp[ii*9+2];
    vx[ii] = inputGradDisp[ii*9+3];
    vy[ii] = inputGradDisp[ii*9+4];
    vz[ii] = inputGradDisp[ii*9+5];
    wx[ii] = inputGradDisp[ii*9+6];
    wy[ii] = inputGradDisp[ii*9+7];
    wz[ii] = inputGradDisp[ii*9+8];

    basis_r = eleBasis_r[ii];
    basis_c = eleBasis_c[ii];
    basis_l = eleBasis_l[ii];
  }

  Interpolater intep( nLocBas );

  for(int ii=0; ii<nqp; ++ii)
  {
    // update fibre direction by the basis vector on an arbitrary node
    //model->update_fibre_dir(basis_r, basis_c, basis_l);
    model->update_fibre_dir(eleBasis_r[ii], eleBasis_c[ii], eleBasis_l[ii]);

    // F
    const Tensor2_3D FF( ux[ii] + 1.0, uy[ii],       uz[ii],
                         vx[ii],       vy[ii] + 1.0, vz[ii],
                         wx[ii],       wy[ii],       wz[ii] + 1.0 );

    Tensor2_3D sigma = model -> get_Cauchy_stress( FF );

    // Principal stress
    double eta1 = 0.0; 
    double eta2 = 0.0;
    double eta3 = 0.0;
    Vector_3 vec1(0.0, 0.0, 0.0); 
    Vector_3 vec2(0.0, 0.0, 0.0); 
    Vector_3 vec3(0.0, 0.0, 0.0); 
    Vector_3 temp_vec(0.0, 0.0, 0.0);
    
    //const int num_dif_eigen = sigma_p.eigen_decomp(eta1, eta2, eta3, vec1, vec2, vec3);
    const int num_dif_eigen = sigma.eigen_decomp(eta1, eta2, eta3, vec1, vec2, vec3);

    double sigma_v = 0.5 * std::sqrt(2.0) * std::sqrt( (eta1 - eta2) * (eta1 - eta2)
                    + (eta2 - eta3) * (eta2 - eta3) + (eta3 - eta1) * (eta3 - eta1) );

    vtkData->InsertComponent( ptid[ii], 0, sigma_v );

  }
}

void VTK_Writer_Smooth_FSI::writeOutput(
    const FEANode * const &fnode_ptr,
    	const ALocal_IEN * const &lien_v,
    	const ALocal_IEN * const &lien_p,
    	const std::vector<int> &sien,
    	const ALocal_Elem * const &lelem_ptr,
    	const IVisDataPrep * const &vdata_ptr,
      IMaterialModel * const &matmodel,
    	FEAElement * const &elemptr,
    	const IQuadPts * const &quad,
      const Tissue_property * const &tp_ptr,
    	const double * const * const &pointArrays,
    	const int &rank, const int &size,
    	const int &num_of_nodes,
    	const double &sol_time,
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
    dataVecs[ii] -> SetNumberOfTuples( num_of_nodes );
  }

  // An integer array for the analysis mesh partition
  vtkIntArray * anaprocId = vtkIntArray::New();
  anaprocId -> SetName("Analysis_Partition");
  anaprocId -> SetNumberOfComponents(1);

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
  {
    if( lelem_ptr->get_elem_tag(ee) == 1)
    {
      const std::vector<int> IEN_v = lien_v -> get_LIEN( ee );
      const std::vector<int> IEN_p = lien_p -> get_LIEN( ee );

      std::vector<double> ectrl_x ( nLocBas, 0.0 );
      std::vector<double> ectrl_y ( nLocBas, 0.0 );
      std::vector<double> ectrl_z ( nLocBas, 0.0 );

      fnode_ptr -> get_ctrlPts_xyz(nLocBas, &IEN_v[0], &ectrl_x[0], &ectrl_y[0], &ectrl_z[0]);

      // get fibre direction basis
      std::vector<Vector_3> ebasis_r(nLocBas);
      std::vector<Vector_3> ebasis_c(nLocBas);
      std::vector<Vector_3> ebasis_l(nLocBas);

      for(int ii=0; ii<nLocBas; ++ii)
      {
        ebasis_r[ii] = tp_ptr -> get_basis_r(IEN_v[ii]);
        ebasis_c[ii] = tp_ptr -> get_basis_c(IEN_v[ii]);
        ebasis_l[ii] = tp_ptr -> get_basis_l(IEN_v[ii]);
      }

      elemptr->buildBasis( quad, &ectrl_x[0], &ectrl_y[0], &ectrl_z[0] );

      std::vector<int> IEN_s(nLocBas, 0.0);

      for(int ii=0; ii<nLocBas; ++ii) IEN_s[ii] = sien[ee * nLocBas + ii];
    
      // Interpolate data and assign to dataVecs
      std::vector<double> inputInfo_d; inputInfo_d.clear();
      int asize = vdata_ptr -> get_arraySizes(7);

      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_v[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo_d.push_back( pointArrays[2][pt_index * asize + kk] );
      }

      // displacement interpolation
      intep.interpolateVTKData( asize, &IEN_s[0], inputInfo_d, elemptr, dataVecs[7] );

      // use displacement to update points
      intep.interpolateVTKPts( &IEN_s[0], &ectrl_x[0], &ectrl_y[0], &ectrl_z[0], inputInfo_d, elemptr, points );

      std::vector<double> inputInfo_p; inputInfo_p.clear();
      asize = vdata_ptr->get_arraySizes(1);  
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_p[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo_p.push_back( pointArrays[1][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, &IEN_s[0], inputInfo_p, elemptr, dataVecs[1] );

      std::vector<double> inputInfo_grad; inputInfo_grad.clear();
      asize = vdata_ptr->get_arraySizes(3);
      for(int jj=0; jj<nLocBas; ++jj)
      {
        int pt_index = IEN_v[jj];
        for(int kk=0; kk<asize; ++kk)
          inputInfo_grad.push_back( pointArrays[0][pt_index * asize + kk ] );
      }
      intep.interpolateVTKData( asize, &IEN_s[0], &inputInfo_grad[0], elemptr, dataVecs[3] );

      // Interpolate F
      interpolateF( &IEN_s[0], inputInfo_grad, elemptr, dataVecs[4] );

      interpolateCauchy( &IEN_s[0], ebasis_r, ebasis_c, ebasis_l, inputInfo_grad, elemptr, matmodel, dataVecs[0] );

      interpolateVonMises( &IEN_s[0], ebasis_r, ebasis_c, ebasis_l, inputInfo_grad, inputInfo_p, elemptr, matmodel, dataVecs[2] );

      interpolateVonMises_nop( &IEN_s[0], ebasis_r, ebasis_c, ebasis_l, inputInfo_grad, elemptr, matmodel, dataVecs[5] );

      interpolateStrain( &IEN_s[0], inputInfo_grad, elemptr, matmodel, dataVecs[6] );

      // Set mesh connectivity
      if( elemptr->get_Type() == 501 )
        VIS_T::setTetraelem( IEN_s[0], IEN_s[1], IEN_s[2], IEN_s[3], gridData );
      else if( elemptr->get_Type() == 601 )
        VIS_T::setHexelem( IEN_s[0], IEN_s[1], IEN_s[2], IEN_s[3], 
          IEN_s[4], IEN_s[5], IEN_s[6], IEN_s[7], gridData );
      else SYS_T::print_fatal("Error: unknown element type.\n");

      // Mesh partition info
      anaprocId->InsertNextValue( epart_map[ lelem_ptr->get_elem_loc(ee) ] );
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

// EOF