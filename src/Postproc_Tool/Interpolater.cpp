#include "Interpolater.hpp"

Interpolater::Interpolater( const int &in_nlocbas, const bool &in_isDer )
: isDer(in_isDer), nLocBas( in_nlocbas )
{
  R.resize(nLocBas);

  if(isDer)
  {
    dR_dx.resize(nLocBas);
    dR_dy.resize(nLocBas);
    dR_dz.resize(nLocBas);
  }
  else
  {
    VEC_T::clean(dR_dx);
    VEC_T::clean(dR_dy);
    VEC_T::clean(dR_dz);
  }
}

Interpolater::~Interpolater()
{
  VEC_T::clean(R);
  if(isDer)
  {
    VEC_T::clean(dR_dx);
    VEC_T::clean(dR_dy);
    VEC_T::clean(dR_dz);
  }
}

void Interpolater::print_info() const
{
  std::cout<<"Interpolater: \n";
  std::cout<<" -- nLocBas = "<<nLocBas<<std::endl;
  if(isDer) std::cout<<" -- derivatives are needed. \n";
  else std::cout<<" -- derivatives are unneeded. \n";

  std::cout<<"R: \n";
  VEC_T::print(R);
  std::cout<<"dR_dx: \n";
  VEC_T::print(dR_dx);
  std::cout<<"dR_dy: \n";
  VEC_T::print(dR_dy);
  std::cout<<"dR_dz: \n";
  VEC_T::print(dR_dz);
}

Interpolater::Interpolater()
: isDer(false), nLocBas(0)
{}

void Interpolater::interpolateFE( const double * const &inputVal,
    const FEAElement * const &elem, std::vector<double> &output )
{
  const int nqp = elem->get_numQuapts();
  output.resize( nqp );

  for(int ii=0; ii<nqp; ++ii)
  {
    output[ii] = 0.0;
    elem->get_R( ii, &R[0] );
    for(int jj=0; jj<nLocBas; ++jj) output[ii] += inputVal[jj] * R[jj];
  }
}

void Interpolater::interpolateFE( const double * const &inputVal_1,
    const double * const &inputVal_2,
    const FEAElement * const &elem,
    std::vector<double> &output_1, std::vector<double> &output_2 )
{
  const int nqp = elem->get_numQuapts();

  output_1.resize( nqp );
  output_2.resize( nqp );

  for(int ii=0; ii<nqp; ++ii)
  {
    output_1[ii] = 0.0;
    output_2[ii] = 0.0;
    elem->get_R( ii, &R[0] );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      output_1[ii] += inputVal_1[jj] * R[jj];
      output_2[ii] += inputVal_2[jj] * R[jj];
    }
  }
}

void Interpolater::interpolateFE( const double * const &inputVal_1,
    const double * const &inputVal_2,
    const double * const &inputVal_3,
    const FEAElement * const &elem,
    std::vector<double> &output_1, std::vector<double> &output_2,
    std::vector<double> &output_3 )
{
  const int nqp = elem->get_numQuapts();

  output_1.resize( nqp );
  output_2.resize( nqp );
  output_3.resize( nqp );

  for(int ii=0; ii<nqp; ++ii)
  {
    output_1[ii] = 0.0;
    output_2[ii] = 0.0;
    output_3[ii] = 0.0;
    elem->get_R( ii, &R[0] );
    for(int jj=0; jj<nLocBas; ++jj)
    {
      output_1[ii] += inputVal_1[jj] * R[jj];
      output_2[ii] += inputVal_2[jj] * R[jj];
      output_3[ii] += inputVal_3[jj] * R[jj];
    }
  }
}

void Interpolater::interpolateFE_Grad( const double * const &inputVal,
    const FEAElement * const &elem, std::vector<double> &output_dx,
    std::vector<double> &output_dy, std::vector<double> &output_dz )
{
  if( !isDer ) SYS_T::print_fatal("Interpolater : isDer is false, dR_dx, dR_dy, and dR_dz are not allocated. \n");

  const int nqp = elem->get_numQuapts();

  output_dx.resize( nqp );
  output_dy.resize( nqp );
  output_dz.resize( nqp );

  for( int ii=0; ii<nqp; ++ii )
  {
    output_dx[ii] = 0.0;
    output_dy[ii] = 0.0;
    output_dz[ii] = 0.0;

    elem->get_gradR(ii, &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    for( int jj=0; jj<nLocBas; ++jj )
    {
      output_dx[ii] += inputVal[jj] * dR_dx[jj];
      output_dy[ii] += inputVal[jj] * dR_dy[jj];
      output_dz[ii] += inputVal[jj] * dR_dz[jj];
    }
  }
}

void Interpolater::interpolateVTKPts( const int &ptoffset,
    const double * const &ctrlPts_x,
    const double * const &ctrlPts_y,
    const double * const &ctrlPts_z,
    const FEAElement * const &elem,
    vtkPoints * const &vtkpts )
{
  const int nqp = elem->get_numQuapts();
  std::vector<double> out_x, out_y, out_z;

  interpolateFE(ctrlPts_x, ctrlPts_y, ctrlPts_z, elem, out_x, out_y, out_z);
  
  double xyz[3];

  for(int ii=0; ii<nqp; ++ii)
  {
    xyz[0] = out_x[ii];
    xyz[1] = out_y[ii];
    xyz[2] = out_z[ii];

    vtkpts->InsertPoint(ptoffset+ii, xyz);
  }
}

void Interpolater::interpolateVTKPts( const int &ptoffset,
    const double * const &ctrlPts_x,
    const double * const &ctrlPts_y,
    const FEAElement * const &elem,
    vtkPoints * const &vtkpts )
{
  const int nqp = elem->get_numQuapts();
  std::vector<double> out_x, out_y;

  interpolateFE(ctrlPts_x, ctrlPts_y, elem, out_x, out_y);
  
  double xyz[2];

  for(int ii=0; ii<nqp; ++ii)
  {
    xyz[0] = out_x[ii];
    xyz[1] = out_y[ii];

    vtkpts->InsertPoint(ptoffset+ii, xyz);
  }
}


void Interpolater::interpolateVTKPts( const int * const &ptid,
    const double * const &ctrlPts_x,
    const double * const &ctrlPts_y,
    const double * const &ctrlPts_z,
    const FEAElement * const &elem,
    vtkPoints * const &vtkpts )
{
  const int nqp = elem->get_numQuapts();
  
  std::vector<double> out_x, out_y, out_z;

  interpolateFE(ctrlPts_x, ctrlPts_y, ctrlPts_z, elem, out_x, out_y, out_z);
  
  double xyz[3];

  for(int ii=0; ii<nqp; ++ii)
  {
    xyz[0] = out_x[ii];
    xyz[1] = out_y[ii];
    xyz[2] = out_z[ii];

    vtkpts->InsertPoint(ptid[ii], xyz);
  }
}

void Interpolater::interpolateVTKData( const int &size, const int &ptoffset,
    const double * const &inputData, const FEAElement * const &elem,
    vtkDoubleArray * const &vtkData )
{
  int nqp = elem->get_numQuapts();

  double * compData = new double [nLocBas];
  std::vector<double> outData;

  for(int comp=0; comp<size; ++comp)
  {
    for(int ii=0; ii<nLocBas; ++ii)
      compData[ii] = inputData[ii*size + comp]; // extract the comp-th component
    
    // FE Interpolation of the component 
    interpolateFE(compData, elem, outData);
    
    // Insert the comp-th component into vtkData
    for(int ii=0; ii<nqp; ++ii)
      vtkData->InsertComponent(ptoffset+ii, comp, outData[ii]);
  }

  delete [] compData;
}


void Interpolater::interpolateVTKData( const int &size, 
    const int * const &ptid, const double * const &inputData, 
    const FEAElement * const &elem, vtkDoubleArray * const &vtkData )
{
  int nqp = elem->get_numQuapts();

  double * compData = new double [nLocBas];
  std::vector<double> outData;

  for(int comp=0; comp<size; ++comp)
  {
    for(int ii=0; ii<nLocBas; ++ii)
      compData[ii] = inputData[ii*size + comp]; // extract the comp-th component
    
    // FE Interpolation of the component 
    interpolateFE(compData, elem, outData);
    
    // Insert the comp-th component into vtkData
    for(int ii=0; ii<nqp; ++ii)
      vtkData->InsertComponent(ptid[ii], comp, outData[ii]);
  }

  delete [] compData;
}


void Interpolater::interpolateData( const int &size,
    const double * const &inputData, const FEAElement * const &elem,
    std::vector< std::vector<double> > &outData )
{
  outData.resize( size );

  double * compData = new double [nLocBas];

  for(int comp=0; comp<size; ++comp)
  {
    // extract the comp-th component
    for(int ii=0; ii<nLocBas; ++ii)
      compData[ii] = inputData[ii*size + comp];
    
    // FE Interpolation of the component 
    interpolateFE( compData, elem, outData[comp] );
  }

  delete [] compData; compData = NULL;
}


void Interpolater::interpolateVTKPts( const int &ptoffset,
    const double * const &ctrlPts_x,
    const double * const &ctrlPts_y,
    const double * const &ctrlPts_z,
    const double * const &disp_vect,
    const FEAElement * const &elem,
    vtkPoints * const &vtkpts )
{
  const int nqp = elem->get_numQuapts();
  std::vector<double> disp_x, disp_y, disp_z, ref_x, ref_y, ref_z;

  interpolateFE(ctrlPts_x, ctrlPts_y, ctrlPts_z, elem, ref_x, ref_y, ref_z);
 
  double * compData = new double [nLocBas];

  // Get the x-component of displacement
  for(int ii=0; ii<nLocBas; ++ii) compData[ii] = disp_vect[ii*3+0];
  interpolateFE(compData, elem, disp_x);

  // Get the y-component of displacement  
  for(int ii=0; ii<nLocBas; ++ii) compData[ii] = disp_vect[ii*3+1];
  interpolateFE(compData, elem, disp_y);

  // Get the z-component of displacement
  for(int ii=0; ii<nLocBas; ++ii) compData[ii] = disp_vect[ii*3+2];
  interpolateFE(compData, elem, disp_z);

  delete [] compData; compData = NULL;

  double xyz[3];

  for(int ii=0; ii<nqp; ++ii)
  {
    xyz[0] = ref_x[ii] + disp_x[ii];
    xyz[1] = ref_y[ii] + disp_y[ii];
    xyz[2] = ref_z[ii] + disp_z[ii];

    vtkpts->InsertPoint(ptoffset+ii, xyz);
  }
}


void Interpolater::interpolateVTKPts( const int * const &ptid,
    const double * const &ctrlPts_x,
    const double * const &ctrlPts_y,
    const double * const &ctrlPts_z,
    const double * const &disp_vect,
    const FEAElement * const &elem,
    vtkPoints * const &vtkpts )
{
  const int nqp = elem->get_numQuapts();
  std::vector<double> disp_x, disp_y, disp_z, ref_x, ref_y, ref_z;

  interpolateFE(ctrlPts_x, ctrlPts_y, ctrlPts_z, elem, ref_x, ref_y, ref_z);
 
  double * compData = new double [nLocBas];

  // Get the x-component of displacement
  for(int ii=0; ii<nLocBas; ++ii) compData[ii] = disp_vect[ii*3+0];
  interpolateFE(compData, elem, disp_x);

  // Get the y-component of displacement  
  for(int ii=0; ii<nLocBas; ++ii) compData[ii] = disp_vect[ii*3+1];
  interpolateFE(compData, elem, disp_y);

  // Get the z-component of displacement
  for(int ii=0; ii<nLocBas; ++ii) compData[ii] = disp_vect[ii*3+2];
  interpolateFE(compData, elem, disp_z);

  delete [] compData; compData = NULL;

  double xyz[3];

  for(int ii=0; ii<nqp; ++ii)
  {
    xyz[0] = ref_x[ii] + disp_x[ii];
    xyz[1] = ref_y[ii] + disp_y[ii];
    xyz[2] = ref_z[ii] + disp_z[ii];

    vtkpts->InsertPoint(ptid[ii], xyz);
  }
}

// EOF
