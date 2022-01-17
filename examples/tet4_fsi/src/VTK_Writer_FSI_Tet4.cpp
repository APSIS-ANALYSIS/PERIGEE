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

// EOF
