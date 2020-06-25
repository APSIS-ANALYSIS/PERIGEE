#include "PLocAssem_Tet_Heat.hpp"

PLocAssem_Tet_Heat::PLocAssem_Tet_Heat( const int &in_nlocbas, 
    const int &in_nqp, const int &in_snlocbas, 
    const int &elemtype ) : nqp(in_nqp)
{
  if(elemtype == 501)
  {
    // 501 is linear element
    nLocBas = 4; snLocBas = 3;
  }
  else if(elemtype == 502)
  {
    // 502 is quadratic element
    nLocBas = 10; snLocBas = 6;
  }
  else SYS_T::print_fatal("Error: unknown elem type.\n");

  vec_size = nLocBas * 1;
  sur_size = snLocBas * 1;

  R.resize(nLocBas);
  dR_dx.resize(nLocBas);
  dR_dy.resize(nLocBas);
  dR_dz.resize(nLocBas);

  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  sur_Tangent = new PetscScalar[sur_size * sur_size];
  sur_Residual = new PetscScalar[sur_size];

  Zero_Tangent_Residual();

  Zero_sur_Tangent_Residual();

  print_info();
}


PLocAssem_Tet_Heat::~PLocAssem_Tet_Heat()
{
  delete [] Tangent; Tangent = nullptr;
  delete [] Residual; Residual = nullptr;
  delete [] sur_Tangent; sur_Tangent = nullptr;
  delete [] sur_Residual; sur_Residual = nullptr;
}


void PLocAssem_Tet_Heat::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("  Three-dimensional Laplace equation. \n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
}


void PLocAssem_Tet_Heat::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const double curr = time + dt;

  double f;

  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_Tangent_Residual();

  for(int qua = 0; qua<nqp; ++qua)
  {
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    double d_x = 0.0, d_y = 0.0, d_z = 0.0;

    element->get_R_gradR(qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0]);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      d_x += sol[ii] * dR_dx[ii];
      d_y += sol[ii] * dR_dy[ii];
      d_z += sol[ii] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    get_f(coor_x, coor_y, coor_z, curr, f);

    for(int A=0; A<nLocBas; ++A)
    {
      Residual[A] += gwts * ( d_x * dR_dx[A] 
          + d_y * dR_dy[A] + d_z * dR_dz[A] - f * R[A] );

      for(int B=0; B<nLocBas; ++B)
      {
        Tangent[A*nLocBas + B] += gwts * ( dR_dx[A] * dR_dx[B]
            + dR_dy[A] * dR_dy[B] + dR_dz[A] * dR_dz[B] );
      }
    }
  }
}


void PLocAssem_Tet_Heat::Assem_Residual_EBC(
    const int &ebc_id,
    const double &time, const double &dt,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  const double curr = time + dt;

  double g, nx, ny, nz, surface_area;

  Zero_Residual();

  for(int qua = 0; qua < face_nqp; ++qua )
  {
    element->get_R(qua, &R[0]);
    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);
  
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    for(int ii=0; ii<snLocBas; ++ii)
    {
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }

    get_ebc_fun( ebc_id, coor_x, coor_y, coor_z, curr, nx, ny, nz, g );

    const double gwts = surface_area * quad -> get_qw( qua );

    for(int A=0; A<snLocBas; ++A)
      Residual[A] -= gwts * R[A] * g;
  }
}


void PLocAssem_Tet_Heat::Assem_Residual_EBC_Resistance(
    const int &ebc_id, const double &val,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  double nx, ny, nz, surface_area;

  Zero_Residual();

  for(int qua = 0; qua < face_nqp; ++qua)
  {
    element->get_R(qua, &R[0]);

    element->get_2d_normal_out(qua, nx, ny, nz, surface_area);

    const double gwts = surface_area * quad -> get_qw(qua);

    for(int A=0; A<snLocBas; ++A)
      Residual[A] += gwts * R[A] * val;
  }
}

// EOF
