#include "FEAElement_Tet4.hpp"

FEAElement_Tet4::FEAElement_Tet4( const int &in_nqua ) : numQuapts( in_nqua )
{
  R = new double [4 * numQuapts];
}

FEAElement_Tet4::~FEAElement_Tet4()
{
  delete [] R; R = nullptr;
}

void FEAElement_Tet4::print_info() const
{
  SYS_T::commPrint("Tet4: ");
  SYS_T::commPrint("4-node tetrahedral element with up to 2nd derivatives. \n");
  SYS_T::commPrint("elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}

double FEAElement_Tet4::get_memory_usage() const
{
  const double double_size = 4 * numQuapts + 31;
  const double int_size = 1;
  return double_size * 8.0 + int_size * 4.0;
}

void FEAElement_Tet4::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  ASSERT( quad -> get_dim() == 4, "FEAElement_Tet4::buildBasis function error.\n" );

  // area coordinates, the rest one is  qua_u = 1.0 - qua_r - qua_s - qua_t
  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );
    const double qua_t = quad -> get_qp( qua, 2 );

    R[qua*4+0] = 1.0 - qua_r - qua_s - qua_t;
    R[qua*4+1] = qua_r;
    R[qua*4+2] = qua_s;
    R[qua*4+3] = qua_t;
  }
  
  Jac[0] = -ctrl_x[0] + ctrl_x[1]; // dx_dr
  Jac[1] = -ctrl_x[0] + ctrl_x[2]; // dx_ds
  Jac[2] = -ctrl_x[0] + ctrl_x[3]; // dx_dt
  
  Jac[3] = -ctrl_y[0] + ctrl_y[1]; // dy_dr
  Jac[4] = -ctrl_y[0] + ctrl_y[2]; // dy_ds
  Jac[5] = -ctrl_y[0] + ctrl_y[3]; // dy_dt

  Jac[6] = -ctrl_z[0] + ctrl_z[1]; // dz_dr
  Jac[7] = -ctrl_z[0] + ctrl_z[2]; // dz_ds
  Jac[8] = -ctrl_z[0] + ctrl_z[3]; // dz_dt

  // Make use of the existing 3x3 matrix tool
  FE_T::Matrix_double_3by3_Array mdrdx( Jac[0], Jac[1], Jac[2], Jac[3], Jac[4],
      Jac[5], Jac[6], Jac[7], Jac[8] );

  // detJac = det(dx/dr)
  detJac = mdrdx.det();

  mdrdx.inverse();

  // dr_dx
  Jac[9]  = mdrdx(0); Jac[10] = mdrdx(1); Jac[11] = mdrdx(2);
  
  // ds_dx
  Jac[12] = mdrdx(3); Jac[13] = mdrdx(4); Jac[14] = mdrdx(5);
  
  // dt_dx
  Jac[15] = mdrdx(6); Jac[16] = mdrdx(7); Jac[17] = mdrdx(8);

  dR_dx[0] = -Jac[9] - Jac[12] - Jac[15];
  dR_dx[1] = Jac[9];
  dR_dx[2] = Jac[12];
  dR_dx[3] = Jac[15];

  dR_dy[0] = -Jac[10] - Jac[13] - Jac[16];
  dR_dy[1] = Jac[10];
  dR_dy[2] = Jac[13];
  dR_dy[3] = Jac[16];

  dR_dz[0] = -Jac[11] - Jac[14] - Jac[17];
  dR_dz[1] = Jac[11];
  dR_dz[2] = Jac[14];
  dR_dz[3] = Jac[17];
}

void FEAElement_Tet4::get_R( const int &quaindex, double * const &basis ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_R function error.\n" );
  const int offset = quaindex * 4;
  basis[0] = R[offset];
  basis[1] = R[offset+1];
  basis[2] = R[offset+2];
  basis[3] = R[offset+3];
}

std::vector<double> FEAElement_Tet4::get_R( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_R function error.\n" );
  const int offset = quaindex * 4;
  return { R[offset], R[offset+1], R[offset+2], R[offset+3] };
}

void FEAElement_Tet4::get_gradR( const int &quaindex, double * const &basis_x,
    double * const &basis_y, double * const &basis_z ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_gradR function error.\n" );
  for( int ii=0; ii<4; ++ii )
  {
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
    basis_z[ii] = dR_dz[ii];
  }
}

std::vector<double> FEAElement_Tet4::get_dR_dx( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_dR_dx function error.\n" );
  return { dR_dx[0], dR_dx[1], dR_dx[2], dR_dx[3] };
}

std::vector<double> FEAElement_Tet4::get_dR_dy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_dR_dy function error.\n" );
  return { dR_dy[0], dR_dy[1], dR_dy[2], dR_dy[3] };
}

std::vector<double> FEAElement_Tet4::get_dR_dz( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_dR_dz function error.\n" );
  return { dR_dz[0], dR_dz[1], dR_dz[2], dR_dz[3] };
}

void FEAElement_Tet4::get_R_gradR( const int &quaindex, double * const &basis,
    double * const &basis_x, double * const &basis_y,
    double * const &basis_z ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_R_gradR function error.\n" );
  const int offset = quaindex * 4;
  for( int ii=0; ii<4; ++ii )
  {
    basis[ii] = R[offset + ii];
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
    basis_z[ii] = dR_dz[ii];
  }
}

void FEAElement_Tet4::get_3D_R_dR_d2R( const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y, double * const &basis_z,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_zz, double * const &basis_xy,
    double * const &basis_xz, double * const &basis_yz ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_3D_R_dR_d2R function error.\n" );
  const int offset = quaindex * 4;
  for( int ii=0; ii<4; ++ii )
  {
    basis[ii] = R[offset + ii];
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
    basis_z[ii] = dR_dz[ii];
    basis_xx[ii] = 0.0;
    basis_yy[ii] = 0.0;
    basis_zz[ii] = 0.0;
    basis_xy[ii] = 0.0;
    basis_xz[ii] = 0.0;
    basis_yz[ii] = 0.0;
  }
}

void FEAElement_Tet4::get_3D_R_gradR_LaplacianR( const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y, double * const &basis_z,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_zz ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet4::get_3D_R_gradR_LaplacianR function error.\n" );
  const int offset = quaindex * 4;
  for( int ii=0; ii<4; ++ii )
  {
    basis[ii] = R[offset + ii];
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
    basis_z[ii] = dR_dz[ii];
    basis_xx[ii] = 0.0;
    basis_yy[ii] = 0.0;
    basis_zz[ii] = 0.0;
  }
}

void FEAElement_Tet4::get_Jacobian(const int &quaindex,
    double * const &jac_value) const
{
  for( int ii=0; ii<9; ++ii ) jac_value[ii] = Jac[ii];
}

void FEAElement_Tet4::get_invJacobian(const int &quaindex,
    double * const &jac_value) const
{
  for(int ii=0; ii<9; ++ii) jac_value[ii] = Jac[9+ii];
}

double FEAElement_Tet4::get_h( const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z ) const
{
  double x,y,z,r;

  FE_T::get_tet_sphere_info(
      ctrl_x[0], ctrl_x[1], ctrl_x[2], ctrl_x[3],
      ctrl_y[0], ctrl_y[1], ctrl_y[2], ctrl_y[3],
      ctrl_z[0], ctrl_z[1], ctrl_z[2], ctrl_z[3],
      x, y, z, r );

  return 2.0 * r;
}

// EOF
