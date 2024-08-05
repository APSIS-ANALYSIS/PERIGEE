#include "FEAElement_MINI_P1_3D.hpp"

FEAElement_MINI_P1_3D::FEAElement_MINI_P1_3D( const int &in_nqp )
: numQuapts( in_nqp )
{
  R = new double [5*numQuapts];
  dB_dx = new double [numQuapts];
  dB_dy = new double [numQuapts];
  dB_dz = new double [numQuapts];
}

FEAElement_MINI_P1_3D::~FEAElement_MINI_P1_3D()
{
  delete [] R; R = nullptr;
  delete [] dB_dx; dB_dx = nullptr;
  delete [] dB_dy; dB_dy = nullptr;
  delete [] dB_dz; dB_dz = nullptr;
}

void FEAElement_MINI_P1_3D::print_info() const
{
  SYS_T::commPrint("MINI_P1: ");
  SYS_T::commPrint("4-node tet element with quartic bubble. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}

double FEAElement_MINI_P1_3D::get_memory_usage() const
{
  const double double_size = 8 * numQuapts + 31;
  const double int_size = 1;
  return double_size * 8.0 + int_size * 4.0;
}

double FEAElement_MINI_P1_3D::get_h( const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z ) const
{
  return 2.0 * FE_T::get_tet_sphere_radius(
      ctrl_x[0], ctrl_x[1], ctrl_x[2], ctrl_x[3],
      ctrl_y[0], ctrl_y[1], ctrl_y[2], ctrl_y[3],
      ctrl_z[0], ctrl_z[1], ctrl_z[2], ctrl_z[3] );  
}

void FEAElement_MINI_P1_3D::get_R( const int &quaindex, 
    double * const &basis ) const
{
  const int offset = quaindex * 5;
  basis[0] = R[offset];
  basis[1] = R[offset+1];
  basis[2] = R[offset+2];
  basis[3] = R[offset+3];
  basis[4] = R[offset+4];
}

void FEAElement_MINI_P1_3D::get_gradR( const int &quaindex, 
    double * const &basis_x, double * const &basis_y, 
    double * const &basis_z ) const
{
  for( int ii=0; ii<4; ++ii )
  {
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
    basis_z[ii] = dR_dz[ii];
  }

  basis_x[4] = dB_dx[quaindex];
  basis_y[4] = dB_dy[quaindex];
  basis_z[4] = dB_dz[quaindex];
}

void FEAElement_MINI_P1_3D::get_R_gradR( const int &quaindex, 
    double * const &basis, double * const &basis_x, 
    double * const &basis_y, double * const &basis_z ) const
{
  const int offset = quaindex * 5;
  basis[0] = R[offset];
  basis[1] = R[offset+1];
  basis[2] = R[offset+2];
  basis[3] = R[offset+3];
  basis[4] = R[offset+4];

  for( int ii=0; ii<4; ++ii )
  {
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
    basis_z[ii] = dR_dz[ii];
  }

  basis_x[4] = dB_dx[quaindex];
  basis_y[4] = dB_dy[quaindex];
  basis_z[4] = dB_dz[quaindex];
}

void FEAElement_MINI_P1_3D::get_Jacobian(const int &quaindex,
    double * const &jac_value) const
{
  for( int ii=0; ii<9; ++ii ) jac_value[ii] = Jac[ii];
}

void FEAElement_MINI_P1_3D::get_invJacobian(const int &quaindex,
    double * const &jac_value) const
{
  for(int ii=0; ii<9; ++ii) jac_value[ii] = Jac[9+ii];
}

void FEAElement_MINI_P1_3D::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  ASSERT( quad -> get_dim() == 4, "FEAElement_MINI_P1_3D::buildBasis function error.\n");

  // Calculate the geometrical info first
  Jac[0] = - ctrl_x[0] + ctrl_x[1]; // dx_dr
  Jac[1] = - ctrl_x[0] + ctrl_x[2]; // dx_ds
  Jac[2] = - ctrl_x[0] + ctrl_x[3]; // dx_dt

  Jac[3] = - ctrl_y[0] + ctrl_y[1]; // dy_dr
  Jac[4] = - ctrl_y[0] + ctrl_y[2]; // dy_ds
  Jac[5] = - ctrl_y[0] + ctrl_y[3]; // dy_dt

  Jac[6] = - ctrl_z[0] + ctrl_z[1]; // dz_dr
  Jac[7] = - ctrl_z[0] + ctrl_z[2]; // dz_ds
  Jac[8] = - ctrl_z[0] + ctrl_z[3]; // dz_dt

  // Generate 3x3 matrix
  FE_T::Matrix_double_3by3_Array mdrdx( Jac[0], Jac[1], Jac[2], Jac[3], Jac[4],
      Jac[5], Jac[6], Jac[7], Jac[8] );

  detJac = mdrdx.det();

  mdrdx.inverse();

  // dr_dx
  Jac[9]  = mdrdx(0); Jac[10] = mdrdx(1); Jac[11] = mdrdx(2);
  
  // ds_dx
  Jac[12] = mdrdx(3); Jac[13] = mdrdx(4); Jac[14] = mdrdx(5);

  // dt_dx
  Jac[15] = mdrdx(6); Jac[16] = mdrdx(7); Jac[17] = mdrdx(8);

  dR_dx[0] = (-1.0) * Jac[9] - Jac[12] - Jac[15];
  dR_dx[1] = Jac[9];
  dR_dx[2] = Jac[12];
  dR_dx[3] = Jac[15];

  dR_dy[0] = (-1.0) * Jac[10] - Jac[13] - Jac[16];
  dR_dy[1] = Jac[10];
  dR_dy[2] = Jac[13];
  dR_dy[3] = Jac[16];

  dR_dz[0] = (-1.0) * Jac[11] - Jac[14] - Jac[17];
  dR_dz[1] = Jac[11];
  dR_dz[2] = Jac[14];
  dR_dz[3] = Jac[17];

  // area coordinates
  double qua_r, qua_s, qua_t, qua_u;
  for(int qua=0; qua < numQuapts; ++qua)
  {
    qua_r = quad -> get_qp( qua, 0 );
    qua_s = quad -> get_qp( qua, 1 );
    qua_t = quad -> get_qp( qua, 2 );
    qua_u = 1.0 - qua_r - qua_s - qua_t;

    R[qua*5+0] = qua_u;
    R[qua*5+1] = qua_r;
    R[qua*5+2] = qua_s;
    R[qua*5+3] = qua_t;
    R[qua*5+4] = 256.0 * qua_r * qua_s * qua_t * qua_u;
  
    const double N_r = 256.0 * qua_s * qua_t * (qua_u - qua_r);
    const double N_s = 256.0 * qua_r * qua_t * (qua_u - qua_s);
    const double N_t = 256.0 * qua_r * qua_s * (qua_u - qua_t);
 
    dB_dx[qua] = N_r * Jac[9] + N_s * Jac[12] + N_t * Jac[15];
    dB_dy[qua] = N_r * Jac[10] + N_s * Jac[13] + N_t * Jac[16];
    dB_dz[qua] = N_r * Jac[11] + N_s * Jac[14] + N_t * Jac[17];
  }
}

// EOF
