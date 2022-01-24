#include "FEAElement_MINI_P1_2D.hpp"

FEAElement_MINI_P1_2D::FEAElement_MINI_P1_2D( const int &in_nqp )
: numQuapts( in_nqp )
{
  R = new double [4*numQuapts];
  dB_dx = new double [numQuapts];
  dB_dy = new double [numQuapts];
}


FEAElement_MINI_P1_2D::~FEAElement_MINI_P1_2D()
{
  delete [] R;     R     = NULL;
  delete [] dB_dx; dB_dx = NULL;
  delete [] dB_dy; dB_dy = NULL;
}

void FEAElement_MINI_P1_2D::print_info() const
{
  SYS_T::commPrint("MINI_P1: ");
  SYS_T::commPrint("4-node triangle element with cubic bubble. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}


double FEAElement_MINI_P1_2D::get_memory_usage() const
{
  const double double_size = 6 * numQuapts + 15;
  const double int_size = 1;
  return double_size * 8.0 + int_size * 4.0;
}


void FEAElement_MINI_P1_2D::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x, const double * const &ctrl_y )
{
  assert( quad -> get_dim() == 3 );
 
  // Because of the constant strain, we can calculate the 
  // geometrical info first 
  Jac[0] = ctrl_x[0] * (-1.0) + ctrl_x[1]; // dx_dr
  Jac[1] = ctrl_x[0] * (-1.0) + ctrl_x[2]; // dx_ds
  Jac[2] = ctrl_y[0] * (-1.0) + ctrl_y[1]; // dy_dr
  Jac[3] = ctrl_y[0] * (-1.0) + ctrl_y[2]; // dy_ds

  detJac = Jac[0] * Jac[3] - Jac[1] * Jac[2];

  const double inv_detJac = 1.0 / detJac;

  Jac[4] = Jac[3] * inv_detJac;        // dr_dx
  Jac[5] = -1.0 * Jac[1] * inv_detJac; // dr_dy
  Jac[6] = -1.0 * Jac[2] * inv_detJac; // ds_dx
  Jac[7] = Jac[0] * inv_detJac;        // ds_dy

  // dR_dx = dR_dr * dr_dx + dR_ds * ds_dx
  dR_dx[0] = (-1.0) * Jac[4] - Jac[6];
  dR_dx[1] = Jac[4];
  dR_dx[2] = Jac[6];

  // dR_dy = dR_dr * dr_dy + dR_ds * ds_dy
  dR_dy[0] = (-1.0) * Jac[5] - Jac[7];
  dR_dy[1] = Jac[5];
  dR_dy[2] = Jac[7];
  
  double qua_r, qua_s;
  for( int qua = 0; qua < numQuapts; ++qua )
  {
    qua_r = quad -> get_qp( qua, 0 );
    qua_s = quad -> get_qp( qua, 1 );
    
    R[qua*4 + 0] = 1.0 - qua_r - qua_s;
    R[qua*4 + 1] = qua_r;
    R[qua*4 + 2] = qua_s;
    R[qua*4 + 3] = 27.0 * qua_r * qua_s * (1.0 - qua_r - qua_s);
  
    const double N3_r = 27.0 * (qua_s - 2.0 * qua_r * qua_s - qua_s * qua_s);
    const double N3_s = 27.0 * (qua_r - 2.0 * qua_r * qua_s - qua_r * qua_r);
    dB_dx[qua] = N3_r * Jac[4] + N3_s * Jac[6];
    dB_dy[qua] = N3_r * Jac[5] + N3_s * Jac[7];
  }
}


double FEAElement_MINI_P1_2D::get_h( const double * const &ctrl_x,
    const double * const &ctrl_y ) const
{
  const double a = 2.0 * ( ctrl_x[1] - ctrl_x[0] );
  const double b = 2.0 * ( ctrl_y[1] - ctrl_y[0] );
  const double c = 2.0 * ( ctrl_x[2] - ctrl_x[0] );
  const double d = 2.0 * ( ctrl_y[2] - ctrl_y[0] );

  const double m = ctrl_x[1] * ctrl_x[1] + ctrl_y[1] * ctrl_y[1]
    - ctrl_x[0] * ctrl_x[0] - ctrl_y[0] * ctrl_y[0];

  const double n = ctrl_x[2] * ctrl_x[2] + ctrl_y[2] * ctrl_y[2]
    - ctrl_x[0] * ctrl_x[0] - ctrl_y[0] * ctrl_y[0];

  const double xc = (d * m - b * n) / (a*d-b*c);
  const double yc = (a * n - c * m) / (a*d-b*c);

  const double radius = std::sqrt( (xc-ctrl_x[0])*(xc-ctrl_x[0])
      + (yc - ctrl_y[0])*(yc - ctrl_y[0]) );

  return 2.0 * radius;
}


void FEAElement_MINI_P1_2D::get_R( const int &quaindex, 
    double * const &basis ) const
{
  const int offset = quaindex * 4;
  basis[0] = R[offset];
  basis[1] = R[offset + 1];
  basis[2] = R[offset + 2];
  basis[3] = R[offset + 3];
}


void FEAElement_MINI_P1_2D::get_gradR( const int &quaindex, 
    double * const &basis_x, double * const &basis_y ) const
{
  basis_x[0] = dR_dx[0];
  basis_x[1] = dR_dx[1];
  basis_x[2] = dR_dx[2];
  basis_x[3] = dB_dx[quaindex];

  basis_y[0] = dR_dy[0];
  basis_y[1] = dR_dy[1];
  basis_y[2] = dR_dy[2];
  basis_y[3] = dB_dy[quaindex];
}


void FEAElement_MINI_P1_2D::get_R_gradR( const int &quaindex, 
    double * const &basis,
    double * const &basis_x, double * const &basis_y ) const
{
  const int offset = quaindex * 4;
  basis[0] = R[offset];
  basis[1] = R[offset + 1];
  basis[2] = R[offset + 2];
  basis[3] = R[offset + 3];

  basis_x[0] = dR_dx[0];
  basis_x[1] = dR_dx[1];
  basis_x[2] = dR_dx[2];
  basis_x[3] = dB_dx[quaindex];

  basis_y[0] = dR_dy[0];
  basis_y[1] = dR_dy[1];
  basis_y[2] = dR_dy[2];
  basis_y[3] = dB_dy[quaindex];
}


void FEAElement_MINI_P1_2D::get_Jacobian(const int &quaindex,
    double * const &jac_value) const
{
  jac_value[0] = Jac[0];
  jac_value[1] = Jac[1];
  jac_value[2] = Jac[2];
  jac_value[3] = Jac[3];
}


void FEAElement_MINI_P1_2D::get_invJacobian(const int &quaindex,
    double * const &jac_value) const
{
  jac_value[0] = Jac[4];
  jac_value[1] = Jac[5];
  jac_value[2] = Jac[6];
  jac_value[3] = Jac[7];
}


// EOF
