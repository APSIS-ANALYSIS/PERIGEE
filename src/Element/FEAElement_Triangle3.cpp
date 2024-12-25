#include "FEAElement_Triangle3.hpp"

FEAElement_Triangle3::FEAElement_Triangle3( const int &in_nqua )
: numQuapts( in_nqua )
{
  R.resize(3 * numQuapts);
}

void FEAElement_Triangle3::print_info() const
{
  SYS_T::commPrint("Tri3: ");
  SYS_T::commPrint("3-node triangle element with up to 2nd derivatives. \n");
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}

void FEAElement_Triangle3::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x, const double * const &ctrl_y )
{
  ASSERT(quad -> get_dim() == 3, "FEAElement_Triangle3::buildBasis function error.\n" );

  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );

    R[qua * 3 + 0] = 1.0 - qua_r - qua_s;
    R[qua * 3 + 1] = qua_r;
    R[qua * 3 + 2] = qua_s;
  }

  Jac[0] = - ctrl_x[0] + ctrl_x[1]; // dx_dr
  Jac[1] = - ctrl_x[0] + ctrl_x[2]; // dx_ds

  Jac[2] = - ctrl_y[0] + ctrl_y[1]; // dy_dr
  Jac[3] = - ctrl_y[0] + ctrl_y[2]; // dy_ds

  detJac = Jac[0] * Jac[3] - Jac[1] * Jac[2];

  double inv_detJac = 1.0 / detJac;

  Jac[4] =   Jac[3] * inv_detJac;  // dr_dx
  Jac[5] = - Jac[1] * inv_detJac;  // dr_dy
  Jac[6] = - Jac[2] * inv_detJac;  // ds_dx
  Jac[7] =   Jac[0] * inv_detJac;  // ds_dy

  dR_dx[0] = - Jac[4] - Jac[6];
  dR_dx[1] = Jac[4];
  dR_dx[2] = Jac[6];

  dR_dy[0] = - Jac[5] - Jac[7];
  dR_dy[1] = Jac[5];
  dR_dy[2] = Jac[7];
}

double FEAElement_Triangle3::get_h( const double * const &ctrl_x,
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

void FEAElement_Triangle3::get_R( const int &quaindex, 
    double * const &basis ) const
{
  const int offset = quaindex * 3;
  basis[0] = R[offset];
  basis[1] = R[offset+1];
  basis[2] = R[offset+2];
}

std::vector<double> FEAElement_Triangle3::get_R( const int &quaindex ) const
{
  const int offset = quaindex * 3;

  return { R[offset], R[offset+1], R[offset+2] };
}

void FEAElement_Triangle3::get_gradR( const int &quaindex, 
    double * const &basis_x, double * const &basis_y ) const
{
  for(int ii=0; ii<3; ++ii)
  {
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
  } 
}

void FEAElement_Triangle3::get_R_gradR( const int &quaindex, 
    double * const &basis,
    double * const &basis_x, double * const &basis_y ) const
{
  const int offset = quaindex * 3;
  for(int ii=0; ii < 3; ++ ii)
  {
    basis[ii] = R[offset + ii];
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
  }
}

void FEAElement_Triangle3::get_2D_R_dR_d2R( const int &quaindex,
    double * const &basis,
    double * const &basis_x, double * const &basis_y,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_xy ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle3::get_2D_R_dR_d2R function error.\n" );
  const int offset = quaindex * 3;
  for(int ii=0; ii < 3; ++ ii)
  {
    basis[ii] = R[offset + ii];
    basis_x[ii] = dR_dx[ii];
    basis_y[ii] = dR_dy[ii];
    basis_xx[ii] = 0.0;
    basis_xy[ii] = 0.0;
    basis_yy[ii] = 0.0;
  }
}

std::array<double,4> FEAElement_Triangle3::get_Jacobian_2D(const int &quaindex) const
{
  return {{ Jac[0], Jac[1], Jac[2], Jac[3] }};
}

std::array<double,4> FEAElement_Triangle3::get_invJacobian_2D(const int &quaindex) const
{
  return {{ Jac[4], Jac[5], Jac[6], Jac[7] }};
}

// EOF
