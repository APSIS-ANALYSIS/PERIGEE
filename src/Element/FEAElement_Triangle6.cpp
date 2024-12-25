#include "FEAElement_Triangle6.hpp"

FEAElement_Triangle6::FEAElement_Triangle6( const int &in_nqua )
: numQuapts( in_nqua )
{
  R.resize(6 * numQuapts);

  dR_dx.resize(6 * numQuapts);
  dR_dy.resize(6 * numQuapts);

  d2R_dxx.resize(6 * numQuapts);
  d2R_dyy.resize(6 * numQuapts);
  d2R_dxy.resize(6 * numQuapts);
  
  Jac.resize(9 * numQuapts);
}

void FEAElement_Triangle6::print_info() const
{
  SYS_T::commPrint("Tri6: ");
  SYS_T::commPrint("6-node triangle element with up to 2nd derivatives. \n");
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}

void FEAElement_Triangle6::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y )
{
  ASSERT(quad->get_dim() == 3, "FEAElement_Triangle6::buildBasis function error.\n" );

  const double d2R_drr[6] { 4.0, 4.0, 0.0, -8.0, 0.0,  0.0 };
  const double d2R_dss[6] { 4.0, 0.0, 4.0,  0.0, 0.0, -8.0 };
  const double d2R_drs[6] { 4.0, 0.0, 0.0, -4.0, 4.0, -4.0 };

  // Caclulate second derivative of geometry
  // Here, second derivatives d2R_drr, etc are constant. We can calculate
  // xrr, etc. out of the quadrature loop.
  double xrr = 0.0, xss = 0.0, xrs = 0.0;
  double yrr = 0.0, yss = 0.0, yrs = 0.0;

  for(int ii=0; ii<6; ++ii)
  {
    xrr += ctrl_x[ii] * d2R_drr[ii];
    xss += ctrl_x[ii] * d2R_dss[ii];
    xrs += ctrl_x[ii] * d2R_drs[ii];

    yrr += ctrl_y[ii] * d2R_drr[ii];
    yss += ctrl_y[ii] * d2R_dss[ii];
    yrs += ctrl_y[ii] * d2R_drs[ii];
  }

  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );
    const double qua_t = 1.0 - qua_r - qua_s;

    const int offset = 6 * qua;

    R[offset + 0] = qua_t * (2.0 * qua_t - 1.0);
    R[offset + 1] = qua_r * (2.0 * qua_r - 1.0);
    R[offset + 2] = qua_s * (2.0 * qua_s - 1.0);
    R[offset + 3] = 4.0 * qua_r * qua_t;
    R[offset + 4] = 4.0 * qua_r * qua_s;
    R[offset + 5] = 4.0 * qua_s * qua_t;

    const double dR_dr[6] { 4.0 * qua_r + 4.0 * qua_s - 3.0,
      4.0 * qua_r - 1.0, 0.0, 4.0 - 8.0 * qua_r - 4.0 * qua_s,
      4.0 * qua_s, -4.0 * qua_s };
    
    const double dR_ds[6] { 4.0 * qua_r + 4.0 * qua_s - 3.0,
      0.0, 4.0 * qua_s - 1.0, -4.0 * qua_r, 4.0 * qua_r,
      4.0 - 4.0 * qua_r - 8.0 * qua_s };  
    
    double dx_dr = 0.0, dx_ds = 0.0, dy_dr = 0.0, dy_ds = 0.0;
    for(int ii=0; ii<6; ++ii)
    {
      dx_dr += ctrl_x[ii] * dR_dr[ii];
      dx_ds += ctrl_x[ii] * dR_ds[ii];
      dy_dr += ctrl_y[ii] * dR_dr[ii];
      dy_ds += ctrl_y[ii] * dR_ds[ii];
    }
    
    Jac[4*qua+0] = dx_dr;
    Jac[4*qua+1] = dx_ds;
    Jac[4*qua+2] = dy_dr;
    Jac[4*qua+3] = dy_ds;

    Jac[8*numQuapts + qua] = dx_dr * dy_ds - dx_ds * dy_dr;

    const double inv_detJac = 1.0 / Jac[8*numQuapts + qua];

    const double dr_dx = dy_ds * inv_detJac;
    const double dr_dy = (-1.0) * dx_ds * inv_detJac;
    const double ds_dx = (-1.0) * dy_dr * inv_detJac;
    const double ds_dy = dx_dr * inv_detJac;

    Jac[4*numQuapts + 4*qua + 0] = dr_dx;
    Jac[4*numQuapts + 4*qua + 1] = dr_dy;
    Jac[4*numQuapts + 4*qua + 2] = ds_dx;
    Jac[4*numQuapts + 4*qua + 3] = ds_dy;
  
    for(int ii=0; ii<6; ++ii)
    {
      dR_dx[offset+ii] = dR_dr[ii] * dr_dx + dR_ds[ii] * ds_dx;
      dR_dy[offset+ii] = dR_dr[ii] * dr_dy + dR_ds[ii] * ds_dy;
    }

    // Setup the 3x3 matrix
    const double a11 = dx_dr * dx_dr;
    const double a12 = dy_dr * dy_dr;
    const double a13 = 2 * dx_dr * dy_dr;
    const double a21 = dx_dr * dx_ds;
    const double a22 = dy_dr * dy_ds;
    const double a23 = dx_dr * dy_ds + dx_ds * dy_dr;
    const double a31 = dx_ds * dx_ds;
    const double a32 = dy_ds * dy_ds;
    const double a33 = 2 * dx_ds * dy_ds;

    FE_T::Matrix_double_3by3_Array LHS(a11, a12, a13, a21, a22, a23, a31, a32, a33);

    // LU factorization
    LHS.LU_fac();

    for(int ii=0; ii<6; ++ii)
    {
      const std::array<double, 3> RHS {{ d2R_drr[ii] - dR_dx[offset+ii] * xrr - dR_dy[offset+ii] * yrr,
        d2R_drs[ii] - dR_dx[offset+ii] * xrs - dR_dy[offset+ii] * yrs,
        d2R_dss[ii] - dR_dx[offset+ii] * xss - dR_dy[offset+ii] * yss }};

      const auto sol = LHS.LU_solve(RHS);

      d2R_dxx[offset+ii] = sol[0];
      d2R_dyy[offset+ii] = sol[1];
      d2R_dxy[offset+ii] = sol[2];
    }
  }
}

double FEAElement_Triangle6::get_h( const double * const &ctrl_x,
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

void FEAElement_Triangle6::get_R( const int &quaindex, 
    double * const &basis ) const
{
  const int offset = quaindex * 6;
  basis[0] = R[offset];   basis[1] = R[offset+1];
  basis[2] = R[offset+2]; basis[3] = R[offset+3];
  basis[4] = R[offset+4]; basis[5] = R[offset+5];
}

std::vector<double> FEAElement_Triangle6::get_R( const int &quaindex ) const
{
  const int offset = quaindex * 6;
  return { R[offset], R[offset+1], R[offset+2], R[offset+3], R[offset+4], R[offset+5] };
}

void FEAElement_Triangle6::get_gradR( const int &quaindex, 
    double * const &basis_x, double * const &basis_y ) const
{
  const int offset = quaindex * 6;
  for(int ii=0; ii<6; ++ii)
  {
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
  }
}

void FEAElement_Triangle6::get_R_gradR( const int &quaindex, 
    double * const &basis, double * const &basis_x, 
    double * const &basis_y ) const
{
  const int offset = quaindex * 6;
  for(int ii=0; ii<6; ++ii)
  {
    basis[ii]   = R[offset + ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
  }
}

void FEAElement_Triangle6::get_2D_R_dR_d2R( const int &quaindex,
    double * const &basis,
    double * const &basis_x, double * const &basis_y,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_xy ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle6::get_2D_R_dR_d2R function error.\n" );
  const int offset = quaindex * 6;
  for( int ii=0; ii<6; ++ii )
  {
    basis[ii]   = R[offset + ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_xx[ii] = d2R_dxx[offset + ii];
    basis_yy[ii] = d2R_dyy[offset + ii];
    basis_xy[ii] = d2R_dxy[offset + ii];
  }
}

std::array<double,4> FEAElement_Triangle6::get_Jacobian_2D(const int &quaindex) const
{
  return {{ Jac[4*quaindex], Jac[4*quaindex+1],
    Jac[4*quaindex+2], Jac[4*quaindex+3] }};
}

std::array<double,4> FEAElement_Triangle6::get_invJacobian_2D(const int &quaindex) const
{
  const int offset = 4 * numQuapts + 4 * quaindex;
  return {{ Jac[offset], Jac[offset+1], Jac[offset+2], Jac[offset+3] }};
}

// EOF
