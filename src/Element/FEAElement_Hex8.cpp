#include "FEAElement_Hex8.hpp"

FEAElement_Hex8::FEAElement_Hex8( const int &in_nqua ) : numQuapts( in_nqua )
{
  R = new double [8 * numQuapts];

  dR_dx = new double [8 * numQuapts];
  dR_dy = new double [8 * numQuapts];
  dR_dz = new double [8 * numQuapts];

  d2R_dxx = new double [8 * numQuapts];
  d2R_dyy = new double [8 * numQuapts];
  d2R_dzz = new double [8 * numQuapts];
  d2R_dxy = new double [8 * numQuapts];
  d2R_dxz = new double [8 * numQuapts];
  d2R_dyz = new double [8 * numQuapts];

  dx_dr = new double [9*numQuapts];
  dr_dx = new double [9*numQuapts];
  detJac = new double [numQuapts];
}

FEAElement_Hex8::~FEAElement_Hex8()
{
  delete [] R;             R = nullptr;
  delete [] dR_dx;     dR_dx = nullptr;
  delete [] dR_dy;     dR_dy = nullptr;
  delete [] dR_dz;     dR_dz = nullptr;
  delete [] d2R_dxx; d2R_dxx = nullptr;
  delete [] d2R_dyy; d2R_dyy = nullptr;
  delete [] d2R_dzz; d2R_dzz = nullptr;
  delete [] d2R_dxy; d2R_dxy = nullptr;
  delete [] d2R_dxz; d2R_dxz = nullptr;
  delete [] d2R_dyz; d2R_dyz = nullptr;

  delete [] dx_dr;     dx_dr = nullptr;
  delete [] dr_dx;     dr_dx = nullptr;
  delete [] detJac;   detJac = nullptr;
}

void FEAElement_Hex8::print_info() const
{
  SYS_T::commPrint("Hex8: ");
  SYS_T::commPrint("8-node hexagon element with up to 2nd derivatives. \n");
  SYS_T::commPrint("elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}

double FEAElement_Hex8::get_memory_usage() const
{
  const double double_size = 99.0 * numQuapts;
  const double int_size = 1.0;
  return double_size * 8.0 + int_size * 4.0;
}

void FEAElement_Hex8::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  ASSERT( quad -> get_dim() == 3, "FEAElement_Hex8::buildBasis function error.\n" );

  for(int qua=0; qua<numQuapts; ++qua)
  {
    const int q8 = qua * 8;

    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );
    const double qua_t = quad -> get_qp( qua, 2 );
    
    R[q8  ] = (1.0 - qua_r) * (1.0 - qua_s) * (1.0 - qua_t);
    R[q8+1] = qua_r * (1.0 - qua_s) * (1.0 - qua_t);
    R[q8+2] = qua_r * qua_s * (1.0 - qua_t);
    R[q8+3] = (1.0 - qua_r) * qua_s * (1.0 - qua_t);
    R[q8+4] = (1.0 - qua_r) * (1.0 - qua_s) * qua_t;
    R[q8+5] = qua_r * (1.0 - qua_s) * qua_t;
    R[q8+6] = qua_r * qua_s * qua_t;
    R[q8+7] = (1.0 - qua_r) * qua_s * qua_t;

    const double dR_dr[8] { (qua_s - 1.0) * (1.0 - qua_t),
                            (1.0 - qua_s) * (1.0 - qua_t),
                             qua_s * (1.0 - qua_t),
                             qua_s * (qua_t - 1.0),
                            (qua_s - 1.0) * qua_t,
                            (1.0 - qua_s) * qua_t,
                             qua_s * qua_t,
                            -qua_s * qua_t };
    const double dR_ds[8] { (qua_r - 1.0) * (1.0 - qua_t),
                             qua_r * (qua_t - 1.0),
                             qua_r * (1.0 - qua_t),
                            (1.0 - qua_r) * (1.0 - qua_t),
                            (qua_r - 1.0) * qua_t,
                            -qua_r * qua_t,
                             qua_r * qua_t,
                            (1.0 - qua_r) * qua_t };
    const double dR_dt[8] { (qua_s - 1.0) * (1.0 - qua_r),
                            (qua_s - 1.0) * qua_r,
                            -qua_s * qua_r,
                             qua_s * (qua_r - 1.0),
                            (1.0 - qua_s) * (1.0 - qua_r),
                            (1.0 - qua_s) * qua_r,
                             qua_s * qua_r,
                             qua_s * (1.0 - qua_r) };
    
    const double d2R_drs[8] = { 1.0 - qua_t, qua_t - 1.0, 1.0 - qua_t, qua_t - 1.0, 
                                qua_t, -qua_t, qua_t, -qua_t };
    const double d2R_drt[8] = { 1.0 - qua_s, qua_s - 1.0, -qua_s, qua_s,
                                qua_s - 1.0, 1.0 - qua_s, qua_s, -qua_s};
    const double d2R_dst[8] = { 1.0 - qua_r, qua_r, -qua_r, qua_r - 1.0,
                                qua_r - 1.0, -qua_r, qua_r, 1.0 - qua_r};
    
    double xr = 0.0, xs = 0.0, xt = 0.0;
    double yr = 0.0, ys = 0.0, yt = 0.0;
    double zr = 0.0, zs = 0.0, zt = 0.0;
    
    double xrs = 0.0, xrt = 0.0, xst = 0.0;
    double yrs = 0.0, yrt = 0.0, yst = 0.0;
    double zrs = 0.0, zrt = 0.0, zst = 0.0;

    for(int ii=0; ii<8; ++ii)
    {
      xr += ctrl_x[ii] * dR_dr[ii];
      xs += ctrl_x[ii] * dR_ds[ii];
      xt += ctrl_x[ii] * dR_dt[ii];
      
      yr += ctrl_y[ii] * dR_dr[ii];
      ys += ctrl_y[ii] * dR_ds[ii];
      yt += ctrl_y[ii] * dR_dt[ii];
      
      zr += ctrl_z[ii] * dR_dr[ii];
      zs += ctrl_z[ii] * dR_ds[ii];
      zt += ctrl_z[ii] * dR_dt[ii];

      xrs += ctrl_x[ii] * d2R_drs[ii];
      xrt += ctrl_x[ii] * d2R_drt[ii];
      xst += ctrl_x[ii] * d2R_dst[ii];
      
      yrs += ctrl_y[ii] * d2R_drs[ii];
      yrt += ctrl_y[ii] * d2R_drt[ii];
      yst += ctrl_y[ii] * d2R_dst[ii];
      
      zrs += ctrl_z[ii] * d2R_drs[ii];
      zrt += ctrl_z[ii] * d2R_drt[ii];
      zst += ctrl_z[ii] * d2R_dst[ii];
    }

    FE_T::Matrix_double_3by3_Array mdrdx(xr, xs, xt, yr, ys, yt, zr, zs, zt);

    detJac[qua] = mdrdx.det(); // detJac = |dx/dr|
 
    mdrdx.inverse();

    const int q9 = qua * 9;
    dx_dr[q9  ] = xr; dx_dr[q9+1] = xs; dx_dr[q9+2] = xt;
    dx_dr[q9+3] = yr; dx_dr[q9+4] = ys; dx_dr[q9+5] = yt;
    dx_dr[q9+6] = zr; dx_dr[q9+7] = zs; dx_dr[q9+8] = zt;

    dr_dx[q9  ] = mdrdx(0); // dr_dx 
    dr_dx[q9+1] = mdrdx(1); // dr_dy
    dr_dx[q9+2] = mdrdx(2); // dr_dz
    dr_dx[q9+3] = mdrdx(3); // ds_dx
    dr_dx[q9+4] = mdrdx(4); // ds_dy
    dr_dx[q9+5] = mdrdx(5); // ds_dz
    dr_dx[q9+6] = mdrdx(6); // dt_dx
    dr_dx[q9+7] = mdrdx(7); // dt_dy
    dr_dx[q9+8] = mdrdx(8); // dt_dz

    for(int ii=0; ii<8; ++ii)
    {
      dR_dx[q8+ii] = dR_dr[ii]*mdrdx(0) + dR_ds[ii]*mdrdx(3) + dR_dt[ii]*mdrdx(6); 
      dR_dy[q8+ii] = dR_dr[ii]*mdrdx(1) + dR_ds[ii]*mdrdx(4) + dR_dt[ii]*mdrdx(7); 
      dR_dz[q8+ii] = dR_dr[ii]*mdrdx(2) + dR_ds[ii]*mdrdx(5) + dR_dt[ii]*mdrdx(8); 
    }

    // Setup the 6x6 matrix
    FE_T::Matrix_double_6by6_Array LHS(xr, xs, xt, yr, ys, yt, zr, zs, zt);

    // LU factorization
    LHS.LU_fac();

    for(int ii=0; ii<8; ++ii)
    {
      const std::array<double, 6> RHS {{ 0.0,
        d2R_drs[ii] - dR_dx[q8+ii] * xrs - dR_dy[q8+ii] * yrs - dR_dz[q8+ii] * zrs,
        d2R_drt[ii] - dR_dx[q8+ii] * xrt - dR_dy[q8+ii] * yrt - dR_dz[q8+ii] * zrt,
        0.0,
        d2R_dst[ii] - dR_dx[q8+ii] * xst - dR_dy[q8+ii] * yst - dR_dz[q8+ii] * zst,
        0.0 }};

      const auto sol = LHS.LU_solve(RHS);

      d2R_dxx[q8+ii] = sol[0];
      d2R_dyy[q8+ii] = sol[1];
      d2R_dzz[q8+ii] = sol[2];
      d2R_dxy[q8+ii] = sol[3];
      d2R_dxz[q8+ii] = sol[4];
      d2R_dyz[q8+ii] = sol[5];
    }
  }
}

double FEAElement_Hex8::get_h( const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z ) const
{
  const double diag[4] { std::pow((ctrl_x[0] - ctrl_x[6]), 2.0)
      + std::pow((ctrl_y[0] - ctrl_y[6]), 2.0) + std::pow((ctrl_z[0] - ctrl_z[6]), 2.0),
      std::pow((ctrl_x[1] - ctrl_x[7]), 2.0)
      + std::pow((ctrl_y[1] - ctrl_y[7]), 2.0) + std::pow((ctrl_z[1] - ctrl_z[7]), 2.0),
      std::pow((ctrl_x[2] - ctrl_x[4]), 2.0)
      + std::pow((ctrl_y[2] - ctrl_y[4]), 2.0) + std::pow((ctrl_z[2] - ctrl_z[4]), 2.0),
      std::pow((ctrl_x[3] - ctrl_x[5]), 2.0)
      + std::pow((ctrl_y[3] - ctrl_y[5]), 2.0) + std::pow((ctrl_z[3] - ctrl_z[5]), 2.0) };
    
  double d = diag[0];
  for(int ii = 1; ii<4; ++ii)
  {
    if(diag[ii] > d) d = diag[ii];
  }

  return std::sqrt(d);
}

void FEAElement_Hex8::get_R( const int &quaindex, double * const &basis ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_R function error.\n" );
  const int offset = quaindex * 8;
  for(int ii=0; ii<8; ++ii) basis[ii] = R[offset+ii];
}

std::vector<double> FEAElement_Hex8::get_R( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_R function error.\n" );
  const int offset = quaindex * 8;
  return { R[offset], R[offset+1], R[offset+2], R[offset+3],
         R[offset+4], R[offset+5], R[offset+6], R[offset+7] };
}

void FEAElement_Hex8::get_gradR( const int &quaindex, double * const &basis_x,
    double * const &basis_y, double * const &basis_z ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_gradR function error.\n" );
  const int offset = quaindex * 8;
  for( int ii=0; ii<8; ++ii )
  {
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_z[ii] = dR_dz[offset + ii];
  }
}

void FEAElement_Hex8::get_R_gradR( const int &quaindex, double * const &basis,
    double * const &basis_x, double * const &basis_y,
    double * const &basis_z ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_R_gradR function error.\n" );
  const int offset = quaindex * 8;
  for( int ii=0; ii<8; ++ii )
  {
    basis[ii]   = R[offset+ ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_z[ii] = dR_dz[offset + ii];
  }
}

std::vector<double> FEAElement_Hex8::get_dR_dx( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_dR_dx function error.\n" );
  const int offset = quaindex * 8;
  return { dR_dx[offset], dR_dx[offset+1], dR_dx[offset+2], dR_dx[offset+3],
    dR_dx[offset+4], dR_dx[offset+5], dR_dx[offset+6], dR_dx[offset+7] };
}

std::vector<double> FEAElement_Hex8::get_dR_dy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_dR_dy function error.\n" );
  const int offset = quaindex * 8;
  return { dR_dy[offset], dR_dy[offset+1], dR_dy[offset+2], dR_dy[offset+3],
    dR_dy[offset+4], dR_dy[offset+5], dR_dy[offset+6], dR_dy[offset+7] };
}

std::vector<double> FEAElement_Hex8::get_dR_dz( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_dR_dz function error.\n" );
  const int offset = quaindex * 8;
  return { dR_dz[offset], dR_dz[offset+1], dR_dz[offset+2], dR_dz[offset+3],
    dR_dz[offset+4], dR_dz[offset+5], dR_dz[offset+6], dR_dz[offset+7] };
}

void FEAElement_Hex8::get_3D_R_dR_d2R( const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y, double * const &basis_z,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_zz, double * const &basis_xy,
    double * const &basis_xz, double * const &basis_yz ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_3D_R_dR_d2R function error.\n" );
  const int offset = quaindex * 8;
  for( int ii=0; ii<8; ++ii )
  {
    basis[ii]   = R[offset + ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_z[ii] = dR_dz[offset + ii];
    basis_xx[ii] = d2R_dxx[offset + ii];
    basis_yy[ii] = d2R_dyy[offset + ii];
    basis_zz[ii] = d2R_dzz[offset + ii];
    basis_xy[ii] = d2R_dxy[offset + ii];
    basis_xz[ii] = d2R_dxz[offset + ii];
    basis_yz[ii] = d2R_dyz[offset + ii];
  }
}

void FEAElement_Hex8::get_3D_R_gradR_LaplacianR( const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y, double * const &basis_z,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_zz ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_3D_R_gradR_LaplacianR function error.\n" );
  const int offset = quaindex * 8;
  for( int ii=0; ii<8; ++ii )
  {
    basis[ii]   = R[offset + ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_z[ii] = dR_dz[offset + ii];
    basis_xx[ii] = d2R_dxx[offset + ii];
    basis_yy[ii] = d2R_dyy[offset + ii];
    basis_zz[ii] = d2R_dzz[offset + ii];
  }
}

std::vector<double> FEAElement_Hex8::get_d2R_dxx( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_d2R_dxx function error.\n" );
  const int offset = quaindex * 8;
  return { d2R_dxx[offset], d2R_dxx[offset+1], d2R_dxx[offset+2], d2R_dxx[offset+3],
    d2R_dxx[offset+4], d2R_dxx[offset+5], d2R_dxx[offset+6], d2R_dxx[offset+7] };
}

std::vector<double> FEAElement_Hex8::get_d2R_dyy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_d2R_dyy function error.\n" );
  const int offset = quaindex * 8;
  return { d2R_dyy[offset], d2R_dyy[offset+1], d2R_dyy[offset+2], d2R_dyy[offset+3],
    d2R_dyy[offset+4], d2R_dyy[offset+5], d2R_dyy[offset+6], d2R_dyy[offset+7] };
}

std::vector<double> FEAElement_Hex8::get_d2R_dzz( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_d2R_dzz function error.\n" );
  const int offset = quaindex * 8;
  return { d2R_dzz[offset], d2R_dzz[offset+1], d2R_dzz[offset+2], d2R_dzz[offset+3],
    d2R_dzz[offset+4], d2R_dzz[offset+5], d2R_dzz[offset+6], d2R_dzz[offset+7] };
}

std::vector<double> FEAElement_Hex8::get_d2R_dxy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_d2R_dxy function error.\n" );
  const int offset = quaindex * 8;
  return { d2R_dxy[offset], d2R_dxy[offset+1], d2R_dxy[offset+2], d2R_dxy[offset+3],
    d2R_dxy[offset+4], d2R_dxy[offset+5], d2R_dxy[offset+6], d2R_dxy[offset+7] };
}

std::vector<double> FEAElement_Hex8::get_d2R_dxz( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_d2R_dxz function error.\n" );
  const int offset = quaindex * 8;
  return { d2R_dxz[offset], d2R_dxz[offset+1], d2R_dxz[offset+2], d2R_dxz[offset+3],
    d2R_dxz[offset+4], d2R_dxz[offset+5], d2R_dxz[offset+6], d2R_dxz[offset+7] };
}

std::vector<double> FEAElement_Hex8::get_d2R_dyz( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_d2R_dyz function error.\n" );
  const int offset = quaindex * 8;
  return { d2R_dyz[offset], d2R_dyz[offset+1], d2R_dyz[offset+2], d2R_dyz[offset+3],
    d2R_dyz[offset+4], d2R_dyz[offset+5], d2R_dyz[offset+6], d2R_dyz[offset+7] };
}

void FEAElement_Hex8::get_Jacobian(const int &quaindex,
    double * const &jac_value) const
{
  for(int ii=0; ii<9; ++ii) jac_value[ii] = dx_dr[9*quaindex + ii];
}

std::array<double,9> FEAElement_Hex8::get_Jacobian(const int &quaindex) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_Jacobian function error.\n" );
  return {{ dx_dr[9*quaindex], dx_dr[9*quaindex+1], dx_dr[9*quaindex+2],
    dx_dr[9*quaindex+3], dx_dr[9*quaindex+4], dx_dr[9*quaindex+5],
    dx_dr[9*quaindex+6], dx_dr[9*quaindex+7], dx_dr[9*quaindex+8] }};
}

void FEAElement_Hex8::get_invJacobian(const int &quaindex,
    double * const &jac_value) const
{
  for(int ii=0; ii<9; ++ii) jac_value[ii] = dr_dx[9*quaindex + ii];
}

std::array<double,9> FEAElement_Hex8::get_invJacobian(const int &quaindex) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex8::get_invJacobian function error.\n" );
  return {{ dr_dx[9*quaindex], dr_dx[9*quaindex+1], dr_dx[9*quaindex+2],
    dr_dx[9*quaindex+3], dr_dx[9*quaindex+4], dr_dx[9*quaindex+5],
    dr_dx[9*quaindex+6], dr_dx[9*quaindex+7], dr_dx[9*quaindex+8] }};
}

// EOF
