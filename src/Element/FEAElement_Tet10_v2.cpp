#include "FEAElement_Tet10_v2.hpp"

FEAElement_Tet10_v2::FEAElement_Tet10_v2( const int &in_nqua )
: numQuapts( in_nqua )
{
  R = new double [10 * numQuapts];

  dR_dx = new double [10 * numQuapts];
  dR_dy = new double [10 * numQuapts];
  dR_dz = new double [10 * numQuapts];

  d2R_dxx = new double [10 * numQuapts];
  d2R_dyy = new double [10 * numQuapts];
  d2R_dzz = new double [10 * numQuapts];
  d2R_dxy = new double [10 * numQuapts];
  d2R_dxz = new double [10 * numQuapts];
  d2R_dyz = new double [10 * numQuapts];

  dx_dr = new double [9*numQuapts];
  dr_dx = new double [9*numQuapts];
  detJac = new double [numQuapts];

  // second derivative w.r.t r-s-t is constant
  for(int ii=0; ii<10; ++ii)
  {
    d2R_drr[ii] = 0.0; d2R_dss[ii] = 0.0; d2R_dtt[ii] = 0.0;
    d2R_drs[ii] = 0.0; d2R_drt[ii] = 0.0; d2R_dst[ii] = 0.0;
  }

  d2R_drr[0] = 4.0; d2R_drr[1] = 4.0; d2R_drr[4] = -8.0;
  d2R_dss[0] = 4.0; d2R_dss[2] = 4.0; d2R_dss[6] = -8.0;
  d2R_dtt[0] = 4.0; d2R_dtt[3] = 4.0; d2R_dtt[7] = -8.0;
  d2R_drs[0] = 4.0; d2R_drs[4] = -4.0; d2R_drs[5] = 4.0; d2R_drs[6] = -4.0;
  d2R_drt[0] = 4.0; d2R_drt[4] = -4.0; d2R_drt[7] = -4.0; d2R_drt[8] = 4.0;
  d2R_dst[0] = 4.0; d2R_dst[6] = -4.0; d2R_dst[7] = -4.0; d2R_dst[9] = 4.0;
}

FEAElement_Tet10_v2::~FEAElement_Tet10_v2()
{
  delete [] R; R = NULL;
  delete [] dR_dx; dR_dx = NULL;
  delete [] dR_dy; dR_dy = NULL;
  delete [] dR_dz; dR_dz = NULL;
  delete [] d2R_dxx; d2R_dxx = NULL;
  delete [] d2R_dyy; d2R_dyy = NULL;
  delete [] d2R_dzz; d2R_dzz = NULL;
  delete [] d2R_dxy; d2R_dxy = NULL;
  delete [] d2R_dxz; d2R_dxz = NULL;
  delete [] d2R_dyz; d2R_dyz = NULL;

  delete [] dx_dr; dx_dr = NULL;
  delete [] dr_dx; dr_dx = NULL;
  delete [] detJac; detJac = NULL;
}

void FEAElement_Tet10_v2::print_info() const
{
  SYS_T::commPrint("Tet10_v2: ");
  SYS_T::commPrint("10-node tetrahedral element with up to 2nd derivatives. Compatible with vtk format \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}

double FEAElement_Tet10_v2::get_memory_usage() const
{
  const double d_size = 119 * numQuapts + 90;
  const double i_size = 1;
  return d_size * 8.0 + i_size * 4.0;
}

void FEAElement_Tet10_v2::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  assert( quad -> get_dim() == 4 );

  // Caclulate second derivative of geometry
  // Here, second derivatives d2R_drr, etc are constant. We can calculate
  // xrr, etc. out of the quadrature loop.
  double xrr = 0.0, xss = 0.0, xtt = 0.0, xrs = 0.0, xrt = 0.0, xst = 0.0;
  double yrr = 0.0, yss = 0.0, ytt = 0.0, yrs = 0.0, yrt = 0.0, yst = 0.0;
  double zrr = 0.0, zss = 0.0, ztt = 0.0, zrs = 0.0, zrt = 0.0, zst = 0.0;

  for(int ii=0; ii<10; ++ii)
  {
    xrr += ctrl_x[ii] * d2R_drr[ii];
    xss += ctrl_x[ii] * d2R_dss[ii];
    xtt += ctrl_x[ii] * d2R_dtt[ii];
    xrs += ctrl_x[ii] * d2R_drs[ii];
    xrt += ctrl_x[ii] * d2R_drt[ii];
    xst += ctrl_x[ii] * d2R_dst[ii];
    
    yrr += ctrl_y[ii] * d2R_drr[ii];
    yss += ctrl_y[ii] * d2R_dss[ii];
    ytt += ctrl_y[ii] * d2R_dtt[ii];
    yrs += ctrl_y[ii] * d2R_drs[ii];
    yrt += ctrl_y[ii] * d2R_drt[ii];
    yst += ctrl_y[ii] * d2R_dst[ii];
    
    zrr += ctrl_z[ii] * d2R_drr[ii];
    zss += ctrl_z[ii] * d2R_dss[ii];
    ztt += ctrl_z[ii] * d2R_dtt[ii];
    zrs += ctrl_z[ii] * d2R_drs[ii];
    zrt += ctrl_z[ii] * d2R_drt[ii];
    zst += ctrl_z[ii] * d2R_dst[ii];
  }

  for(int qua=0; qua<numQuapts; ++qua)
  {
    const int q10 = qua * 10;

    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );
    const double qua_t = quad -> get_qp( qua, 2 );
    const double qua_u = quad -> get_qp( qua, 3 ); 

    R[q10+0] = qua_u * (2.0*qua_u - 1.0);
    R[q10+1] = qua_r * (2.0*qua_r - 1.0);
    R[q10+2] = qua_s * (2.0*qua_s - 1.0);
    R[q10+3] = qua_t * (2.0*qua_t - 1.0);
    R[q10+4] = 4.0 * qua_u * qua_r;
    R[q10+5] = 4.0 * qua_r * qua_s;
    R[q10+6] = 4.0 * qua_s * qua_u;
    R[q10+7] = 4.0 * qua_t * qua_u;
    R[q10+8] = 4.0 * qua_r * qua_t;
    R[q10+9] = 4.0 * qua_s * qua_t;

    dR_dr[0] = 1.0 - 4.0 * qua_u;
    dR_dr[1] = 4.0 * qua_r - 1.0;
    dR_dr[2] = 0.0;
    dR_dr[3] = 0.0;
    dR_dr[4] = 4.0 * (qua_u - qua_r);
    dR_dr[5] = 4.0 * qua_s;
    dR_dr[6] = -4.0 * qua_s;
    dR_dr[7] = -4.0 * qua_t;
    dR_dr[8] = 4.0 * qua_t;
    dR_dr[9] = 0.0; 

    dR_ds[0] = 1.0 - 4.0 * qua_u;
    dR_ds[1] = 0.0;
    dR_ds[2] = 4.0 * qua_s - 1.0;
    dR_ds[3] = 0.0;
    dR_ds[4] = -4.0 * qua_r;
    dR_ds[5] = 4.0 * qua_r;
    dR_ds[6] = 4.0 * (qua_u - qua_s);
    dR_ds[7] = -4.0 * qua_t;
    dR_ds[8] = 0.0;
    dR_ds[9] = 4.0 * qua_t;
    
    dR_dt[0] = 1.0 - 4.0 * qua_u;
    dR_dt[1] = 0.0;
    dR_dt[2] = 0.0;
    dR_dt[3] = 4.0 * qua_t - 1.0;
    dR_dt[4] = -4.0 * qua_r;
    dR_dt[5] = 0.0;
    dR_dt[6] = -4.0 * qua_s;
    dR_dt[7] = 4.0 * (qua_u - qua_t);
    dR_dt[8] = 4.0 * qua_r;
    dR_dt[9] = 4.0 * qua_s;
    
    double xr = 0.0, xs = 0.0, xt = 0.0;
    double yr = 0.0, ys = 0.0, yt = 0.0;
    double zr = 0.0, zs = 0.0, zt = 0.0;
    for(int ii=0; ii<10; ++ii)
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
    }
  
    Matrix_double_3by3_Array mdrdx(xr, xs, xt, yr, ys, yt, zr, zs, zt);

    detJac[qua] = mdrdx.det(); // detJac = |dx/dr|
 
    mdrdx.inverse();
    
    dx_dr[9*qua + 0] = xr; dx_dr[9*qua + 1] = xs; dx_dr[9*qua + 2] = xt;
    dx_dr[9*qua + 3] = yr; dx_dr[9*qua + 4] = ys; dx_dr[9*qua + 5] = yt;
    dx_dr[9*qua + 6] = zr; dx_dr[9*qua + 7] = zs; dx_dr[9*qua + 8] = zt;
    
    dr_dx[9*qua + 0] = mdrdx(0); // dr_dx 
    dr_dx[9*qua + 1] = mdrdx(1); // dr_dy
    dr_dx[9*qua + 2] = mdrdx(2); // dr_dz
    dr_dx[9*qua + 3] = mdrdx(3); // ds_dx
    dr_dx[9*qua + 4] = mdrdx(4); // ds_dy
    dr_dx[9*qua + 5] = mdrdx(5); // ds_dz
    dr_dx[9*qua + 6] = mdrdx(6); // dt_dx
    dr_dx[9*qua + 7] = mdrdx(7); // dt_dy
    dr_dx[9*qua + 8] = mdrdx(8); // dt_dz

    for(int ii=0; ii<10; ++ii)
    {
      dR_dx[q10+ii] = dR_dr[ii]*mdrdx(0) + dR_ds[ii]*mdrdx(3) + dR_dt[ii]*mdrdx(6); 
      dR_dy[q10+ii] = dR_dr[ii]*mdrdx(1) + dR_ds[ii]*mdrdx(4) + dR_dt[ii]*mdrdx(7); 
      dR_dz[q10+ii] = dR_dr[ii]*mdrdx(2) + dR_ds[ii]*mdrdx(5) + dR_dt[ii]*mdrdx(8); 
    }

    // Setup the 6x6 matrix
    Matrix_double_6by6_Array LHS(xr, xs, xt, yr, ys, yt, zr, zs, zt);

    // LU factorization
    LHS.LU_fac();

    for(int ii=0; ii<10; ++ii)
    {
      RHS[0] = d2R_drr[ii] - dR_dx[q10+ii] * xrr - dR_dy[q10+ii] * yrr - dR_dz[q10+ii] * zrr;
      RHS[1] = d2R_drs[ii] - dR_dx[q10+ii] * xrs - dR_dy[q10+ii] * yrs - dR_dz[q10+ii] * zrs;
      RHS[2] = d2R_drt[ii] - dR_dx[q10+ii] * xrt - dR_dy[q10+ii] * yrt - dR_dz[q10+ii] * zrt;
      RHS[3] = d2R_dss[ii] - dR_dx[q10+ii] * xss - dR_dy[q10+ii] * yss - dR_dz[q10+ii] * zss;
      RHS[4] = d2R_dst[ii] - dR_dx[q10+ii] * xst - dR_dy[q10+ii] * yst - dR_dz[q10+ii] * zst;
      RHS[5] = d2R_dtt[ii] - dR_dx[q10+ii] * xtt - dR_dy[q10+ii] * ytt - dR_dz[q10+ii] * ztt;

      LHS.LU_solve(RHS, sol);

      d2R_dxx[q10+ii] = sol[0];
      d2R_dyy[q10+ii] = sol[1];
      d2R_dzz[q10+ii] = sol[2];
      d2R_dxy[q10+ii] = sol[3];
      d2R_dxz[q10+ii] = sol[4];
      d2R_dyz[q10+ii] = sol[5];
    }
  }
}

double FEAElement_Tet10_v2::get_h( const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z ) const
{
  double x,y,z,r;

  MATH_T::get_tet_sphere_info(
      ctrl_x[0], ctrl_x[1], ctrl_x[2], ctrl_x[3],
      ctrl_y[0], ctrl_y[1], ctrl_y[2], ctrl_y[3],
      ctrl_z[0], ctrl_z[1], ctrl_z[2], ctrl_z[3],
      x, y, z, r );

  return 2.0 * r;
}

void FEAElement_Tet10_v2::get_R( const int &quaindex, double * const &basis ) const
{
  assert( quaindex >= 0 && quaindex < numQuapts );
  const int offset = quaindex * 10;
  for(int ii=0; ii<10; ++ii) basis[ii] = R[offset+ii];
}

void FEAElement_Tet10_v2::get_gradR( const int &quaindex, double * const &basis_x,
    double * const &basis_y, double * const &basis_z ) const
{
  assert( quaindex >= 0 && quaindex < numQuapts );
  const int offset = quaindex * 10;
  for( int ii=0; ii<10; ++ii )
  {
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_z[ii] = dR_dz[offset + ii];
  }
}

void FEAElement_Tet10_v2::get_R_gradR( const int &quaindex, double * const &basis,
    double * const &basis_x, double * const &basis_y,
    double * const &basis_z ) const
{
  assert( quaindex >= 0 && quaindex < numQuapts );
  const int offset = quaindex * 10;
  for( int ii=0; ii<10; ++ii )
  {
    basis[ii]   = R[offset+ ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_z[ii] = dR_dz[offset + ii];
  }
}

void FEAElement_Tet10_v2::get_3D_R_dR_d2R( const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y, double * const &basis_z,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_zz, double * const &basis_xy,
    double * const &basis_xz, double * const &basis_yz ) const
{
  assert( quaindex >= 0 && quaindex < numQuapts );
  const int offset = quaindex * 10;
  for( int ii=0; ii<10; ++ii )
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

void FEAElement_Tet10_v2::get_3D_R_gradR_LaplacianR( const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y, double * const &basis_z,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_zz ) const
{
  assert( quaindex >= 0 && quaindex < numQuapts );
  const int offset = quaindex * 10;
  for( int ii=0; ii<10; ++ii )
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

void FEAElement_Tet10_v2::get_Jacobian(const int &quaindex,
    double * const &jac_value) const
{
  for(int ii=0; ii<9; ++ii) jac_value[ii] = dx_dr[9*quaindex + ii];
}

void FEAElement_Tet10_v2::get_invJacobian(const int &quaindex,
    double * const &jac_value) const
{
  for(int ii=0; ii<9; ++ii) jac_value[ii] = dr_dx[9*quaindex + ii];
}

// EOF
