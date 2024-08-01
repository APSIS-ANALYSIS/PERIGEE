#include "FEAElement_Tet_P2_P1.hpp"

FEAElement_Tet_P2_P1::FEAElement_Tet_P2_P1( const int &in_nqua )
: numQuapts( in_nqua )
{
  R = new double [14 * numQuapts];
  
  dR_dx = new double [14 * numQuapts];
  dR_dy = new double [14 * numQuapts];
  dR_dz = new double [14 * numQuapts];

  dx_dr = new double [9*numQuapts];
  dr_dx = new double [9*numQuapts];
  detJac = new double [numQuapts];
}


FEAElement_Tet_P2_P1::~FEAElement_Tet_P2_P1()
{
  delete [] R; R = NULL;
  delete [] dR_dx; dR_dx = NULL;
  delete [] dR_dy; dR_dy = NULL;
  delete [] dR_dz; dR_dz = NULL;
  delete [] dx_dr; dx_dr = NULL;
  delete [] dr_dx; dr_dx = NULL;
  delete [] detJac; detJac = NULL;
}


void FEAElement_Tet_P2_P1::print_info() const
{
  SYS_T::commPrint("Tet_P2_P1: First 10 nodes for disp field, additional 4 nodes for pressure field. \n");
  SYS_T::commPrint("10-node tetrahedral element with up to 1st derivatives. \n");
  SYS_T::commPrint("4-node tetrahedral element with up to 1st derivatives. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}


double FEAElement_Tet_P2_P1::get_memory_usage() const
{
  const double dsize = (14*4+19)*numQuapts;
  return dsize * 8.0 + 4.0;
}


double FEAElement_Tet_P2_P1::get_h( const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z ) const
{
  return 2.0 * FE_T::get_tet_sphere_radius(
      ctrl_x[0], ctrl_x[1], ctrl_x[2], ctrl_x[3],
      ctrl_y[0], ctrl_y[1], ctrl_y[2], ctrl_y[3],
      ctrl_z[0], ctrl_z[1], ctrl_z[2], ctrl_z[3]
      );  
}


void FEAElement_Tet_P2_P1::get_R( const int &quaindex, 
    double * const &basis ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet_P2_P1::get_R function error.\n" );
  const int offset = quaindex * 14;
  for(int ii=0; ii<14; ++ii) basis[ii] = R[offset+ii];
}


void FEAElement_Tet_P2_P1::get_gradR( const int &quaindex, 
    double * const &basis_x, double * const &basis_y, 
    double * const &basis_z ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet_P2_P1::get_gradR function error.\n" );
  const int offset = quaindex * 14;
  for( int ii=0; ii<14; ++ii )
  {
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_z[ii] = dR_dz[offset + ii];
  }
}


void FEAElement_Tet_P2_P1::get_R_gradR( const int &quaindex, 
    double * const &basis, double * const &basis_x, 
    double * const &basis_y, double * const &basis_z ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Tet_P2_P1::get_R_gradR function error.\n" );
  const int offset = quaindex * 14;
  for( int ii=0; ii<14; ++ii )
  {
    basis[ii]   = R[offset+ ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_z[ii] = dR_dz[offset + ii];
  }
}


void FEAElement_Tet_P2_P1::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  ASSERT( quad -> get_dim() == 4, "FEAElement_Tet_P2_P1::buildBasis function error.\n" );

  double qua_r, qua_s, qua_t, qua_u;
  for(int qua=0; qua<numQuapts; ++qua)
  {
    const int q14 = qua * 14;

    qua_r = quad -> get_qp( qua, 0 );
    qua_s = quad -> get_qp( qua, 1 );
    qua_t = quad -> get_qp( qua, 2 );
    qua_u = quad -> get_qp( qua, 3 );

    // P2 part
    R[q14+0] = qua_u * (2.0*qua_u - 1.0);
    R[q14+1] = qua_r * (2.0*qua_r - 1.0);
    R[q14+2] = qua_s * (2.0*qua_s - 1.0);
    R[q14+3] = qua_t * (2.0*qua_t - 1.0);
    R[q14+4] = 4.0 * qua_u * qua_r;
    R[q14+5] = 4.0 * qua_r * qua_s;
    R[q14+6] = 4.0 * qua_s * qua_u;
    R[q14+7] = 4.0 * qua_t * qua_u;
    R[q14+8] = 4.0 * qua_s * qua_t;
    R[q14+9] = 4.0 * qua_r * qua_t;

    // P1 part
    R[q14+10] = qua_u;
    R[q14+11] = qua_r;
    R[q14+12] = qua_s;
    R[q14+13] = qua_t;

    // 1st derivative
    dR_dr[0] = 1.0 - 4.0 * qua_u;
    dR_dr[1] = 4.0 * qua_r - 1.0;
    dR_dr[2] = 0.0;
    dR_dr[3] = 0.0;
    dR_dr[4] = 4.0 * (qua_u - qua_r);
    dR_dr[5] = 4.0 * qua_s;
    dR_dr[6] = -4.0 * qua_s;
    dR_dr[7] = -4.0 * qua_t;
    dR_dr[8] = 0.0;
    dR_dr[9] = 4.0 * qua_t;

    dR_dr[10] = -1.0;
    dR_dr[11] = 1.0;
    dR_dr[12] = 0.0;
    dR_dr[13] = 0.0;

    dR_ds[0] = 1.0 - 4.0 * qua_u;
    dR_ds[1] = 0.0;
    dR_ds[2] = 4.0 * qua_s - 1.0;
    dR_ds[3] = 0.0;
    dR_ds[4] = -4.0 * qua_r;
    dR_ds[5] = 4.0 * qua_r;
    dR_ds[6] = 4.0 * (qua_u - qua_s);
    dR_ds[7] = -4.0 * qua_t;
    dR_ds[8] = 4.0 * qua_t;
    dR_ds[9] = 0.0;

    dR_ds[10] = -1.0;
    dR_ds[11] = 0.0;
    dR_ds[12] = 1.0;
    dR_ds[13] = 0.0;

    dR_dt[0] = 1.0 - 4.0 * qua_u;
    dR_dt[1] = 0.0;
    dR_dt[2] = 0.0;
    dR_dt[3] = 4.0 * qua_t - 1.0;
    dR_dt[4] = -4.0 * qua_r;
    dR_dt[5] = 0.0;
    dR_dt[6] = -4.0 * qua_s;
    dR_dt[7] = 4.0 * (qua_u - qua_t);
    dR_dt[8] = 4.0 * qua_s;
    dR_dt[9] = 4.0 * qua_r;

    dR_dt[10] = -1.0;
    dR_dt[11] = 0.0;
    dR_dt[12] = 0.0;
    dR_dt[13] = 1.0;

    double xr = 0.0, xs = 0.0, xt = 0.0;
    double yr = 0.0, ys = 0.0, yt = 0.0;
    double zr = 0.0, zs = 0.0, zt = 0.0;
    // First ten node define geometry
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

    FE_T::Matrix_double_3by3_Array mdrdx(xr, xs, xt, yr, ys, yt, zr, zs, zt);

    detJac[qua] = mdrdx.det();

    mdrdx.inverse();

    // Define dx_dr and dr_dx arrays
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

    // Now physical derivatives
    for(int ii=0; ii<14; ++ii)
    {
      dR_dx[q14+ii] = dR_dr[ii]*mdrdx(0) + dR_ds[ii]*mdrdx(3) + dR_dt[ii]*mdrdx(6);
      dR_dy[q14+ii] = dR_dr[ii]*mdrdx(1) + dR_ds[ii]*mdrdx(4) + dR_dt[ii]*mdrdx(7);
      dR_dz[q14+ii] = dR_dr[ii]*mdrdx(2) + dR_ds[ii]*mdrdx(5) + dR_dt[ii]*mdrdx(8);
    }
  } // End qua-loop
}

// EOF
