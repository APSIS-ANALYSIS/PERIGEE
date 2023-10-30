#include "FEAElement_Hex27.hpp"

FEAElement_Hex27::FEAElement_Hex27( const int &in_nqua ) : numQuapts( in_nqua )
{
  R = new double [27 * numQuapts];

  dR_dx = new double [27 * numQuapts];
  dR_dy = new double [27 * numQuapts];
  dR_dz = new double [27 * numQuapts];

  d2R_dxx = new double [27 * numQuapts];
  d2R_dyy = new double [27 * numQuapts];
  d2R_dzz = new double [27 * numQuapts];
  d2R_dxy = new double [27 * numQuapts];
  d2R_dxz = new double [27 * numQuapts];
  d2R_dyz = new double [27 * numQuapts];

  dx_dr = new double [9*numQuapts];
  dr_dx = new double [9*numQuapts];
  detJac = new double [numQuapts];

  quadrilateral_face = new FEAElement_Quad9_3D_der0( numQuapts );
}

FEAElement_Hex27::~FEAElement_Hex27()
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

  delete quadrilateral_face; quadrilateral_face = nullptr;
}

void FEAElement_Hex27::print_info() const
{
  SYS_T::commPrint("Hex27: ");
  SYS_T::commPrint("27-node hexagon element with up to 2nd derivatives. \n");
  SYS_T::commPrint("elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}

double FEAElement_Hex27::get_memory_usage() const
{
  const double double_size = 289.0 * numQuapts;
  const double int_size = 1.0;
  return double_size * 8.0 + int_size * 4.0;
}

void FEAElement_Hex27::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  ASSERT( quad -> get_dim() == 3, "FEAElement_Hex27::buildBasis function error.\n" );

  for(int qua=0; qua<numQuapts; ++qua)
  {
    const int q27 = qua * 27;

    const double qua_r = quad -> get_qp( qua, 0 );
    const double qua_s = quad -> get_qp( qua, 1 );
    const double qua_t = quad -> get_qp( qua, 2 );
    
    const double Nr[3] = { (2.0 * qua_r - 1.0) * (qua_r - 1.0),
        - 4.0 * qua_r * (qua_r - 1.0), qua_r * (2.0 * qua_r - 1.0) };

    const double Ns[3] = { (2.0 * qua_s - 1.0) * (qua_s - 1.0),
        - 4.0 * qua_s * (qua_s - 1.0), qua_s * (2.0 * qua_s - 1.0) };

    const double Nt[3] = { (2.0 * qua_t - 1.0) * (qua_t - 1.0),
        - 4.0 * qua_t * (qua_t - 1.0), qua_t * (2.0 * qua_t - 1.0) };

    // vertices 0 - 7
    R[q27   ] = Nr[0] * Ns[0] * Nt[0];
    R[q27+1 ] = Nr[2] * Ns[0] * Nt[0];
    R[q27+2 ] = Nr[2] * Ns[2] * Nt[0];
    R[q27+3 ] = Nr[0] * Ns[2] * Nt[0];
    R[q27+4 ] = Nr[0] * Ns[0] * Nt[2];
    R[q27+5 ] = Nr[2] * Ns[0] * Nt[2];
    R[q27+6 ] = Nr[2] * Ns[2] * Nt[2];
    R[q27+7 ] = Nr[0] * Ns[2] * Nt[2];

    // edge 8 - 19
    R[q27+8 ] = Nr[1] * Ns[0] * Nt[0];
    R[q27+9 ] = Nr[2] * Ns[1] * Nt[0];
    R[q27+10] = Nr[1] * Ns[2] * Nt[0];
    R[q27+11] = Nr[0] * Ns[1] * Nt[0];
    R[q27+12] = Nr[1] * Ns[0] * Nt[2];
    R[q27+13] = Nr[2] * Ns[1] * Nt[2];
    R[q27+14] = Nr[1] * Ns[2] * Nt[2];
    R[q27+15] = Nr[0] * Ns[1] * Nt[2];
    R[q27+16] = Nr[0] * Ns[0] * Nt[1];
    R[q27+17] = Nr[2] * Ns[0] * Nt[1];
    R[q27+18] = Nr[2] * Ns[2] * Nt[1];
    R[q27+19] = Nr[0] * Ns[2] * Nt[1];

    // surface 20 - 25
    R[q27+20] = Nr[0] * Ns[1] * Nt[1];
    R[q27+21] = Nr[2] * Ns[1] * Nt[1];
    R[q27+22] = Nr[1] * Ns[0] * Nt[1];
    R[q27+23] = Nr[1] * Ns[2] * Nt[1];
    R[q27+24] = Nr[1] * Ns[1] * Nt[0];
    R[q27+25] = Nr[1] * Ns[1] * Nt[2];

    // center 26
    R[q27+26] = Nr[1] * Ns[1] * Nt[1];

    const double dNr[3] = { 4.0 * qua_r - 3.0, 
        - 8.0 * qua_r + 4.0, 4.0 * qua_r - 1.0 };
    const double dNs[3] = { 4.0 * qua_s - 3.0, 
        - 8.0 * qua_s + 4.0, 4.0 * qua_s - 1.0 };
    const double dNt[3] = { 4.0 * qua_t - 3.0, 
        - 8.0 * qua_t + 4.0, 4.0 * qua_t - 1.0 };

    const double dR_dr[27] { 
    dNr[0] * Ns[0] * Nt[0], dNr[2] * Ns[0] * Nt[0], dNr[2] * Ns[2] * Nt[0],
    dNr[0] * Ns[2] * Nt[0], dNr[0] * Ns[0] * Nt[2], dNr[2] * Ns[0] * Nt[2],
    dNr[2] * Ns[2] * Nt[2], dNr[0] * Ns[2] * Nt[2], dNr[1] * Ns[0] * Nt[0],
    dNr[2] * Ns[1] * Nt[0], dNr[1] * Ns[2] * Nt[0], dNr[0] * Ns[1] * Nt[0],
    dNr[1] * Ns[0] * Nt[2], dNr[2] * Ns[1] * Nt[2], dNr[1] * Ns[2] * Nt[2],
    dNr[0] * Ns[1] * Nt[2], dNr[0] * Ns[0] * Nt[1], dNr[2] * Ns[0] * Nt[1],
    dNr[2] * Ns[2] * Nt[1], dNr[0] * Ns[2] * Nt[1], dNr[0] * Ns[1] * Nt[1],
    dNr[2] * Ns[1] * Nt[1], dNr[1] * Ns[0] * Nt[1], dNr[1] * Ns[2] * Nt[1],
    dNr[1] * Ns[1] * Nt[0], dNr[1] * Ns[1] * Nt[2], dNr[1] * Ns[1] * Nt[1] };
    
    const double dR_ds[27] { 
    Nr[0] * dNs[0] * Nt[0], Nr[2] * dNs[0] * Nt[0], Nr[2] * dNs[2] * Nt[0],
    Nr[0] * dNs[2] * Nt[0], Nr[0] * dNs[0] * Nt[2], Nr[2] * dNs[0] * Nt[2],
    Nr[2] * dNs[2] * Nt[2], Nr[0] * dNs[2] * Nt[2], Nr[1] * dNs[0] * Nt[0],
    Nr[2] * dNs[1] * Nt[0], Nr[1] * dNs[2] * Nt[0], Nr[0] * dNs[1] * Nt[0],
    Nr[1] * dNs[0] * Nt[2], Nr[2] * dNs[1] * Nt[2], Nr[1] * dNs[2] * Nt[2],
    Nr[0] * dNs[1] * Nt[2], Nr[0] * dNs[0] * Nt[1], Nr[2] * dNs[0] * Nt[1],
    Nr[2] * dNs[2] * Nt[1], Nr[0] * dNs[2] * Nt[1], Nr[0] * dNs[1] * Nt[1],
    Nr[2] * dNs[1] * Nt[1], Nr[1] * dNs[0] * Nt[1], Nr[1] * dNs[2] * Nt[1],
    Nr[1] * dNs[1] * Nt[0], Nr[1] * dNs[1] * Nt[2], Nr[1] * dNs[1] * Nt[1] };
    
    const double dR_dt[27] { 
    Nr[0] * Ns[0] * dNt[0], Nr[2] * Ns[0] * dNt[0], Nr[2] * Ns[2] * dNt[0],
    Nr[0] * Ns[2] * dNt[0], Nr[0] * Ns[0] * dNt[2], Nr[2] * Ns[0] * dNt[2],
    Nr[2] * Ns[2] * dNt[2], Nr[0] * Ns[2] * dNt[2], Nr[1] * Ns[0] * dNt[0],
    Nr[2] * Ns[1] * dNt[0], Nr[1] * Ns[2] * dNt[0], Nr[0] * Ns[1] * dNt[0],
    Nr[1] * Ns[0] * dNt[2], Nr[2] * Ns[1] * dNt[2], Nr[1] * Ns[2] * dNt[2],
    Nr[0] * Ns[1] * dNt[2], Nr[0] * Ns[0] * dNt[1], Nr[2] * Ns[0] * dNt[1],
    Nr[2] * Ns[2] * dNt[1], Nr[0] * Ns[2] * dNt[1], Nr[0] * Ns[1] * dNt[1],
    Nr[2] * Ns[1] * dNt[1], Nr[1] * Ns[0] * dNt[1], Nr[1] * Ns[2] * dNt[1],
    Nr[1] * Ns[1] * dNt[0], Nr[1] * Ns[1] * dNt[2], Nr[1] * Ns[1] * dNt[1] };
    
    const double d2R_drs[27] { 
    dNr[0] * dNs[0] * Nt[0], dNr[2] * dNs[0] * Nt[0], dNr[2] * dNs[2] * Nt[0],
    dNr[0] * dNs[2] * Nt[0], dNr[0] * dNs[0] * Nt[2], dNr[2] * dNs[0] * Nt[2],
    dNr[2] * dNs[2] * Nt[2], dNr[0] * dNs[2] * Nt[2], dNr[1] * dNs[0] * Nt[0],
    dNr[2] * dNs[1] * Nt[0], dNr[1] * dNs[2] * Nt[0], dNr[0] * dNs[1] * Nt[0],
    dNr[1] * dNs[0] * Nt[2], dNr[2] * dNs[1] * Nt[2], dNr[1] * dNs[2] * Nt[2],
    dNr[0] * dNs[1] * Nt[2], dNr[0] * dNs[0] * Nt[1], dNr[2] * dNs[0] * Nt[1],
    dNr[2] * dNs[2] * Nt[1], dNr[0] * dNs[2] * Nt[1], dNr[0] * dNs[1] * Nt[1],
    dNr[2] * dNs[1] * Nt[1], dNr[1] * dNs[0] * Nt[1], dNr[1] * dNs[2] * Nt[1],
    dNr[1] * dNs[1] * Nt[0], dNr[1] * dNs[1] * Nt[2], dNr[1] * dNs[1] * Nt[1] };
    
    const double d2R_drt[27] { 
    dNr[0] * Ns[0] * dNt[0], dNr[2] * Ns[0] * dNt[0], dNr[2] * Ns[2] * dNt[0],
    dNr[0] * Ns[2] * dNt[0], dNr[0] * Ns[0] * dNt[2], dNr[2] * Ns[0] * dNt[2],
    dNr[2] * Ns[2] * dNt[2], dNr[0] * Ns[2] * dNt[2], dNr[1] * Ns[0] * dNt[0],
    dNr[2] * Ns[1] * dNt[0], dNr[1] * Ns[2] * dNt[0], dNr[0] * Ns[1] * dNt[0],
    dNr[1] * Ns[0] * dNt[2], dNr[2] * Ns[1] * dNt[2], dNr[1] * Ns[2] * dNt[2],
    dNr[0] * Ns[1] * dNt[2], dNr[0] * Ns[0] * dNt[1], dNr[2] * Ns[0] * dNt[1],
    dNr[2] * Ns[2] * dNt[1], dNr[0] * Ns[2] * dNt[1], dNr[0] * Ns[1] * dNt[1],
    dNr[2] * Ns[1] * dNt[1], dNr[1] * Ns[0] * dNt[1], dNr[1] * Ns[2] * dNt[1],
    dNr[1] * Ns[1] * dNt[0], dNr[1] * Ns[1] * dNt[2], dNr[1] * Ns[1] * dNt[1] };
    
    const double d2R_dst[27] {
    Nr[0] * dNs[0] * dNt[0], Nr[2] * dNs[0] * dNt[0], Nr[2] * dNs[2] * dNt[0],
    Nr[0] * dNs[2] * dNt[0], Nr[0] * dNs[0] * dNt[2], Nr[2] * dNs[0] * dNt[2],
    Nr[2] * dNs[2] * dNt[2], Nr[0] * dNs[2] * dNt[2], Nr[1] * dNs[0] * dNt[0],
    Nr[2] * dNs[1] * dNt[0], Nr[1] * dNs[2] * dNt[0], Nr[0] * dNs[1] * dNt[0],
    Nr[1] * dNs[0] * dNt[2], Nr[2] * dNs[1] * dNt[2], Nr[1] * dNs[2] * dNt[2],
    Nr[0] * dNs[1] * dNt[2], Nr[0] * dNs[0] * dNt[1], Nr[2] * dNs[0] * dNt[1],
    Nr[2] * dNs[2] * dNt[1], Nr[0] * dNs[2] * dNt[1], Nr[0] * dNs[1] * dNt[1],
    Nr[2] * dNs[1] * dNt[1], Nr[1] * dNs[0] * dNt[1], Nr[1] * dNs[2] * dNt[1],
    Nr[1] * dNs[1] * dNt[0], Nr[1] * dNs[1] * dNt[2], Nr[1] * dNs[1] * dNt[1] };

    const double d2Nr[3] { 4.0, -8.0, 4.0 };
    const double d2Ns[3] { 4.0, -8.0, 4.0 };
    const double d2Nt[3] { 4.0, -8.0, 4.0 };

    const double d2R_drr[27] {
    d2Nr[0] * Ns[0] * Nt[0], d2Nr[2] * Ns[0] * Nt[0], d2Nr[2] * Ns[2] * Nt[0],
    d2Nr[0] * Ns[2] * Nt[0], d2Nr[0] * Ns[0] * Nt[2], d2Nr[2] * Ns[0] * Nt[2],
    d2Nr[2] * Ns[2] * Nt[2], d2Nr[0] * Ns[2] * Nt[2], d2Nr[1] * Ns[0] * Nt[0],
    d2Nr[2] * Ns[1] * Nt[0], d2Nr[1] * Ns[2] * Nt[0], d2Nr[0] * Ns[1] * Nt[0],
    d2Nr[1] * Ns[0] * Nt[2], d2Nr[2] * Ns[1] * Nt[2], d2Nr[1] * Ns[2] * Nt[2],
    d2Nr[0] * Ns[1] * Nt[2], d2Nr[0] * Ns[0] * Nt[1], d2Nr[2] * Ns[0] * Nt[1],
    d2Nr[2] * Ns[2] * Nt[1], d2Nr[0] * Ns[2] * Nt[1], d2Nr[0] * Ns[1] * Nt[1],
    d2Nr[2] * Ns[1] * Nt[1], d2Nr[1] * Ns[0] * Nt[1], d2Nr[1] * Ns[2] * Nt[1],
    d2Nr[1] * Ns[1] * Nt[0], d2Nr[1] * Ns[1] * Nt[2], d2Nr[1] * Ns[1] * Nt[1] };
    
    const double d2R_dss[27] { 
    Nr[0] * d2Ns[0] * Nt[0], Nr[2] * d2Ns[0] * Nt[0], Nr[2] * d2Ns[2] * Nt[0],
    Nr[0] * d2Ns[2] * Nt[0], Nr[0] * d2Ns[0] * Nt[2], Nr[2] * d2Ns[0] * Nt[2],
    Nr[2] * d2Ns[2] * Nt[2], Nr[0] * d2Ns[2] * Nt[2], Nr[1] * d2Ns[0] * Nt[0],
    Nr[2] * d2Ns[1] * Nt[0], Nr[1] * d2Ns[2] * Nt[0], Nr[0] * d2Ns[1] * Nt[0],
    Nr[1] * d2Ns[0] * Nt[2], Nr[2] * d2Ns[1] * Nt[2], Nr[1] * d2Ns[2] * Nt[2],
    Nr[0] * d2Ns[1] * Nt[2], Nr[0] * d2Ns[0] * Nt[1], Nr[2] * d2Ns[0] * Nt[1],
    Nr[2] * d2Ns[2] * Nt[1], Nr[0] * d2Ns[2] * Nt[1], Nr[0] * d2Ns[1] * Nt[1],
    Nr[2] * d2Ns[1] * Nt[1], Nr[1] * d2Ns[0] * Nt[1], Nr[1] * d2Ns[2] * Nt[1],
    Nr[1] * d2Ns[1] * Nt[0], Nr[1] * d2Ns[1] * Nt[2], Nr[1] * d2Ns[1] * Nt[1] };
    
    const double d2R_dtt[27] { 
    Nr[0] * Ns[0] * d2Nt[0], Nr[2] * Ns[0] * d2Nt[0], Nr[2] * Ns[2] * d2Nt[0],
    Nr[0] * Ns[2] * d2Nt[0], Nr[0] * Ns[0] * d2Nt[2], Nr[2] * Ns[0] * d2Nt[2],
    Nr[2] * Ns[2] * d2Nt[2], Nr[0] * Ns[2] * d2Nt[2], Nr[1] * Ns[0] * d2Nt[0],
    Nr[2] * Ns[1] * d2Nt[0], Nr[1] * Ns[2] * d2Nt[0], Nr[0] * Ns[1] * d2Nt[0],
    Nr[1] * Ns[0] * d2Nt[2], Nr[2] * Ns[1] * d2Nt[2], Nr[1] * Ns[2] * d2Nt[2],
    Nr[0] * Ns[1] * d2Nt[2], Nr[0] * Ns[0] * d2Nt[1], Nr[2] * Ns[0] * d2Nt[1],
    Nr[2] * Ns[2] * d2Nt[1], Nr[0] * Ns[2] * d2Nt[1], Nr[0] * Ns[1] * d2Nt[1],
    Nr[2] * Ns[1] * d2Nt[1], Nr[1] * Ns[0] * d2Nt[1], Nr[1] * Ns[2] * d2Nt[1],
    Nr[1] * Ns[1] * d2Nt[0], Nr[1] * Ns[1] * d2Nt[2], Nr[1] * Ns[1] * d2Nt[1] };
    
    double xr = 0.0, xs = 0.0, xt = 0.0;
    double yr = 0.0, ys = 0.0, yt = 0.0;
    double zr = 0.0, zs = 0.0, zt = 0.0;
    
    double xrr = 0.0, xss = 0.0, xtt = 0.0;
    double yrr = 0.0, yss = 0.0, ytt = 0.0;
    double zrr = 0.0, zss = 0.0, ztt = 0.0;

    double xrs = 0.0, xrt = 0.0, xst = 0.0;
    double yrs = 0.0, yrt = 0.0, yst = 0.0;
    double zrs = 0.0, zrt = 0.0, zst = 0.0;

    for(int ii=0; ii<27; ++ii)
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
      
      xrr += ctrl_x[ii] * d2R_drr[ii];
      xss += ctrl_x[ii] * d2R_dss[ii];
      xtt += ctrl_x[ii] * d2R_dtt[ii];
      
      yrr += ctrl_y[ii] * d2R_drr[ii];
      yss += ctrl_y[ii] * d2R_dss[ii];
      ytt += ctrl_y[ii] * d2R_dtt[ii];
      
      zrr += ctrl_z[ii] * d2R_drr[ii];
      zss += ctrl_z[ii] * d2R_dss[ii];
      ztt += ctrl_z[ii] * d2R_dtt[ii];

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

    for(int ii=0; ii<27; ++ii)
    {
      dR_dx[q27+ii] = dR_dr[ii]*mdrdx(0) + dR_ds[ii]*mdrdx(3) + dR_dt[ii]*mdrdx(6); 
      dR_dy[q27+ii] = dR_dr[ii]*mdrdx(1) + dR_ds[ii]*mdrdx(4) + dR_dt[ii]*mdrdx(7); 
      dR_dz[q27+ii] = dR_dr[ii]*mdrdx(2) + dR_ds[ii]*mdrdx(5) + dR_dt[ii]*mdrdx(8); 
    }

    // Setup the 6x6 matrix
    FE_T::Matrix_double_6by6_Array LHS(xr, xs, xt, yr, ys, yt, zr, zs, zt);

    // LU factorization
    LHS.LU_fac();

    for(int ii=0; ii<27; ++ii)
    {
      const std::array<double, 6> RHS {{ d2R_drr[ii] - dR_dx[q27+ii] * xrr - dR_dy[q27+ii] * yrr - dR_dz[q27+ii] * zrr,
        d2R_drs[ii] - dR_dx[q27+ii] * xrs - dR_dy[q27+ii] * yrs - dR_dz[q27+ii] * zrs,
        d2R_drt[ii] - dR_dx[q27+ii] * xrt - dR_dy[q27+ii] * yrt - dR_dz[q27+ii] * zrt,
        d2R_dss[ii] - dR_dx[q27+ii] * xss - dR_dy[q27+ii] * yss - dR_dz[q27+ii] * zss,
        d2R_dst[ii] - dR_dx[q27+ii] * xst - dR_dy[q27+ii] * yst - dR_dz[q27+ii] * zst,
        d2R_dtt[ii] - dR_dx[q27+ii] * xtt - dR_dy[q27+ii] * ytt - dR_dz[q27+ii] * ztt }};

      const auto sol = LHS.LU_solve(RHS);

      d2R_dxx[q27+ii] = sol[0];
      d2R_dyy[q27+ii] = sol[1];
      d2R_dzz[q27+ii] = sol[2];
      d2R_dxy[q27+ii] = sol[3];
      d2R_dxz[q27+ii] = sol[4];
      d2R_dyz[q27+ii] = sol[5];
    }
  }
}

double FEAElement_Hex27::get_h( const double * const &ctrl_x,
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

void FEAElement_Hex27::get_R( const int &quaindex, double * const &basis ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_R function error.\n" );
  const int offset = quaindex * 27;
  for(int ii=0; ii<27; ++ii) basis[ii] = R[offset+ii];
}

std::vector<double> FEAElement_Hex27::get_R( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_R function error.\n" );
  const int offset = quaindex * 27;
  std::vector<double> vec(R + offset, R + offset + 27);
  return vec;
}

void FEAElement_Hex27::get_gradR( const int &quaindex, double * const &basis_x,
    double * const &basis_y, double * const &basis_z ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_gradR function error.\n" );
  const int offset = quaindex * 27;
  for( int ii=0; ii<27; ++ii )
  {
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_z[ii] = dR_dz[offset + ii];
  }
}

void FEAElement_Hex27::get_R_gradR( const int &quaindex, double * const &basis,
    double * const &basis_x, double * const &basis_y,
    double * const &basis_z ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_R_gradR function error.\n" );
  const int offset = quaindex * 27;
  for( int ii=0; ii<27; ++ii )
  {
    basis[ii]   = R[offset+ ii];
    basis_x[ii] = dR_dx[offset + ii];
    basis_y[ii] = dR_dy[offset + ii];
    basis_z[ii] = dR_dz[offset + ii];
  }
}

std::vector<double> FEAElement_Hex27::get_dR_dx( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_dR_dx function error.\n" );
  const int offset = quaindex * 27;
  std::vector<double> vec(dR_dx + offset, dR_dx + offset + 27);
  return vec;
}

std::vector<double> FEAElement_Hex27::get_dR_dy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_dR_dy function error.\n" );
  const int offset = quaindex * 27;
  std::vector<double> vec(dR_dy + offset, dR_dy + offset + 27);
  return vec;
}

std::vector<double> FEAElement_Hex27::get_dR_dz( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_dR_dz function error.\n" );
  const int offset = quaindex * 27;
  std::vector<double> vec(dR_dz + offset, dR_dz + offset + 27);
  return vec;
}

void FEAElement_Hex27::get_3D_R_dR_d2R( const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y, double * const &basis_z,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_zz, double * const &basis_xy,
    double * const &basis_xz, double * const &basis_yz ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_3D_R_dR_d2R function error.\n" );
  const int offset = quaindex * 27;
  for( int ii=0; ii<27; ++ii )
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

void FEAElement_Hex27::get_3D_R_gradR_LaplacianR( const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y, double * const &basis_z,
    double * const &basis_xx, double * const &basis_yy,
    double * const &basis_zz ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_3D_R_gradR_LaplacianR function error.\n" );
  const int offset = quaindex * 27;
  for( int ii=0; ii<27; ++ii )
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

std::vector<double> FEAElement_Hex27::get_d2R_dxx( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_d2R_dxx function error.\n" );
  const int offset = quaindex * 27;
  std::vector<double> vec(d2R_dxx + offset, d2R_dxx + offset + 27);
  return vec;
}

std::vector<double> FEAElement_Hex27::get_d2R_dyy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_d2R_dyy function error.\n" );
  const int offset = quaindex * 27;
  std::vector<double> vec(d2R_dyy + offset, d2R_dyy + offset + 27);
  return vec;
}

std::vector<double> FEAElement_Hex27::get_d2R_dzz( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_d2R_dzz function error.\n" );
  const int offset = quaindex * 27;
  std::vector<double> vec(d2R_dzz + offset, d2R_dzz + offset + 27);
  return vec;
}

std::vector<double> FEAElement_Hex27::get_d2R_dxy( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_d2R_dxy function error.\n" );
  const int offset = quaindex * 27;
  std::vector<double> vec(d2R_dxy + offset, d2R_dxy + offset + 27);
  return vec;
}

std::vector<double> FEAElement_Hex27::get_d2R_dxz( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_d2R_dxz function error.\n" );
  const int offset = quaindex * 27;
  std::vector<double> vec(d2R_dxz + offset, d2R_dxz + offset + 27);
  return vec;
}

std::vector<double> FEAElement_Hex27::get_d2R_dyz( const int &quaindex ) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_d2R_dyz function error.\n" );
  const int offset = quaindex * 27;
  std::vector<double> vec(d2R_dyz + offset, d2R_dyz + offset + 27);
  return vec;
}

void FEAElement_Hex27::get_Jacobian(const int &quaindex,
    double * const &jac_value) const
{
  for(int ii=0; ii<9; ++ii) jac_value[ii] = dx_dr[9*quaindex + ii];
}

std::array<double,9> FEAElement_Hex27::get_Jacobian(const int &quaindex) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_Jacobian function error.\n" );
  return {{ dx_dr[9*quaindex], dx_dr[9*quaindex+1], dx_dr[9*quaindex+2],
    dx_dr[9*quaindex+3], dx_dr[9*quaindex+4], dx_dr[9*quaindex+5],
    dx_dr[9*quaindex+6], dx_dr[9*quaindex+7], dx_dr[9*quaindex+8] }};
}

void FEAElement_Hex27::get_invJacobian(const int &quaindex,
    double * const &jac_value) const
{
  for(int ii=0; ii<9; ++ii) jac_value[ii] = dr_dx[9*quaindex + ii];
}

std::array<double,9> FEAElement_Hex27::get_invJacobian(const int &quaindex) const
{
  ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Hex27::get_invJacobian function error.\n" );
  return {{ dr_dx[9*quaindex], dr_dx[9*quaindex+1], dr_dx[9*quaindex+2],
    dr_dx[9*quaindex+3], dr_dx[9*quaindex+4], dr_dx[9*quaindex+5],
    dr_dx[9*quaindex+6], dr_dx[9*quaindex+7], dr_dx[9*quaindex+8] }};
}

void FEAElement_Hex27::buildBasis( const IQuadPts * const &quad_s, const int &face_id,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  // Build the volume element
  const auto quad_v = FE_T::QuadPts_Gauss_on_boundary( this->get_Type(), face_id, quad_s );
  this->buildBasis( &quad_v, ctrl_x, ctrl_y, ctrl_z );

  std::vector<double> face_ctrl_x( 9, 0.0 ), face_ctrl_y( 9, 0.0 ), face_ctrl_z( 9, 0.0 );

  switch( face_id )
  {
    case 0:
      face_ctrl_x = std::vector<double> { ctrl_x[0], ctrl_x[3], ctrl_x[2], ctrl_x[1],
                                          ctrl_x[11], ctrl_x[10],ctrl_x[9], ctrl_x[8], ctrl_x[24] };
      face_ctrl_y = std::vector<double> { ctrl_y[0], ctrl_y[3], ctrl_y[2], ctrl_y[1],
                                          ctrl_y[11], ctrl_y[10],ctrl_y[9], ctrl_y[8], ctrl_y[24] };
      face_ctrl_z = std::vector<double> { ctrl_z[0], ctrl_z[3], ctrl_z[2], ctrl_z[1],
                                          ctrl_z[11], ctrl_z[10],ctrl_z[9], ctrl_z[8], ctrl_z[24] };
      break;

    case 1:
      face_ctrl_x = std::vector<double> { ctrl_x[4], ctrl_x[5], ctrl_x[6], ctrl_x[7],
                                          ctrl_x[12], ctrl_x[13],ctrl_x[14], ctrl_x[15], ctrl_x[25] };
      face_ctrl_y = std::vector<double> { ctrl_y[4], ctrl_y[5], ctrl_y[6], ctrl_y[7],
                                          ctrl_y[12], ctrl_y[13],ctrl_y[14], ctrl_y[15], ctrl_y[25] };
      face_ctrl_z = std::vector<double> { ctrl_z[4], ctrl_z[5], ctrl_z[6], ctrl_z[7],
                                          ctrl_z[12], ctrl_z[13],ctrl_z[14], ctrl_z[15], ctrl_z[25] };
      break;

    case 2:
      face_ctrl_x = std::vector<double> { ctrl_x[0], ctrl_x[1], ctrl_x[5], ctrl_x[4],
                                          ctrl_x[8], ctrl_x[17],ctrl_x[12], ctrl_x[16], ctrl_x[22] };
      face_ctrl_y = std::vector<double> { ctrl_y[0], ctrl_y[1], ctrl_y[5], ctrl_y[4],
                                          ctrl_y[8], ctrl_y[17],ctrl_y[12], ctrl_y[16], ctrl_y[22] };
      face_ctrl_z = std::vector<double> { ctrl_z[0], ctrl_z[1], ctrl_z[5], ctrl_z[4],
                                          ctrl_z[8], ctrl_z[17],ctrl_z[12], ctrl_z[16], ctrl_z[22] };
      break;

    case 3:
      face_ctrl_x = std::vector<double> { ctrl_x[1], ctrl_x[2], ctrl_x[6], ctrl_x[5],
                                          ctrl_x[9], ctrl_x[18],ctrl_x[13], ctrl_x[17], ctrl_x[21] };
      face_ctrl_y = std::vector<double> { ctrl_y[1], ctrl_y[2], ctrl_y[6], ctrl_y[5],
                                          ctrl_y[9], ctrl_y[18],ctrl_y[13], ctrl_y[17], ctrl_y[21] };
      face_ctrl_z = std::vector<double> { ctrl_z[1], ctrl_z[2], ctrl_z[6], ctrl_z[5],
                                          ctrl_z[9], ctrl_z[18],ctrl_z[13], ctrl_z[17], ctrl_z[21] };
      break;

    case 4:
      face_ctrl_x = std::vector<double> { ctrl_x[3], ctrl_x[7], ctrl_x[6], ctrl_x[2],
                                          ctrl_x[19], ctrl_x[14],ctrl_x[18], ctrl_x[10], ctrl_x[23] };
      face_ctrl_y = std::vector<double> { ctrl_y[3], ctrl_y[7], ctrl_y[6], ctrl_y[2],
                                          ctrl_y[19], ctrl_y[14],ctrl_y[18], ctrl_y[10], ctrl_y[23] };
      face_ctrl_z = std::vector<double> { ctrl_z[3], ctrl_z[7], ctrl_z[6], ctrl_z[2],
                                          ctrl_z[19], ctrl_z[14],ctrl_z[18], ctrl_z[10], ctrl_z[23] };
      break;

    case 5:
      face_ctrl_x = std::vector<double> { ctrl_x[0], ctrl_x[4], ctrl_x[7], ctrl_x[3],
                                          ctrl_x[16], ctrl_x[15],ctrl_x[19], ctrl_x[11], ctrl_x[20] };
      face_ctrl_y = std::vector<double> { ctrl_y[0], ctrl_y[4], ctrl_y[7], ctrl_y[3],
                                          ctrl_y[16], ctrl_y[15],ctrl_y[19], ctrl_y[11], ctrl_y[20] };
      face_ctrl_z = std::vector<double> { ctrl_z[0], ctrl_z[4], ctrl_z[7], ctrl_z[3],
                                          ctrl_z[16], ctrl_z[15],ctrl_z[19], ctrl_z[11], ctrl_z[20] };
      break;

    default:
      SYS_T::print_fatal("Error: FEAElement_Hex27::buildBoundaryBasis, wrong face id.\n");
      break;
  }

  quadrilateral_face->buildBasis( quad_s, &face_ctrl_x[0], &face_ctrl_y[0], &face_ctrl_z[0] );
}

// EOF
