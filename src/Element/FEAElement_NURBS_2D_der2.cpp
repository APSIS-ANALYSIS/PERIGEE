#include "FEAElement_NURBS_2D_der2.hpp"

FEAElement_NURBS_2D_der2::FEAElement_NURBS_2D_der2( const int &in_eIndex,
    const IALocal_meshSize * const &mSize,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const FEANode * const &feaNode,
    const IAExtractor * const &extractor,
    const ALocal_IEN * const &locIEN )
{
  elem_index = in_eIndex;

  int sDegree = bs->get_degree();
  int tDegree = bt->get_degree();

  int sDp1 = sDegree + 1; int tDp1 = tDegree + 1;

  nLocBas = sDp1 * tDp1;
  
  int num_qua_s = bs->get_nQuapts();
  int num_qua_t = bt->get_nQuapts();

  numQuapts = num_qua_s * num_qua_t;

  is_sNonzero = true;
  double refSize = mSize->get_hx(elem_index) * mSize->get_hy(elem_index);

  if(refSize == 0.0)
    is_sNonzero = false;

  R = NULL; dR_dx = NULL; dR_dy = NULL; detJac = NULL;
  d2R_dxx = NULL; d2R_dxy = NULL; d2R_dyy = NULL;

  // Calculate the quadrature of basis functions
  // memory layout is
  // basis0 qua0, basis1, qua0, ... basis_nLocbas-1, qua0, ...
  // basis_nLocbas-1, qua_numquapts-1
  //
  // Jacobian: dx_dxi, dx_deta, dy_dxi, dy_deta at qua0, ..., at last qua
  if(is_sNonzero == true)
  {
    buildBasis( sDp1, tDp1, num_qua_s, num_qua_t, mSize, bs, bt, 
        feaNode, extractor, locIEN );
  }
}


FEAElement_NURBS_2D_der2::~FEAElement_NURBS_2D_der2()
{
  clearBasisCache();
}


void FEAElement_NURBS_2D_der2::clearBasisCache()
{
  if(is_sNonzero)
  {
    delete [] R;        delete [] dR_dx;   delete [] dR_dy;
    delete [] d2R_dxx;  delete [] d2R_dxy; delete [] d2R_dyy;
    delete [] detJac;
  }
  R = NULL; dR_dx = NULL; dR_dy = NULL; detJac = NULL;
  d2R_dxx = NULL; d2R_dxy = NULL; d2R_dyy = NULL;
}


void FEAElement_NURBS_2D_der2::buildBasis( const int &ssdp1, const int &ttdp1,
    const int &num_qs, const int &num_qt,
    const IALocal_meshSize * const &mSize,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const FEANode * const &feaNode,
    const IAExtractor * const &extractor,
    const ALocal_IEN * const &locIEN )
{
  int basis_size = nLocBas * numQuapts;

  R      = new double [basis_size];
  dR_dx  = new double [basis_size];
  dR_dy  = new double [basis_size];
  d2R_dxx = new double [basis_size];
  d2R_dxy = new double [basis_size];
  d2R_dyy = new double [basis_size];
  detJac = new double [basis_size];

  double hx = mSize->get_hx(elem_index);
  double hy = mSize->get_hy(elem_index);

  double * ext_x; double * ext_y;
  extractor->get_EXT_x(elem_index, ext_x);
  extractor->get_EXT_y(elem_index, ext_y);

  double * ctrl_x = new double [nLocBas];
  double * ctrl_y = new double [nLocBas];
  double * ctrl_z = new double [nLocBas];
  double * ctrl_w = new double [nLocBas];

  for(int ii=0; ii<nLocBas; ++ii)
  {
    int index = locIEN->get_LIEN(elem_index, ii);
    ctrl_x[ii] = feaNode->get_ctrlPts_x(index);
    ctrl_y[ii] = feaNode->get_ctrlPts_y(index);
    ctrl_z[ii] = feaNode->get_ctrlPts_z(index);
    ctrl_w[ii] = feaNode->get_ctrlPts_w(index);
  }

  // Initialize B-Spline functions
  double * Nn     = new double [nLocBas];
  double * dNn_ds = new double [nLocBas];
  double * dNn_dt = new double [nLocBas];

  double * d2Nn_dss = new double [nLocBas];
  double * d2Nn_dst = new double [nLocBas];
  double * d2Nn_dtt = new double [nLocBas];

  double * dRr_ds = new double [nLocBas];
  double * dRr_dt = new double [nLocBas];
  
  double * d2Rr_dss = new double [nLocBas];
  double * d2Rr_dst = new double [nLocBas];
  double * d2Rr_dtt = new double [nLocBas];
  
  double * Nns = new double [ssdp1];
  double * Nnt = new double [ttdp1];

  double * dNns_ds = new double [ssdp1];
  double * dNnt_dt = new double [ttdp1];
  
  double * d2Nns_dss = new double [ssdp1];
  double * d2Nnt_dtt = new double [ttdp1];

  double * Bbs       = new double [ssdp1];
  double * dBbs_ds   = new double [ssdp1];
  double * d2Bbs_dss = new double [ssdp1];
  double * Bbt       = new double [ttdp1];
  double * dBbt_dt   = new double [ttdp1];
  double * d2Bbt_dtt = new double [ttdp1];

  int ii, jj, ll;
  int counter = -1;
  const double invhx = 1.0 / hx; 
  const double invhy = 1.0 / hy;
  double extor;

  for(int qua_y=0; qua_y < num_qt; ++qua_y)
  {
    for(ii=0; ii<ttdp1; ++ii)
    {
      Bbt[ii]       = bt->get_der0(ii, qua_y);
      dBbt_dt[ii]   = bt->get_der1(ii, qua_y) * invhy;
      d2Bbt_dtt[ii] = bt->get_der2(ii, qua_y) * invhy * invhy;
    }
    ll = -1;
    for(ii=0; ii<ttdp1; ++ii)
    {
      Nnt[ii] = 0.0; dNnt_dt[ii] = 0.0; d2Nnt_dtt[ii] = 0.0;
      for(jj=0; jj<ttdp1; ++jj)
      {
        ll += 1;
        extor = ext_y[ll];
        Nnt[ii]       += extor * Bbt[jj];
        dNnt_dt[ii]   += extor * dBbt_dt[jj];
        d2Nnt_dtt[ii] += extor * d2Bbt_dtt[jj];
      }
    }
    
    for(int qua_x = 0; qua_x < num_qs; ++qua_x)
    {
      for(ii=0; ii<ssdp1; ++ii)
      {
        Bbs[ii]       = bs->get_der0(ii, qua_x);
        dBbs_ds[ii]   = bs->get_der1(ii, qua_x) * invhx;
        d2Bbs_dss[ii] = bs->get_der2(ii, qua_x) * invhx * invhx;
      }
      ll = -1;
      for(ii=0; ii<ssdp1; ++ii)
      {
        Nns[ii] = 0.0; dNns_ds[ii] = 0.0; d2Nns_dss[ii] = 0.0;
        for(jj=0; jj<ssdp1; ++jj)
        {
          ll += 1;
          extor = ext_x[ll];
          Nns[ii]       += extor * Bbs[jj];
          dNns_ds[ii]   += extor * dBbs_ds[jj];
          d2Nns_dss[ii] += extor * d2Bbs_dss[jj];
        }
      }
      counter += 1;

      BuildShape_atQua( ssdp1, ttdp1, counter, qua_x, qua_y, hx, hy,
          ctrl_x, ctrl_y, ctrl_z, ctrl_w, Nn, dNn_ds, dNn_dt,
          d2Nn_dss, d2Nn_dtt, d2Nn_dst, dRr_ds, dRr_dt, d2Rr_dss, 
          d2Rr_dtt, d2Rr_dst, Nns, Nnt, dNns_ds, dNnt_dt,
          d2Nns_dss, d2Nnt_dtt );
    } 
  }

  // control pts
  delete [] ctrl_x; delete [] ctrl_y; delete [] ctrl_z; delete [] ctrl_w;
  // extraction
  delete [] ext_x; delete [] ext_y;
  // B splines
  delete [] Nn; delete [] dNn_ds; delete [] dNn_dt;
  delete [] d2Nn_dss; delete [] d2Nn_dtt; delete [] d2Nn_dst;
  // NURBS on parametric domain
  delete [] dRr_ds; delete [] dRr_dt; delete [] d2Rr_dss; delete [] d2Rr_dtt;
  delete [] d2Rr_dst; 
  // univariate B-splines
  delete [] Nns; delete [] Nnt;
  delete [] dNns_ds; delete [] dNnt_dt;
  delete [] d2Nns_dss; delete [] d2Nnt_dtt;
  // Bernstein polynomial
  delete [] Bbs; delete [] Bbt;  delete [] dBbs_ds; delete [] dBbt_dt;
  delete [] d2Bbs_dss; delete [] d2Bbt_dtt;
}

void FEAElement_NURBS_2D_der2::get_R(
    const int &quaindex, double * const &basis ) const
{
  for(int ii=0; ii<nLocBas; ++ii)
    basis[ii] = R[quaindex * nLocBas + ii];
}


void FEAElement_NURBS_2D_der2::get_R_gradR(const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y ) const
{
  int offset = quaindex * nLocBas;
  int loc;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    loc = offset + ii;
    basis[ii]   = R[loc];
    basis_x[ii] = dR_dx[loc];
    basis_y[ii] = dR_dy[loc];
  }
}

void FEAElement_NURBS_2D_der2::get_2D_R_gradR_LaplacianR( 
    const int &quaindex, double * const &basis, double * const &basis_x, 
    double * const &basis_y, double * const &basis_xx, double * const &basis_yy ) const
{
  int offset = quaindex * nLocBas;
  int loc;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    loc = offset + ii;
    basis[ii]    = R[loc];
    basis_x[ii]  = dR_dx[loc];
    basis_y[ii]  = dR_dy[loc];
    basis_xx[ii] = d2R_dxx[loc];
    basis_yy[ii] = d2R_dyy[loc];
  }
}

void FEAElement_NURBS_2D_der2::print() const
{
  SYS_T::commPrint("NURBS_2D_der2: ");
  SYS_T::commPrint("2-Dimensional NURBS shape function with 2nd order derivatives. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d. \n", get_Type());
  SYS_T::commPrint("Note: Jacobian matrix is not evaluated. \n");
}



double FEAElement_NURBS_2D_der2::get_memory_usage() const
{
  unsigned int double_size = 0; unsigned int int_size = 3;
  double_size = nLocBas * numQuapts * 6 + numQuapts;
  return (double) double_size * 8.0 + (double) int_size * 4.0;
}


void FEAElement_NURBS_2D_der2::BuildShape_atQua(
    const int &sdp1, const int &tdp1,
    const int &quaindex,
    const int &quaindex_s, const int &quaindex_t,
    const double &hx, const double &hy,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z,
    const double * const &ctrl_w,
    double * &N, double * &dN_ds, double * &dN_dt,
    double * &d2N_dss, double * &d2N_dtt, double * &d2N_dst,
    double * &dR_ds, double * &dR_dt,
    double * &d2R_dss, double * &d2R_dtt, double * &d2R_dst,
    const double * const &Ns,
    const double * const &Nt,
    const double * const &dNs_ds,
    const double * const &dNt_dt,
    const double * const &d2Ns_dss,
    const double * const &d2Nt_dtt )
{
  int ii, jj, ll;

  double detJac_temp = 0.0;

  double w = 0.0; double dw_ds = 0.0; double dw_dt = 0.0;
  double d2w_dss = 0.0; double d2w_dtt = 0.0; double d2w_dst = 0.0;

  ll = -1;

  for(jj=0; jj<tdp1; ++jj)
  {
    for(ii=0; ii<sdp1; ++ii)
    {
      ll += 1;
      N[ll] = Ns[ii] * Nt[jj] * ctrl_w[ll];
      w += N[ll];

      dN_ds[ll] = dNs_ds[ii] * Nt[jj] * ctrl_w[ll];
      dN_dt[ll] = Ns[ii] * dNt_dt[jj] * ctrl_w[ll];

      d2N_dss[ll] = d2Ns_dss[ii] * Nt[jj] * ctrl_w[ll];
      d2N_dtt[ll] = Ns[ii] * d2Nt_dtt[jj] * ctrl_w[ll];
      d2N_dst[ll] = dNs_ds[ii] * dNt_dt[jj] * ctrl_w[ll];

      dw_ds += dN_ds[ll];
      dw_dt += dN_dt[ll];
      d2w_dss += d2N_dss[ll];
      d2w_dtt += d2N_dtt[ll];
      d2w_dst += d2N_dst[ll];
    }
  }

  const double inv_w = 1.0 / w;
  const int sindex = nLocBas * quaindex;

  double m0 = 0.0; double m1 = 0.0; double m2 = 0.0; double m3 = 0.0;
 
  double d2x_dss = 0.0; double d2x_dtt = 0.0; double d2x_dst = 0.0;
  double d2y_dss = 0.0; double d2y_dtt = 0.0; double d2y_dst = 0.0; 

  for(ii=0; ii<nLocBas; ++ii)
  {
    R[sindex+ii] = N[ii] * inv_w;

    dR_ds[ii] = (dN_ds[ii] - R[sindex+ii] * dw_ds) * inv_w; 
    dR_dt[ii] = (dN_dt[ii] - R[sindex+ii] * dw_dt) * inv_w; 
    
    d2R_dss[ii] = (d2N_dss[ii] - R[sindex+ii] * d2w_dss - 2.0 * dw_ds * dR_ds[ii]) * inv_w;
    d2R_dtt[ii] = (d2N_dtt[ii] - R[sindex+ii] * d2w_dtt - 2.0 * dw_dt * dR_dt[ii]) * inv_w;
    d2R_dst[ii] = (d2N_dst[ii] - R[sindex+ii]*d2w_dst - dw_dt*dR_ds[ii] - dw_ds*dR_dt[ii])*inv_w;

    m0 += ctrl_x[ii] * dR_ds[ii];
    m1 += ctrl_x[ii] * dR_dt[ii];
    m2 += ctrl_y[ii] * dR_ds[ii];
    m3 += ctrl_y[ii] * dR_dt[ii];

    d2x_dss += ctrl_x[ii] * d2R_dss[ii];
    d2x_dtt += ctrl_x[ii] * d2R_dtt[ii];
    d2x_dst += ctrl_x[ii] * d2R_dst[ii];
    d2y_dss += ctrl_y[ii] * d2R_dss[ii];
    d2y_dtt += ctrl_y[ii] * d2R_dtt[ii];
    d2y_dst += ctrl_y[ii] * d2R_dst[ii];
  }
  
  detJac_temp = m0 * m3 - m1 * m2;

  detJac[quaindex] = detJac_temp * hx * hy;

  const double inv_detJactmp = 1.0 / detJac_temp;

  const double ds_dx = m3 * inv_detJactmp;
  const double ds_dy = -1.0 * m1 * inv_detJactmp;
  const double dt_dx = -1.0 * m2 * inv_detJactmp;
  const double dt_dy = m0 * inv_detJactmp;

  for(ii=0; ii<nLocBas; ++ii)
  {
    dR_dx[sindex+ii] = dR_ds[ii] * ds_dx + dR_dt[ii] * dt_dx;
    dR_dy[sindex+ii] = dR_ds[ii] * ds_dy + dR_dt[ii] * dt_dy;
  }

  double Rx, Ry, rhs1, rhs2, rhs3;
  const double inv_bottom = inv_detJactmp * inv_detJactmp;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    Rx = dR_dx[sindex + ii]; Ry = dR_dy[sindex + ii];
    rhs1 = d2R_dss[ii] - Rx * d2x_dss - Ry * d2y_dss;
    rhs2 = d2R_dst[ii] - Rx * d2x_dst - Ry * d2y_dst;
    rhs3 = d2R_dtt[ii] - Rx * d2x_dtt - Ry * d2y_dtt;

    d2R_dxx[sindex+ii] = (m3*m3*rhs1 - 2.0*m2*m3*rhs2 + m2*m2*rhs3) * inv_bottom;
    d2R_dxy[sindex+ii] = (-1.0*m1*m3*rhs1 + (m0*m3+m1*m2)*rhs2 - m0*m2*rhs3) * inv_bottom;
    d2R_dyy[sindex+ii] = (m1*m1*rhs1 - 2.0*m0*m1*rhs2 + m0*m0*rhs3) * inv_bottom;
  }
}

//EOF
