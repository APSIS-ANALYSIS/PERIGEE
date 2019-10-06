#include "FEAElement_NURBS_3D_der2_vms.hpp"

FEAElement_NURBS_3D_der2_vms::FEAElement_NURBS_3D_der2_vms(const int &in_eIndex,
    const IALocal_meshSize * const &mSize, const BernsteinBasis_Array * const &Bs,
    const BernsteinBasis_Array * const &Bt, const BernsteinBasis_Array * const &Bu,
    const FEANode * const &feaNode, const IAExtractor * const &extractor,
    const ALocal_IEN * const &locIEN)
{
  // Assign values to single data, initialize quadrature arrays.
  elem_index = in_eIndex;
	
  int sDegree = Bs->get_degree();
	int tDegree = Bt->get_degree();
	int uDegree = Bu->get_degree();

	int sDp1 = sDegree + 1; int tDp1 = tDegree + 1; int uDp1 = uDegree + 1;

  nLocBas = sDp1 * tDp1 * uDp1;

  int num_qua_s = Bs->get_nQuapts();
  int num_qua_t = Bt->get_nQuapts();
  int num_qua_u = Bu->get_nQuapts();

  numQuapts = num_qua_s * num_qua_t * num_qua_u;

  is_sNonzero = true;
  double refSize = mSize->get_hx(elem_index) * mSize->get_hy(elem_index) * 
    mSize->get_hz(elem_index);

  if(refSize == 0.0)
    is_sNonzero = false;

	R = NULL; dR_dx = NULL; dR_dy = NULL; dR_dz = NULL; detJac = NULL;
  d2R_dxx = NULL; d2R_dyy = NULL; d2R_dzz = NULL;
  d2R_dxy = NULL; d2R_dxz = NULL; d2R_dyz = NULL;
  
  ds_dx_qua = NULL; ds_dy_qua = NULL; ds_dz_qua = NULL;
  dt_dx_qua = NULL; dt_dy_qua = NULL; dt_dz_qua = NULL;
  du_dx_qua = NULL; du_dy_qua = NULL; du_dz_qua = NULL;
  
  // Calculate the quadrature of basis functions
  // memory layout of basis function vectors R, dR_dx, dR_dy, dR_dz:
  // (basis0 qua0) ... (basis_nLocBas-1, qua 0) ... (basis 0, qua1) ...
  // detJac: qua0, qua1, ..., qua_numQuapts-1
  // Jacobian: dx/dxi dx/deta dx/dzeta dy/dxi dy/deta dy/dzeta 
  //           dz/dxi dz/deta dz/dzeta at qua 0, ... at qua_numQuapts-1
  if(is_sNonzero == true)
  {
    buildBasis( sDp1, tDp1, uDp1, num_qua_s, num_qua_t, num_qua_u,
        mSize, Bs, Bt, Bu, feaNode, extractor, locIEN);
  }
}

FEAElement_NURBS_3D_der2_vms::~FEAElement_NURBS_3D_der2_vms()
{
  clearBasisCache();
}

void FEAElement_NURBS_3D_der2_vms::clearBasisCache()
{
  if( is_sNonzero )
  {
    delete [] R; 
    delete [] dR_dx;   delete [] dR_dy;   delete [] dR_dz;
    delete [] d2R_dxx; delete [] d2R_dyy; delete [] d2R_dzz;
    delete [] d2R_dxy; delete [] d2R_dxz; delete [] d2R_dyz;
    delete [] detJac;
    delete [] ds_dx_qua; delete [] ds_dy_qua; delete [] ds_dz_qua;
    delete [] dt_dx_qua; delete [] dt_dy_qua; delete [] dt_dz_qua;
    delete [] du_dx_qua; delete [] du_dy_qua; delete [] du_dz_qua;
  }
  R = NULL; dR_dx = NULL; dR_dy = NULL; dR_dz = NULL; detJac = NULL;
  d2R_dxx = NULL; d2R_dyy = NULL; d2R_dzz = NULL;
  d2R_dxy = NULL; d2R_dxz = NULL; d2R_dyz = NULL;
  ds_dx_qua = NULL; ds_dy_qua = NULL; ds_dz_qua = NULL;
  dt_dx_qua = NULL; dt_dy_qua = NULL; dt_dz_qua = NULL;
  du_dx_qua = NULL; du_dy_qua = NULL; du_dz_qua = NULL;
}

void FEAElement_NURBS_3D_der2_vms::buildBasis(
    const int &ssdp1, const int &ttdp1, const int &uudp1,
    const int &num_qs, const int &num_qt, const int &num_qu,
    const IALocal_meshSize * const &mSize, const BernsteinBasis_Array * const &Bs,
    const BernsteinBasis_Array * const &Bt, const BernsteinBasis_Array * const &Bu,
    const FEANode * const &feaNode, const IAExtractor * const &extractor,
    const ALocal_IEN * const &locIEN )
{
  // Allocate memory
  int basis_size = nLocBas * numQuapts;

  R = new double [basis_size];
  dR_dx = new double [basis_size];
  dR_dy = new double [basis_size];
  dR_dz = new double [basis_size];
  d2R_dxx = new double [basis_size];
  d2R_dyy = new double [basis_size];
  d2R_dzz = new double [basis_size];
  d2R_dxy = new double [basis_size];
  d2R_dxz = new double [basis_size];
  d2R_dyz = new double [basis_size];
  detJac  = new double [numQuapts];
  
  ds_dx_qua = new double [numQuapts];
  ds_dy_qua = new double [numQuapts];
  ds_dz_qua = new double [numQuapts];
  dt_dx_qua = new double [numQuapts];
  dt_dy_qua = new double [numQuapts];
  dt_dz_qua = new double [numQuapts];
  du_dx_qua = new double [numQuapts];
  du_dy_qua = new double [numQuapts];
  du_dz_qua = new double [numQuapts];


  // Assign necessary data
  double hx = mSize->get_hx(elem_index);
  double hy = mSize->get_hy(elem_index);
  double hz = mSize->get_hz(elem_index);

  double * ext_x; double * ext_y; double * ext_z;
  extractor->get_EXT_x(elem_index, ext_x);
  extractor->get_EXT_y(elem_index, ext_y);
  extractor->get_EXT_z(elem_index, ext_z);

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
  
  // Initialize three-dimensional b-spline functions
  double * Nn = new double [nLocBas];
  double * dNn_ds = new double [nLocBas];
  double * dNn_dt = new double [nLocBas];
  double * dNn_du = new double [nLocBas];
  
  double * d2Nn_dss = new double [nLocBas];
  double * d2Nn_dtt = new double [nLocBas];
  double * d2Nn_duu = new double [nLocBas];
  double * d2Nn_dst = new double [nLocBas];
  double * d2Nn_dsu = new double [nLocBas];
  double * d2Nn_dtu = new double [nLocBas];

  // Initialize three-dimensional NURBS functions
  double * dRr_ds = new double [nLocBas];
  double * dRr_dt = new double [nLocBas];
  double * dRr_du = new double [nLocBas];
  
  double * d2Rr_dss = new double [nLocBas];
  double * d2Rr_dtt = new double [nLocBas];
  double * d2Rr_duu = new double [nLocBas];
  double * d2Rr_dst = new double [nLocBas];
  double * d2Rr_dsu = new double [nLocBas];
  double * d2Rr_dtu = new double [nLocBas];

  // Store the quadrature info of univariate B-spline 
  double * Nns = new double [ssdp1];
  double * Nnt = new double [ttdp1];
  double * Nnu = new double [uudp1];
  double * dNns_ds = new double [ssdp1];
  double * dNnt_dt = new double [ttdp1];
  double * dNnu_du = new double [uudp1];

  double * d2Nns_dss = new double [ssdp1];
  double * d2Nnt_dtt = new double [ttdp1];
  double * d2Nnu_duu = new double [uudp1];

  // Store the bernstein polynomials
  double * Bbs       = new double [ssdp1];
  double * dBbs_ds   = new double [ssdp1];
  double * d2Bbs_dss = new double [ssdp1];
  double * Bbt       = new double [ttdp1];
  double * dBbt_dt   = new double [ttdp1];
  double * d2Bbt_dtt = new double [ttdp1];
  double * Bbu       = new double [uudp1];
  double * dBbu_du   = new double [uudp1];
  double * d2Bbu_duu = new double [uudp1];
  
  // Allocation for 6-by-6 right hand side
  double * der2rhs = new double [6];
  double * der2sol = new double [6];


  // This function will calculate basis functions at the quadrature point
  // with index qua. And it will fill R, dR_dx, ... at positions
  // [nLocBas * qua, nLocBas * (qua+1) -1]
  int ii, jj, ll; // iterator
  int counter = -1;
  double invhx = 1.0 / hx; double invhy = 1.0 / hy; double invhz = 1.0 / hz;
  double extor;
  for(int qua_z = 0; qua_z<num_qu; ++qua_z)
  {
    for(ii=0; ii<uudp1; ++ii)
    {
      Bbu[ii]       = Bu->get_der0(ii, qua_z);
      dBbu_du[ii]   = Bu->get_der1(ii, qua_z) * invhz;
      d2Bbu_duu[ii] = Bu->get_der2(ii, qua_z) * invhz * invhz;
    }
    ll = -1;
    for(ii=0; ii<uudp1; ++ii)
    {
      Nnu[ii] = 0.0; dNnu_du[ii] = 0.0; d2Nnu_duu[ii] = 0.0;
      for(jj=0; jj<uudp1; ++jj)
      {
        ll += 1;
        extor = ext_z[ll];
        Nnu[ii]       += extor * Bbu[jj]; 
        dNnu_du[ii]   += extor * dBbu_du[jj];
        d2Nnu_duu[ii] += extor * d2Bbu_duu[jj]; 
      }
    }
    for(int qua_y = 0; qua_y<num_qt; ++qua_y)
    {
      for(ii=0; ii<ttdp1; ++ii)
      {
        Bbt[ii]       = Bt->get_der0(ii,qua_y);
        dBbt_dt[ii]   = Bt->get_der1(ii,qua_y) * invhy;
        d2Bbt_dtt[ii] = Bt->get_der2(ii,qua_y) * invhy * invhy;
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
      for(int qua_x = 0; qua_x<num_qs; ++qua_x)
      {
        for(ii=0; ii<ssdp1; ++ii)
        {
          Bbs[ii]       = Bs->get_der0(ii, qua_x);
          dBbs_ds[ii]   = Bs->get_der1(ii, qua_x) * invhx;
          d2Bbs_dss[ii] = Bs->get_der2(ii, qua_x) * invhx * invhx;
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
        BuildShape_atQua( ssdp1, ttdp1, uudp1, counter, qua_x, qua_y, qua_z,
            hx, hy, hz, ctrl_x, ctrl_y, ctrl_z, ctrl_w, 
            Nn, dNn_ds, dNn_dt, dNn_du, d2Nn_dss, d2Nn_dtt, d2Nn_duu,
            d2Nn_dst, d2Nn_dsu, d2Nn_dtu, dRr_ds, dRr_dt, dRr_du,
            d2Rr_dss, d2Rr_dtt, d2Rr_duu, d2Rr_dst, d2Rr_dsu, d2Rr_dtu,
            der2rhs, der2sol,
            Nns, Nnt, Nnu, dNns_ds, dNnt_dt, dNnu_du,
            d2Nns_dss, d2Nnt_dtt, d2Nnu_duu );
      }
    }
  }

  delete [] ctrl_x; delete [] ctrl_y; delete [] ctrl_z; delete [] ctrl_w;
  delete [] ext_x; delete [] ext_y; delete [] ext_z;
  delete [] Nn; delete [] dNn_ds; delete [] dNn_dt; delete [] dNn_du;
  delete [] d2Nn_dss; delete [] d2Nn_dtt; delete [] d2Nn_duu;
  delete [] d2Nn_dst; delete [] d2Nn_dsu; delete [] d2Nn_dtu;
  delete [] dRr_ds; delete [] dRr_dt; delete [] dRr_du;
  delete [] d2Rr_dss; delete [] d2Rr_dtt; delete [] d2Rr_duu;
  delete [] d2Rr_dst; delete [] d2Rr_dsu; delete [] d2Rr_dtu;
  delete [] Nns; delete [] Nnt; delete [] Nnu;
  delete [] dNns_ds; delete [] dNnt_dt; delete [] dNnu_du;
  delete [] d2Nns_dss; delete [] d2Nnt_dtt; delete [] d2Nnu_duu;
  delete [] Bbs; delete [] Bbt; delete [] Bbu;
  delete [] dBbs_ds; delete [] dBbt_dt; delete [] dBbu_du;
  delete [] d2Bbs_dss; delete [] d2Bbt_dtt; delete [] d2Bbu_duu;
  delete [] der2rhs; delete [] der2sol;
}


void FEAElement_NURBS_3D_der2_vms::get_R(const int &quaindex, 
    double * const &basis) const
{
  for(int ii=0; ii<nLocBas; ++ii)
    basis[ii] = R[quaindex * nLocBas + ii];
}

void FEAElement_NURBS_3D_der2_vms::get_gradR(const int &quaindex, 
    double * const &basis_x, double * const &basis_y, double * const &basis_z) const
{
  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis_x[ii] = dR_dx[quaindex * nLocBas + ii];
    basis_y[ii] = dR_dy[quaindex * nLocBas + ii];
    basis_z[ii] = dR_dz[quaindex * nLocBas + ii];
  }
}


void FEAElement_NURBS_3D_der2_vms::get_3D_R_gradR_LaplacianR(
    const int &quaindex, double * const &basis,
    double * const &basis_x, double * const &basis_y,
    double * const &basis_z, double * const &basis_xx, 
    double * const &basis_yy, double * const &basis_zz ) const
{
  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis[ii] = R[quaindex * nLocBas + ii];
    basis_x[ii] = dR_dx[quaindex * nLocBas + ii];
    basis_y[ii] = dR_dy[quaindex * nLocBas + ii];
    basis_z[ii] = dR_dz[quaindex * nLocBas + ii];
    basis_xx[ii] = d2R_dxx[quaindex * nLocBas + ii];
    basis_yy[ii] = d2R_dyy[quaindex * nLocBas + ii];
    basis_zz[ii] = d2R_dzz[quaindex * nLocBas + ii];
  }
}


void FEAElement_NURBS_3D_der2_vms::get_invJacobian(const int &quaindex,
    double * const &dxi_dx) const
{
  dxi_dx[0] = ds_dx_qua[quaindex];
  dxi_dx[1] = ds_dy_qua[quaindex];
  dxi_dx[2] = ds_dz_qua[quaindex];
  dxi_dx[3] = dt_dx_qua[quaindex];
  dxi_dx[4] = dt_dy_qua[quaindex];
  dxi_dx[5] = dt_dz_qua[quaindex];
  dxi_dx[6] = du_dx_qua[quaindex];
  dxi_dx[7] = du_dy_qua[quaindex];
  dxi_dx[8] = du_dz_qua[quaindex];
}


void FEAElement_NURBS_3D_der2_vms::print() const
{
  SYS_T::commPrint("NURBS_3D_der2_vms: ");
  SYS_T::commPrint("3-Dimensional NURBS shape function with 2nd order derivatives. \n");
  SYS_T::commPrint("Note: Jacobian matrix is not evaluated. \n");
  SYS_T::commPrint("      The 9 components of dxi_dx are evaluated at qua points. \n");
}

double FEAElement_NURBS_3D_der2_vms::get_memory_usage() const
{
  unsigned int double_size = 0; unsigned int int_size = 3;
  double_size = nLocBas * numQuapts * 10 + numQuapts * 10;

  return (double) double_size * 8.0 + (double) int_size * 4.0;
}

void FEAElement_NURBS_3D_der2_vms::BuildShape_atQua(
    const int &sDp1, const int &tDp1, const int &uDp1,
    const int &quaindex,
    const int &quaindex_x, 
    const int &quaindex_y, 
    const int &quaindex_z, 
    const double &hx, const double &hy, const double &hz,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z,
    const double * const &ctrl_w,
    double * &N, double * &dN_ds, double * &dN_dt, double * &dN_du,
    double * &d2N_dss, double * &d2N_dtt, double * &d2N_duu,
    double * &d2N_dst, double * &d2N_dsu, double * &d2N_dtu,
    double * &dR_ds, double * &dR_dt, double * &dR_du,
    double * &d2R_dss, double * &d2R_dtt, double * &d2R_duu,
    double * &d2R_dst, double * &d2R_dsu, double * &d2R_dtu,
    double * &RHS, double * &sol,
    const double * const &Ns, 
    const double * const &Nt, 
    const double * const &Nu,
    const double * const &dNs_ds,
    const double * const &dNt_dt,
    const double * const &dNu_du,
    const double * const &d2Ns_dss,
    const double * const &d2Nt_dtt,
    const double * const &d2Nu_duu
    )
{
  // iterator index
  int ii, jj, kk, ll;

  double detJac_temp = 0.0; 

  double w = 0.0; // weight at this qua point
  double dw_ds = 0.0;   double dw_dt = 0.0;   double dw_du = 0.0;
  double d2w_dss = 0.0; double d2w_dtt = 0.0; double d2w_duu = 0.0;
  double d2w_dst = 0.0; double d2w_dsu = 0.0; double d2w_dtu = 0.0;

  ll = -1;
  for(kk=0; kk<uDp1; ++kk)
  {
    for(jj=0; jj<tDp1; ++jj)
    {
      for(ii=0; ii<sDp1; ++ii)
      {
        ll += 1;
        N[ll] = Ns[ii] * Nt[jj] * Nu[kk] * ctrl_w[ll];
        w += N[ll];

        dN_ds[ll] = dNs_ds[ii] * Nt[jj] * Nu[kk] * ctrl_w[ll];
        dN_dt[ll] = Ns[ii] * dNt_dt[jj] * Nu[kk] * ctrl_w[ll];
        dN_du[ll] = Ns[ii] * Nt[jj] * dNu_du[kk] * ctrl_w[ll];

        d2N_dss[ll] = d2Ns_dss[ii] * Nt[jj] * Nu[kk] * ctrl_w[ll];
        d2N_dtt[ll] = Ns[ii] * d2Nt_dtt[jj] * Nu[kk] * ctrl_w[ll];
        d2N_duu[ll] = Ns[ii] * Nt[jj] * d2Nu_duu[kk] * ctrl_w[ll];
        d2N_dst[ll] = dNs_ds[ii] * dNt_dt[jj] * Nu[kk] * ctrl_w[ll];
        d2N_dsu[ll] = dNs_ds[ii] * Nt[jj] * dNu_du[kk] * ctrl_w[ll];
        d2N_dtu[ll] = Ns[ii] * dNt_dt[jj] * dNu_du[kk] * ctrl_w[ll];

        dw_ds += dN_ds[ll];
        dw_dt += dN_dt[ll];
        dw_du += dN_du[ll];

        d2w_dss += d2N_dss[ll];
        d2w_dtt += d2N_dtt[ll];
        d2w_duu += d2N_duu[ll];
        d2w_dst += d2N_dst[ll];
        d2w_dsu += d2N_dsu[ll];
        d2w_dtu += d2N_dtu[ll];
      }
    }
  }

  double inv_w = 1.0 / w;
  int sindex = nLocBas * quaindex;  

  // Compute the geometry mapping and Jacobian
  //vector<double> vecJac(9, 0.0);
  double m0 = 0.0; double m1 = 0.0; double m2 = 0.0;
  double m3 = 0.0; double m4 = 0.0; double m5 = 0.0;
  double m6 = 0.0; double m7 = 0.0; double m8 = 0.0;

  double d2x_dss = 0.0; double d2x_dtt = 0.0; double d2x_duu = 0.0;
  double d2x_dst = 0.0; double d2x_dsu = 0.0; double d2x_dtu = 0.0;

  double d2y_dss = 0.0; double d2y_dtt = 0.0; double d2y_duu = 0.0;
  double d2y_dst = 0.0; double d2y_dsu = 0.0; double d2y_dtu = 0.0;

  double d2z_dss = 0.0; double d2z_dtt = 0.0; double d2z_duu = 0.0;
  double d2z_dst = 0.0; double d2z_dsu = 0.0; double d2z_dtu = 0.0;

  for(ii=0; ii<nLocBas; ++ii)
  {
    R[sindex+ii] = N[ii] * inv_w;
    dR_ds[ii] = (dN_ds[ii] - R[sindex+ii] * dw_ds) * inv_w;
    dR_dt[ii] = (dN_dt[ii] - R[sindex+ii] * dw_dt) * inv_w;
    dR_du[ii] = (dN_du[ii] - R[sindex+ii] * dw_du) * inv_w;

    d2R_dss[ii] = (d2N_dss[ii] - 2.0 * dR_ds[ii] * dw_ds - R[sindex+ii] * d2w_dss ) * inv_w;
    d2R_dtt[ii] = (d2N_dtt[ii] - 2.0 * dR_dt[ii] * dw_dt - R[sindex+ii] * d2w_dtt ) * inv_w;
    d2R_duu[ii] = (d2N_duu[ii] - 2.0 * dR_du[ii] * dw_du - R[sindex+ii] * d2w_duu ) * inv_w;
    d2R_dst[ii] = (d2N_dst[ii] - dR_dt[ii]*dw_ds - dR_ds[ii]*dw_dt - R[sindex+ii]*d2w_dst) * inv_w;
    d2R_dsu[ii] = (d2N_dsu[ii] - dR_du[ii]*dw_ds - dR_ds[ii]*dw_du - R[sindex+ii]*d2w_dsu) * inv_w;
    d2R_dtu[ii] = (d2N_dtu[ii] - dR_du[ii]*dw_dt - dR_dt[ii]*dw_du - R[sindex+ii]*d2w_dtu) * inv_w;


    m0 += ctrl_x[ii] * dR_ds[ii]; // dx_ds
    m1 += ctrl_x[ii] * dR_dt[ii]; // dx_dt
    m2 += ctrl_x[ii] * dR_du[ii]; // dx_du

    m3 += ctrl_y[ii] * dR_ds[ii]; // dy_ds
    m4 += ctrl_y[ii] * dR_dt[ii]; // dy_dt
    m5 += ctrl_y[ii] * dR_du[ii]; // dy_du

    m6 += ctrl_z[ii] * dR_ds[ii]; // dz_ds
    m7 += ctrl_z[ii] * dR_dt[ii]; // dz_dt
    m8 += ctrl_z[ii] * dR_du[ii]; // dz_du

    d2x_dss += ctrl_x[ii] * d2R_dss[ii];
    d2x_dtt += ctrl_x[ii] * d2R_dtt[ii];
    d2x_duu += ctrl_x[ii] * d2R_duu[ii];
    d2x_dst += ctrl_x[ii] * d2R_dst[ii];
    d2x_dsu += ctrl_x[ii] * d2R_dsu[ii];
    d2x_dtu += ctrl_x[ii] * d2R_dtu[ii];

    d2y_dss += ctrl_y[ii] * d2R_dss[ii];
    d2y_dtt += ctrl_y[ii] * d2R_dtt[ii];
    d2y_duu += ctrl_y[ii] * d2R_duu[ii];
    d2y_dst += ctrl_y[ii] * d2R_dst[ii];
    d2y_dsu += ctrl_y[ii] * d2R_dsu[ii];
    d2y_dtu += ctrl_y[ii] * d2R_dtu[ii];

    d2z_dss += ctrl_z[ii] * d2R_dss[ii];
    d2z_dtt += ctrl_z[ii] * d2R_dtt[ii];
    d2z_duu += ctrl_z[ii] * d2R_duu[ii];
    d2z_dst += ctrl_z[ii] * d2R_dst[ii];
    d2z_dsu += ctrl_z[ii] * d2R_dsu[ii];
    d2z_dtu += ctrl_z[ii] * d2R_dtu[ii];
  } 

  // get determinant of jacobian
  detJac_temp = m0 * m4 * m8 + m1 * m5 * m6 + m2 * m3 * m7
    - m2 * m4 * m6 - m0 * m5 * m7 - m1 * m3 * m8 ;
  detJac[quaindex] = detJac_temp * hx * hy * hz;

  // Jacobian invertion get dxi/dx, dxi/dy, dxi/dz; deta/dx, deta/dy,
  // deta/dz; dzeta/dx, dzeta/dy, dzeta/dz.
  double inv_detJactmp = 1.0 / detJac_temp;
  double ds_dx = (m4*m8 - m5*m7) * inv_detJactmp;
  double ds_dy = (m2*m7 - m1*m8) * inv_detJactmp;
  double ds_dz = (m1*m5 - m2*m4) * inv_detJactmp;
  double dt_dx = (m5*m6 - m3*m8) * inv_detJactmp;
  double dt_dy = (m0*m8 - m2*m6) * inv_detJactmp;
  double dt_dz = (m2*m3 - m0*m5) * inv_detJactmp;
  double du_dx = (m3*m7 - m4*m6) * inv_detJactmp;
  double du_dy = (m1*m6 - m0*m7) * inv_detJactmp;
  double du_dz = (m0*m4 - m1*m3) * inv_detJactmp;

  // Now save the inversion of the Jacobian in the nine arrays
  ds_dx_qua[quaindex] = ds_dx;
  ds_dy_qua[quaindex] = ds_dy;
  ds_dz_qua[quaindex] = ds_dz;
  dt_dx_qua[quaindex] = dt_dx;
  dt_dy_qua[quaindex] = dt_dy;
  dt_dz_qua[quaindex] = dt_dz;
  du_dx_qua[quaindex] = du_dx;
  du_dy_qua[quaindex] = du_dy;
  du_dz_qua[quaindex] = du_dz;

  // Now compute the derivatives w.r.t phyiscal coordinates
  for(ii=0; ii<nLocBas; ++ii)
  {
    dR_dx[sindex+ii] = dR_ds[ii] * ds_dx + dR_dt[ii] * dt_dx + dR_du[ii] * du_dx;
    dR_dy[sindex+ii] = dR_ds[ii] * ds_dy + dR_dt[ii] * dt_dy + dR_du[ii] * du_dy;
    dR_dz[sindex+ii] = dR_ds[ii] * ds_dz + dR_dt[ii] * dt_dz + dR_du[ii] * du_dz;
  }

  // Calculate 2nd order derivatives by solving a 6-by-6 matrix problem
  Matrix_double_6by6_Array LHS(m0, m1, m2, m3, m4, m5, m6, m7, m8);

  // LU decomposition of this matrix problem
  LHS.LU_fac();

  double Rx, Ry, Rz;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    Rx = dR_dx[sindex+ii]; Ry = dR_dy[sindex+ii]; Rz = dR_dz[sindex+ii];
    RHS[0] = d2R_dss[ii] - Rx * d2x_dss - Ry * d2y_dss - Rz * d2z_dss;
    RHS[1] = d2R_dst[ii] - Rx * d2x_dst - Ry * d2y_dst - Rz * d2z_dst;
    RHS[2] = d2R_dsu[ii] - Rx * d2x_dsu - Ry * d2y_dsu - Rz * d2z_dsu;
    RHS[3] = d2R_dtt[ii] - Rx * d2x_dtt - Ry * d2y_dtt - Rz * d2z_dtt;
    RHS[4] = d2R_dtu[ii] - Rx * d2x_dtu - Ry * d2y_dtu - Rz * d2z_dtu;
    RHS[5] = d2R_duu[ii] - Rx * d2x_duu - Ry * d2y_duu - Rz * d2z_duu;

    LHS.LU_solve(RHS, sol);

    d2R_dxx[sindex+ii] = sol[0];
    d2R_dyy[sindex+ii] = sol[1];
    d2R_dzz[sindex+ii] = sol[2];
    d2R_dxy[sindex+ii] = sol[3];
    d2R_dxz[sindex+ii] = sol[4];
    d2R_dyz[sindex+ii] = sol[5];
  }
}

// EOF
