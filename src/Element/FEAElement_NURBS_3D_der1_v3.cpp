#include "FEAElement_NURBS_3D_der1_v3.hpp"

FEAElement_NURBS_3D_der1_v3::FEAElement_NURBS_3D_der1_v3(const int &in_eIndex,
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

FEAElement_NURBS_3D_der1_v3::~FEAElement_NURBS_3D_der1_v3()
{
  clearBasisCache();
}

void FEAElement_NURBS_3D_der1_v3::clearBasisCache()
{
  if(is_sNonzero)
  {
    delete [] R; 
    delete [] dR_dx; 
    delete [] dR_dy; 
    delete [] dR_dz;
    delete [] detJac;
  }
	R = NULL; dR_dx = NULL; dR_dy = NULL; dR_dz = NULL; detJac = NULL;
}

void FEAElement_NURBS_3D_der1_v3::buildBasis(
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
  detJac = new double [numQuapts];

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

  // Initialize three-dimensional NURBS functions
  double * dRr_ds = new double [nLocBas];
  double * dRr_dt = new double [nLocBas];
  double * dRr_du = new double [nLocBas];

  // Store the quadrature info of univariate B-spline 
  double * Nns = new double [ssdp1];
  double * Nnt = new double [ttdp1];
  double * Nnu = new double [uudp1];
  double * dNns_ds = new double [ssdp1];
  double * dNnt_dt = new double [ttdp1];
  double * dNnu_du = new double [uudp1];

  // Store the bernstein polynomials
  double * Bbs = new double [ssdp1];
  double * dBbs_ds = new double [ssdp1];
  double * Bbt = new double [ttdp1];
  double * dBbt_dt = new double [ttdp1];
  double * Bbu = new double [uudp1];
  double * dBbu_du = new double [uudp1];

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
      Bbu[ii] = Bu->get_der0(ii, qua_z);
      dBbu_du[ii] = Bu->get_der1(ii, qua_z) * invhz;
    }
    ll = -1;
    for(ii=0; ii<uudp1; ++ii)
    {
      Nnu[ii] = 0.0; dNnu_du[ii] = 0.0;
      for(jj=0; jj<uudp1; ++jj)
      {
        ll += 1;
        extor = ext_z[ll];
        Nnu[ii]     += extor * Bbu[jj]; 
        dNnu_du[ii] += extor * dBbu_du[jj]; 
      }
    }
    for(int qua_y = 0; qua_y<num_qt; ++qua_y)
    {
      for(ii=0; ii<ttdp1; ++ii)
      {
        Bbt[ii] = Bt->get_der0(ii,qua_y);
        dBbt_dt[ii] = Bt->get_der1(ii,qua_y) * invhy;
      }
      ll = -1;
      for(ii=0; ii<ttdp1; ++ii)
      {
        Nnt[ii] = 0.0; dNnt_dt[ii] = 0.0;
        for(jj=0; jj<ttdp1; ++jj)
        {
          ll += 1;
          extor = ext_y[ll];
          Nnt[ii]     += extor * Bbt[jj]; 
          dNnt_dt[ii] += extor * dBbt_dt[jj]; 
        }
      }
      for(int qua_x = 0; qua_x<num_qs; ++qua_x)
      {
        for(ii=0; ii<ssdp1; ++ii)
        {
          Bbs[ii] = Bs->get_der0(ii,qua_x);
          dBbs_ds[ii] = Bs->get_der1(ii, qua_x) * invhx;
        }
        ll = -1;
        for(ii=0; ii<ssdp1; ++ii)
        {
          Nns[ii] = 0.0; dNns_ds[ii] = 0.0;
          for(jj=0; jj<ssdp1; ++jj)
          {
            ll += 1;
            extor = ext_x[ll];
            Nns[ii]     += extor * Bbs[jj]; 
            dNns_ds[ii] += extor * dBbs_ds[jj]; 
          }
        }
        counter += 1; 
        BuildShape_atQua( ssdp1, ttdp1, uudp1, counter, qua_x, qua_y, qua_z,
            hx, hy, hz,
            ctrl_x, ctrl_y, ctrl_z, ctrl_w, Nn, dNn_ds, dNn_dt, dNn_du,
            dRr_ds, dRr_dt, dRr_du, Nns, Nnt, Nnu, dNns_ds, dNnt_dt, dNnu_du);
      }
    }
  }

  delete [] ctrl_x; delete [] ctrl_y; delete [] ctrl_z; delete [] ctrl_w;
  delete [] ext_x; delete [] ext_y; delete [] ext_z;
  delete [] Nn; delete [] dNn_ds; delete [] dNn_dt; delete [] dNn_du;
  delete [] dRr_ds; delete [] dRr_dt; delete [] dRr_du;
  delete [] Nns; delete [] Nnt; delete [] Nnu;
  delete [] dNns_ds; delete [] dNnt_dt; delete [] dNnu_du;
  delete [] Bbs; delete [] Bbt; delete [] Bbu;
  delete [] dBbs_ds; delete [] dBbt_dt; delete [] dBbu_du;
}


void FEAElement_NURBS_3D_der1_v3::get_R(const int &quaindex, 
    double * const &basis) const
{
  for(int ii=0; ii<nLocBas; ++ii)
    basis[ii] = R[quaindex * nLocBas + ii];
}

void FEAElement_NURBS_3D_der1_v3::get_gradR(const int &quaindex, 
    double * const &basis_x, double * const &basis_y, double * const &basis_z) const
{
  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis_x[ii] = dR_dx[quaindex * nLocBas + ii];
    basis_y[ii] = dR_dy[quaindex * nLocBas + ii];
    basis_z[ii] = dR_dz[quaindex * nLocBas + ii];
  }
}

void FEAElement_NURBS_3D_der1_v3::get_R_gradR(const int &quaindex, 
    double * const &basis, double * const &basis_x, 
    double * const &basis_y, double * const &basis_z) const
{
  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis[ii]   = R[quaindex * nLocBas + ii];
    basis_x[ii] = dR_dx[quaindex * nLocBas + ii];
    basis_y[ii] = dR_dy[quaindex * nLocBas + ii];
    basis_z[ii] = dR_dz[quaindex * nLocBas + ii];
  }
}

double FEAElement_NURBS_3D_der1_v3::get_detJac( const int &quaindex ) const
{
  return detJac[quaindex];
}

void FEAElement_NURBS_3D_der1_v3::print() const
{
  SYS_T::commPrint("NURBS_3D_der1_v3: ");
  SYS_T::commPrint("3-Dimensional NURBS shape function with 1st order derivatives. \n");
  SYS_T::commPrint("Note: Jacobian matrix is not evaluated. \n");
}

double FEAElement_NURBS_3D_der1_v3::get_memory_usage() const
{
  unsigned int double_size = 0; unsigned int int_size = 3;
  //double_size += R.size(); double_size += dR_dx.size(); double_size += dR_dy.size();
  //double_size += dR_dz.size(); double_size += detJac.size();
  double_size = nLocBas * numQuapts * 4 + numQuapts;

  return (double) double_size * 8.0 + (double) int_size * 4.0;
}

void FEAElement_NURBS_3D_der1_v3::BuildShape_atQua(
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
    double * &dR_ds, double * &dR_dt, double * &dR_du,
    const double * const &Ns, 
    const double * const &Nt, 
    const double * const &Nu,
    const double * const &dNs_ds,
    const double * const &dNt_dt,
    const double * const &dNu_du
    )
{
  // iterator index
  int ii, jj, kk, ll;

  double detJac_temp = 0.0; 

  double w = 0.0; // weight at this qua point
  double dw_ds = 0.0; double dw_dt = 0.0; double dw_du = 0.0;

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

        dw_ds += dN_ds[ll];
        dw_dt += dN_dt[ll];
        dw_du += dN_du[ll];
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

  for(ii=0; ii<nLocBas; ++ii)
  {
    R[sindex+ii] = N[ii] * inv_w;
    dR_ds[ii] = (dN_ds[ii] - R[sindex+ii] * dw_ds) * inv_w;
    dR_dt[ii] = (dN_dt[ii] - R[sindex+ii] * dw_dt) * inv_w;
    dR_du[ii] = (dN_du[ii] - R[sindex+ii] * dw_du) * inv_w;

    m0 += ctrl_x[ii] * dR_ds[ii];   // dx_dxi
    m1 += ctrl_x[ii] * dR_dt[ii];  // dx_deta
    m2 += ctrl_x[ii] * dR_du[ii]; // dx_dzeta

    m3 += ctrl_y[ii] * dR_ds[ii];   // dy_dxi
    m4 += ctrl_y[ii] * dR_dt[ii];  // dy_deta
    m5 += ctrl_y[ii] * dR_du[ii]; // dy_dzeta

    m6 += ctrl_z[ii] * dR_ds[ii];   // dz_dxi
    m7 += ctrl_z[ii] * dR_dt[ii];  //dz_deta
    m8 += ctrl_z[ii] * dR_du[ii]; // dz_dzeta
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

  // Now compute the derivatives w.r.t phyiscal coordinates
  for(ii=0; ii<nLocBas; ++ii)
  {
    dR_dx[sindex+ii] = dR_ds[ii] * ds_dx + dR_dt[ii] * dt_dx + dR_du[ii] * du_dx;
    dR_dy[sindex+ii] = dR_ds[ii] * ds_dy + dR_dt[ii] * dt_dy + dR_du[ii] * du_dy;
    dR_dz[sindex+ii] = dR_ds[ii] * ds_dz + dR_dt[ii] * dt_dz + dR_du[ii] * du_dz;
  }
}
