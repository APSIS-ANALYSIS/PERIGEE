#include "FEAElement_NURBS_2D_der1.hpp"

FEAElement_NURBS_2D_der1::FEAElement_NURBS_2D_der1( const int &in_eIndex,
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


FEAElement_NURBS_2D_der1::~FEAElement_NURBS_2D_der1()
{
  clearBasisCache();
}


void FEAElement_NURBS_2D_der1::clearBasisCache()
{
  if(is_sNonzero)
  {
    delete [] R;
    delete [] dR_dx;
    delete [] dR_dy;
    delete [] detJac;
  }
  R = NULL; dR_dx = NULL; dR_dy = NULL; detJac = NULL;
}


void FEAElement_NURBS_2D_der1::buildBasis( const int &ssdp1, const int &ttdp1,
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

  double * dRr_ds = new double [nLocBas];
  double * dRr_dt = new double [nLocBas];
  
  double * Nns = new double [ssdp1];
  double * Nnt = new double [ttdp1];

  double * dNns_ds = new double [ssdp1];
  double * dNnt_dt = new double [ttdp1];

  double * Bbs = new double [ssdp1];
  double * dBbs_ds = new double [ssdp1];
  double * Bbt = new double [ttdp1];
  double * dBbt_dt = new double [ttdp1];

  int ii, jj, ll;
  int counter = -1;
  const double invhx = 1.0 / hx; 
  const double invhy = 1.0 / hy;
  double extor;

  for(int qua_y=0; qua_y < num_qt; ++qua_y)
  {
    for(ii=0; ii<ttdp1; ++ii)
    {
      Bbt[ii]     = bt->get_der0(ii, qua_y);
      dBbt_dt[ii] = bt->get_der1(ii, qua_y) * invhy;
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
    
    for(int qua_x = 0; qua_x < num_qs; ++qua_x)
    {
      for(ii=0; ii<ssdp1; ++ii)
      {
        Bbs[ii] = bs->get_der0(ii, qua_x);
        dBbs_ds[ii] = bs->get_der1(ii, qua_x) * invhx;
      }
      ll = -1;
      for(ii=0; ii<ssdp1; ++ii)
      {
        Nns[ii] = 0.0; dNns_ds[ii] = 0.0;
        for(jj=0; jj<ssdp1; ++jj)
        {
          ll += 1;
          extor = ext_x[ll];
          Nns[ii] += extor * Bbs[jj];
          dNns_ds[ii] += extor * dBbs_ds[jj];
        }
      }
      counter += 1;

      BuildShape_atQua(  ssdp1, ttdp1, counter, qua_x, qua_y, hx, hy,
          ctrl_x, ctrl_y, ctrl_z, ctrl_w, Nn, dNn_ds, dNn_dt,
          dRr_ds, dRr_dt, Nns, Nnt, dNns_ds, dNnt_dt);
    } 
  }

  delete [] ctrl_x; delete [] ctrl_y; delete [] ctrl_z; delete [] ctrl_w;
  delete [] ext_x; delete [] ext_y;
  delete [] Nn; delete [] dNn_ds; delete [] dNn_dt;
  delete [] dRr_ds; delete [] dRr_dt;
  delete [] Nns; delete [] Nnt;
  delete [] dNns_ds; delete [] dNnt_dt;
  delete [] Bbs; delete [] Bbt;
  delete [] dBbs_ds; delete [] dBbt_dt;
}


void FEAElement_NURBS_2D_der1::get_R(
    const int &quaindex, double * const &basis ) const
{
  for(int ii=0; ii<nLocBas; ++ii)
    basis[ii] = R[quaindex * nLocBas + ii];
}


void FEAElement_NURBS_2D_der1::get_R_gradR(const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y ) const
{
  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis[ii]   = R[quaindex * nLocBas + ii];
    basis_x[ii] = dR_dx[quaindex * nLocBas + ii];
    basis_y[ii] = dR_dy[quaindex * nLocBas + ii];
  }
}


double FEAElement_NURBS_2D_der1::get_detJac(const int &quaindex) const
{
  return detJac[quaindex];
}


void FEAElement_NURBS_2D_der1::print() const
{
  SYS_T::commPrint("NURBS_2D_der1: ");
  SYS_T::commPrint("2-Dimensional NURBS shape function with 1st order derivatives. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d. \n", get_Type());
  SYS_T::commPrint("Note: Jacobian matrix is not evaluated. \n");
}


double FEAElement_NURBS_2D_der1::get_memory_usage() const
{
  unsigned int double_size = 0; unsigned int int_size = 3;
  double_size = nLocBas * numQuapts * 3 + numQuapts;
  return (double) double_size * 8.0 + (double) int_size * 4.0;
}


void FEAElement_NURBS_2D_der1::BuildShape_atQua(
    const int &sdp1, const int &tdp1,
    const int &quaindex,
    const int &quaindex_s, const int &quaindex_t,
    const double &hx, const double &hy,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z,
    const double * const &ctrl_w,
    double * &N, double * &dN_ds, double * &dN_dt,
    double * &dR_ds, double * &dR_dt,
    const double * const &Ns,
    const double * const &Nt,
    const double * const &dNs_ds,
    const double * const &dNt_dt )
{
  int ii, jj, ll;

  double detJac_temp = 0.0;
  
  double w = 0.0; double dw_ds = 0.0; double dw_dt = 0.0;

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

      dw_ds += dN_ds[ll];
      dw_dt += dN_dt[ll];
    }
  }

  const double inv_w = 1.0 / w;
  int sindex = nLocBas * quaindex;

  double m0 = 0.0; double m1 = 0.0; double m2 = 0.0; double m3 = 0.0;

  for(ii=0; ii<nLocBas; ++ii)
  {
    R[sindex+ii] = N[ii] * inv_w;

    dR_ds[ii] = (dN_ds[ii] - R[sindex+ii] * dw_ds) * inv_w; 
    dR_dt[ii] = (dN_dt[ii] - R[sindex+ii] * dw_dt) * inv_w; 
  
    m0 += ctrl_x[ii] * dR_ds[ii];
    m1 += ctrl_x[ii] * dR_dt[ii];
    m2 += ctrl_y[ii] * dR_ds[ii];
    m3 += ctrl_y[ii] * dR_dt[ii];
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
}

//EOF
