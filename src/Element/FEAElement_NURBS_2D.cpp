#include "FEAElement_NURBS_2D.hpp"

FEAElement_NURBS_2D::FEAElement_NURBS_2D( const int &in_sdeg, const int &in_tdeg,
    const int &in_nquas, const int &in_nquat, const bool &in_is2ndder )
{
  sdeg = in_sdeg;
  tdeg = in_tdeg;
  sdp1 = sdeg + 1;
  tdp1 = tdeg + 1;

  nLocBas = sdp1 * tdp1;

  num_qua_s = in_nquas;
  num_qua_t = in_nquat;

  numQuapts = num_qua_s * num_qua_t;

  rlength = nLocBas * numQuapts;

  is2ndder = in_is2ndder;

  if(is2ndder)
  {
    R   = new double [6*rlength];
    dRr = new double [5*nLocBas];
    Jac = new double [9*numQuapts];
    Nns = new double [3*sdp1];
    Nnt = new double [3*tdp1];
  }
  else
  {
    R   = new double [3*rlength];
    dRr = new double [2*nLocBas];
    Jac = new double [9*numQuapts];
    Nns = new double [2*sdp1];
    Nnt = new double [2*tdp1];
  }
}


FEAElement_NURBS_2D::~FEAElement_NURBS_2D()
{
  clearBasisCache();
}


void FEAElement_NURBS_2D::clearBasisCache()
{
  delete [] R;     R = NULL;
  delete [] dRr; dRr = NULL;
  delete [] Jac; Jac = NULL;
  delete [] Nns; Nns = NULL;
  delete [] Nnt; Nnt = NULL;
}


void FEAElement_NURBS_2D::resize_container()
{
  clearBasisCache();

  if(is2ndder)
  {
    R   = new double [6*rlength];
    dRr = new double [5*nLocBas];
    Jac = new double [9*numQuapts];
    Nns = new double [3*sdp1];
    Nnt = new double [3*tdp1];
  }
  else
  {
    R   = new double [3*rlength];
    dRr = new double [2*nLocBas];
    Jac = new double [9*numQuapts];
    Nns = new double [2*sdp1];
    Nnt = new double [2*tdp1];
  }
}


void FEAElement_NURBS_2D::reset_degree( const int &new_sdeg, const int &new_tdeg )
{
  sdeg = new_sdeg;
  tdeg = new_tdeg;
  sdp1 = sdeg + 1;
  tdp1 = tdeg + 1;
  nLocBas = sdp1 * tdp1;
  rlength = nLocBas * numQuapts;
  resize_container();
}

void FEAElement_NURBS_2D::reset_numQua( const int &new_squa, const int &new_tqua )
{
  num_qua_s = new_squa;
  num_qua_t = new_tqua;
  numQuapts = num_qua_s * num_qua_t;
  rlength = nLocBas * numQuapts;
  resize_container();
}

void FEAElement_NURBS_2D::print() const
{
  SYS_T::commPrint("NURBS_2D : ");
  if(is2ndder)
    SYS_T::commPrint("Two-dimensional NURBS shape function with 2nd derivatives. \n");
  else
    SYS_T::commPrint("Two-dimensional NURBS shape function with 1st derivatives. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d. \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}


double FEAElement_NURBS_2D::get_memory_usage() const
{
  double double_size = numQuapts * 9;
  if(is2ndder)
    double_size += rlength * 6 + nLocBas * 5 + sdp1 * 3 + tdp1 * 3;
  else
    double_size += rlength * 3 + nLocBas * 2 + sdp1 * 2 + tdp1 * 2;

  double int_size = 9;
  return double_size * 8.0 + int_size * 4.0 + 2.0;
}


void FEAElement_NURBS_2D::buildBasis( const double &hx, const double &hy,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_w,
    const double * const &ext_x,
    const double * const &ext_y )
{
  const double invhx = 1.0 / hx;
  const double invhy = 1.0 / hy;

  int ii, jj, ll;
  int counter = -1;

  if(is2ndder)
  {
    for(int qua_y = 0; qua_y < num_qua_t; ++qua_y)
    {
      ll = -1;
      for(ii=0; ii<tdp1; ++ii)
      {
        Nnt[ii] = 0.0; Nnt[tdp1+ii] = 0.0; Nnt[2*tdp1+ii] = 0.0;
        for(jj=0; jj<tdp1; ++jj)
        {
          ll += 1;
          Nnt[ii]        += ext_y[ll] * bt->get_der0(jj,qua_y);
          Nnt[tdp1+ii]   += ext_y[ll] * bt->get_der1(jj,qua_y) * invhy;
          Nnt[2*tdp1+ii] += ext_y[ll] * bt->get_der2(jj,qua_y) * invhy * invhy;
        }
      }

      for(int qua_x = 0; qua_x < num_qua_s; ++qua_x)
      {
        ll = -1;
        for(ii=0; ii<sdp1; ++ii)
        {
          Nns[ii] = 0.0; Nns[sdp1+ii] = 0.0; Nns[2*sdp1+ii] = 0.0;
          for(jj=0; jj<sdp1; ++jj)
          {
            ll += 1;
            Nns[ii]        += ext_x[ll] * bs->get_der0(jj, qua_x);
            Nns[sdp1+ii]   += ext_x[ll] * bs->get_der1(jj, qua_x) * invhx;
            Nns[2*sdp1+ii] += ext_x[ll] * bs->get_der2(jj, qua_x) * invhx * invhx;
          }
        }
        counter += 1;
        BuildShape_atQua2(counter, hx, hy, ctrl_x, ctrl_y, ctrl_w);
      }
    }
  }
  else
  {
    for(int qua_y = 0; qua_y < num_qua_t; ++qua_y)
    {
      ll = -1;
      for(ii=0; ii<tdp1; ++ii)
      {
        Nnt[ii] = 0.0; Nnt[tdp1+ii] = 0.0;
        for(jj=0; jj<tdp1; ++jj)
        {
          ll += 1;
          Nnt[ii]        += ext_y[ll] * bt->get_der0(jj,qua_y);
          Nnt[tdp1+ii]   += ext_y[ll] * bt->get_der1(jj,qua_y) * invhy;
        }
      }

      for(int qua_x = 0; qua_x < num_qua_s; ++qua_x)
      {
        ll = -1;
        for(ii=0; ii<sdp1; ++ii)
        {
          Nns[ii] = 0.0; Nns[sdp1+ii] = 0.0;
          for(jj=0; jj<sdp1; ++jj)
          {
            ll += 1;
            Nns[ii]        += ext_x[ll] * bs->get_der0(jj, qua_x);
            Nns[sdp1+ii]   += ext_x[ll] * bs->get_der1(jj, qua_x) * invhx;
          }
        }
        counter += 1;
        BuildShape_atQua1(counter, hx, hy, ctrl_x, ctrl_y, ctrl_w);
      }
    }
  }
}


void FEAElement_NURBS_2D::BuildShape_atQua1( const int &quaindex, 
    const double &hx, const double &hy, const double * const &ctrl_x,
    const double * const &ctrl_y, const double * const &ctrl_w )
{
  const int offset = nLocBas * quaindex;
  int ii, jj, ll;

  double w = 0.0, dw_ds = 0.0, dw_dt = 0.0;

  ll = -1;

  for(jj=0; jj<tdp1; ++jj)
  {
    for(ii=0; ii<sdp1; ++ii)
    {
      ll += 1;
      R[offset + ll] = Nns[ii] * Nnt[jj] * ctrl_w[ll];
      w += R[offset + ll];

      dRr[ll]           = Nns[sdp1+ii] * Nnt[jj] * ctrl_w[ll];
      dw_ds += dRr[ll];

      dRr[nLocBas + ll] = Nns[ii] * Nnt[tdp1+jj] * ctrl_w[ll];
      dw_dt += dRr[nLocBas + ll];
    }
  }

  const double inv_w = 1.0 / w;
  double dx_ds = 0.0, dx_dt = 0.0, dy_ds = 0.0, dy_dt = 0.0;

  for(ii=0; ii<nLocBas; ++ii)
  {
    R[offset+ii] = R[offset+ii] * inv_w;

    dRr[ii] = (dRr[ii] - R[offset+ii] * dw_ds) * inv_w;
    dRr[nLocBas + ii] = (dRr[nLocBas+ii] - R[offset+ii] * dw_dt) * inv_w;

    dx_ds += ctrl_x[ii] * dRr[ii];
    dx_dt += ctrl_x[ii] * dRr[nLocBas + ii];
    dy_ds += ctrl_y[ii] * dRr[ii];
    dy_dt += ctrl_y[ii] * dRr[nLocBas + ii];
  }

  const double detJac_temp = dx_ds * dy_dt - dx_dt * dy_ds;

  Jac[8*numQuapts + quaindex] = detJac_temp * hx * hy;

  const double inv_detJactmp = 1.0 / detJac_temp;

  const int offset_dxds = 4 * quaindex;

  Jac[offset_dxds + 0] = dx_ds;
  Jac[offset_dxds + 1] = dx_dt;
  Jac[offset_dxds + 2] = dy_ds;
  Jac[offset_dxds + 3] = dy_dt;

  const double ds_dx = dy_dt * inv_detJactmp;
  const double ds_dy = -1.0 * dx_dt * inv_detJactmp;
  const double dt_dx = -1.0 * dy_ds * inv_detJactmp;
  const double dt_dy = dx_ds * inv_detJactmp;

  const int offset_dsdx = offset_dxds + 4 * numQuapts;

  Jac[offset_dsdx + 0] = ds_dx;
  Jac[offset_dsdx + 1] = ds_dy;
  Jac[offset_dsdx + 2] = dt_dx;
  Jac[offset_dsdx + 3] = dt_dy;

  for(ii=0; ii<nLocBas; ++ii)
  {
    R[offset +  rlength +ii] = dRr[ii] * ds_dx + dRr[ii + nLocBas] * dt_dx;
    R[offset + 2*rlength+ii] = dRr[ii] * ds_dy + dRr[ii + nLocBas] * dt_dy;
  }
}


void FEAElement_NURBS_2D::BuildShape_atQua2( const int &quaindex, const double &hx,
    const double &hy, const double * const &ctrl_x,
    const double * const &ctrl_y, const double * const &ctrl_w )
{
  const int offset = nLocBas * quaindex;
  int ii, jj, ll;

  double w = 0.0, dw_ds = 0.0, dw_dt = 0.0; 
  double d2w_dss = 0.0, d2w_dst = 0.0, d2w_dtt = 0.0;

  ll = -1;

  for(jj=0; jj<tdp1; ++jj)
  {
    for(ii=0; ii<sdp1; ++ii)
    {
      ll += 1;
      R[offset + ll] = Nns[ii] * Nnt[jj] * ctrl_w[ll];
      w += R[offset + ll];

      dRr[ll]           = Nns[sdp1+ii] * Nnt[jj] * ctrl_w[ll];
      dw_ds += dRr[ll];

      dRr[nLocBas + ll] = Nns[ii] * Nnt[tdp1+jj] * ctrl_w[ll];
      dw_dt += dRr[nLocBas + ll];

      dRr[2*nLocBas + ll] = Nns[2*sdp1+ii] * Nnt[jj] * ctrl_w[ll];
      d2w_dss += dRr[2*nLocBas + ll];

      dRr[3*nLocBas + ll] = Nns[ii] * Nnt[2*tdp1+jj] * ctrl_w[ll];
      d2w_dtt += dRr[3*nLocBas + ll];

      dRr[4*nLocBas + ll] = Nns[sdp1+ii] * Nnt[tdp1+jj] * ctrl_w[ll];
      d2w_dst += dRr[4*nLocBas + ll];
    }
  }

  const double inv_w = 1.0 / w;
  double dx_ds = 0.0, dx_dt = 0.0, dy_ds = 0.0, dy_dt = 0.0;

  double d2x_dss = 0.0, d2x_dtt = 0.0, d2x_dst = 0.0;
  double d2y_dss = 0.0, d2y_dtt = 0.0, d2y_dst = 0.0;

  for(ii=0; ii<nLocBas; ++ii)
  {
    R[offset+ii] = R[offset+ii] * inv_w;

    dRr[ii] = (dRr[ii] - R[offset+ii] * dw_ds) * inv_w;
    dRr[nLocBas + ii] = (dRr[nLocBas+ii] - R[offset+ii] * dw_dt) * inv_w;

    dRr[2*nLocBas + ii] = (dRr[2*nLocBas+ii] - R[offset+ii] * d2w_dss - 2.0*dw_ds*dRr[ii]) * inv_w;
    dRr[3*nLocBas + ii] = (dRr[3*nLocBas+ii] - R[offset+ii] * d2w_dtt - 2.0*dw_dt*dRr[nLocBas+ii]) * inv_w;
    dRr[4*nLocBas + ii] = (dRr[4*nLocBas+ii] - R[offset+ii] * d2w_dst - dw_dt * dRr[ii] - dw_ds*dRr[nLocBas+ii]) * inv_w;

    dx_ds += ctrl_x[ii] * dRr[ii];
    dx_dt += ctrl_x[ii] * dRr[nLocBas + ii];
    dy_ds += ctrl_y[ii] * dRr[ii];
    dy_dt += ctrl_y[ii] * dRr[nLocBas + ii];

    d2x_dss += ctrl_x[ii] * dRr[2*nLocBas + ii];
    d2x_dtt += ctrl_x[ii] * dRr[3*nLocBas + ii];
    d2x_dst += ctrl_x[ii] * dRr[4*nLocBas + ii];

    d2y_dss += ctrl_y[ii] * dRr[2*nLocBas + ii];
    d2y_dtt += ctrl_y[ii] * dRr[3*nLocBas + ii];
    d2y_dst += ctrl_y[ii] * dRr[4*nLocBas + ii];
  }

  const double detJac_temp = dx_ds * dy_dt - dx_dt * dy_ds;

  Jac[8*numQuapts + quaindex] = detJac_temp * hx * hy;

  const double inv_detJactmp = 1.0 / detJac_temp;

  const int offset_dxds = 4 * quaindex;

  Jac[offset_dxds + 0] = dx_ds;
  Jac[offset_dxds + 1] = dx_dt;
  Jac[offset_dxds + 2] = dy_ds;
  Jac[offset_dxds + 3] = dy_dt;

  const double ds_dx = dy_dt * inv_detJactmp;
  const double ds_dy = -1.0 * dx_dt * inv_detJactmp;
  const double dt_dx = -1.0 * dy_ds * inv_detJactmp;
  const double dt_dy = dx_ds * inv_detJactmp;

  const int offset_dsdx = offset_dxds + 4 * numQuapts;

  Jac[offset_dsdx + 0] = ds_dx;
  Jac[offset_dsdx + 1] = ds_dy;
  Jac[offset_dsdx + 2] = dt_dx;
  Jac[offset_dsdx + 3] = dt_dy;

  for(ii=0; ii<nLocBas; ++ii)
  {
    R[offset +  rlength +ii] = dRr[ii] * ds_dx + dRr[ii + nLocBas] * dt_dx;
    R[offset + 2*rlength+ii] = dRr[ii] * ds_dy + dRr[ii + nLocBas] * dt_dy;
  }

  double Rx, Ry, rhs1, rhs2, rhs3;

  Matrix_double_3by3_Array LHS(dx_ds*dx_ds, dy_ds*dy_ds, 2.0*dx_ds*dy_ds,
      dx_dt*dx_ds, dy_ds*dy_dt, dx_ds*dy_dt + dx_dt*dy_ds,
      dx_dt*dx_dt, dy_dt*dy_dt, 2.0*dx_dt*dy_dt);

  LHS.LU_fac();
  
  for(int ii=0; ii<nLocBas; ++ii)
  {
    Rx = R[offset + rlength + ii];
    Ry = R[offset + 2*rlength + ii];
    rhs1 = dRr[2*nLocBas + ii] - Rx * d2x_dss - Ry * d2y_dss;
    rhs2 = dRr[4*nLocBas + ii] - Rx * d2x_dst - Ry * d2y_dst;
    rhs3 = dRr[3*nLocBas + ii] - Rx * d2x_dtt - Ry * d2y_dtt;
    
    LHS.LU_solve(rhs1, rhs2, rhs3, R[offset+3*rlength+ii],
        R[offset+4*rlength+ii], R[offset+5*rlength+ii]);
  }
}


void FEAElement_NURBS_2D::get_2d_normal_left( const int &quaindex,
    double &nx, double &ny, double &line ) const
{
  const int offset = quaindex * 4;
  const double x_s = Jac[offset];
  const double x_t = Jac[offset+1];
  const double y_s = Jac[offset+2];
  const double y_t = Jac[offset+3];

  line = sqrt(x_s * x_s + y_s * y_s);

  const double tx = x_s / line;
  const double ty = y_s / line;

  nx = - x_t + (x_t*tx + y_t*ty) * tx;
  ny = - y_t + (x_t*tx + y_t*ty) * ty;
  normalize_2d_vector(nx, ny);
}


void FEAElement_NURBS_2D::get_2d_normal_right( const int &quaindex,
    double &nx, double &ny, double &line ) const
{
  const int offset = quaindex * 4;
  const double x_s = Jac[offset];
  const double x_t = Jac[offset+1];
  const double y_s = Jac[offset+2];
  const double y_t = Jac[offset+3];

  line = sqrt(x_s * x_s + y_s * y_s);

  const double tx = x_s / line;
  const double ty = y_s / line;

  nx = x_t - (x_t*tx + y_t*ty) * tx;
  ny = y_t - (x_t*tx + y_t*ty) * ty;
  normalize_2d_vector(nx, ny);
}


void FEAElement_NURBS_2D::get_2d_normal_back( const int &quaindex,
    double &nx, double &ny, double &line ) const
{
  const int offset = quaindex * 4;
  const double x_s = Jac[offset];
  const double x_t = Jac[offset+1];
  const double y_s = Jac[offset+2];
  const double y_t = Jac[offset+3];

  line = sqrt(x_t * x_t + y_t * y_t);

  const double tx = x_t / line;
  const double ty = y_t / line;

  nx = -x_s + (x_s*tx + y_s * ty) * tx;
  ny = -y_s + (x_s*tx + y_s * ty) * ty;
  normalize_2d_vector(nx, ny);
}


void FEAElement_NURBS_2D::get_2d_normal_front( const int &quaindex,
    double &nx, double &ny, double &line ) const
{
  const int offset = quaindex * 4;
  const double x_s = Jac[offset];
  const double x_t = Jac[offset+1];
  const double y_s = Jac[offset+2];
  const double y_t = Jac[offset+3];

  line = sqrt(x_t * x_t + y_t * y_t);

  const double tx = x_t / line;
  const double ty = y_t / line;

  nx = x_s - (x_s*tx + y_s * ty) * tx;
  ny = y_s - (x_s*tx + y_s * ty) * ty;
  normalize_2d_vector(nx, ny);
}


void FEAElement_NURBS_2D::get_R( const int &quaindex, double * const &basis ) const
{
  const int offset = quaindex * nLocBas;
  for(int ii=0; ii<nLocBas; ++ii)
    basis[ii] = R[offset + ii];
}


double FEAElement_NURBS_2D::get_detJac(const int &quaindex) const
{
  return Jac[8*numQuapts + quaindex];
}


void FEAElement_NURBS_2D::get_R_gradR( const int &quaindex, double * const &basis,
    double * const &basis_x, double * const &basis_y ) const
{
  const int offset1 = quaindex * nLocBas;
  const int offset2 = offset1 + rlength;
  const int offset3 = offset2 + rlength;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis[ii]   = R[offset1 + ii];
    basis_x[ii] = R[offset2 + ii];
    basis_y[ii] = R[offset3 + ii];
  }
}


void FEAElement_NURBS_2D::get_2D_R_dR_d2R( const int &quaindex, 
    double * const &basis, double * const &basis_x, double * const &basis_y,
    double * const &basis_xx, double * const &basis_yy, 
    double * const &basis_xy ) const
{
  const int offset1 = quaindex * nLocBas;
  const int offset2 = offset1 + rlength;
  const int offset3 = offset2 + rlength;
  const int offset4 = offset3 + rlength;
  const int offset5 = offset4 + rlength;
  const int offset6 = offset5 + rlength;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis[ii]   = R[offset1 + ii];
    basis_x[ii] = R[offset2 + ii];
    basis_y[ii] = R[offset3 + ii];
    basis_xx[ii] = R[offset4 + ii];
    basis_yy[ii] = R[offset5 + ii];
    basis_xy[ii] = R[offset6 + ii];
  }
}


void FEAElement_NURBS_2D::get_2D_R_gradR_LaplacianR( const int &quaindex,
    double * const &basis, double * const &basis_x, double * const &basis_y,
    double * const &basis_xx, double * const &basis_yy ) const
{
  const int offset1 = quaindex * nLocBas;
  const int offset2 = offset1 + rlength;
  const int offset3 = offset2 + rlength;
  const int offset4 = offset3 + rlength;
  const int offset5 = offset4 + rlength;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis[ii]    = R[offset1 + ii];
    basis_x[ii]  = R[offset2 + ii];
    basis_y[ii]  = R[offset3 + ii];
    basis_xx[ii] = R[offset4 + ii];
    basis_yy[ii] = R[offset5 + ii];
  }
}


void FEAElement_NURBS_2D::get_Jacobian(const int &quaindex, double * const &jac_val) const
{
  const int offset = 4 * quaindex;
  jac_val[0] = Jac[offset];
  jac_val[1] = Jac[offset + 1];
  jac_val[2] = Jac[offset + 2];
  jac_val[3] = Jac[offset + 3];
}


void FEAElement_NURBS_2D::get_invJacobian(const int &quaindex, double * const &jac_val) const
{
  const int offset = 4 * quaindex + 4 * numQuapts;
  jac_val[0] = Jac[offset];
  jac_val[1] = Jac[offset+1];
  jac_val[2] = Jac[offset+2];
  jac_val[3] = Jac[offset+3];
}

// EOF
