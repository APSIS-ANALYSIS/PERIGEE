#include "FEAElement_NURBS_3D_der2.hpp"

FEAElement_NURBS_3D_der2::FEAElement_NURBS_3D_der2(
    const int &in_sdeg, const int &in_tdeg,
    const int &in_udeg, const int &in_nquas,
    const int &in_nquat, const int &in_nquau )
{
  sdeg = in_sdeg;
  tdeg = in_tdeg;
  udeg = in_udeg;

  sdp1 = sdeg + 1;
  tdp1 = tdeg + 1;
  udp1 = udeg + 1;

  nLocBas = sdp1 * tdp1 * udp1;

  num_qua_s = in_nquas;
  num_qua_t = in_nquat;
  num_qua_u = in_nquau;

  numQuapts = num_qua_s * num_qua_t * num_qua_u;

  rlength = nLocBas * numQuapts;

  R   = new double [10*rlength];
  Jac = new double [19*numQuapts];
  dRr = new double [9 * nLocBas];
  Nns = new double [3 * sdp1];
  Nnt = new double [3 * tdp1];
  Nnu = new double [3 * udp1];
}


FEAElement_NURBS_3D_der2::~FEAElement_NURBS_3D_der2()
{
  clearBasisCache();
}


void FEAElement_NURBS_3D_der2::clearBasisCache()
{
  delete [] R;     R = NULL;
  delete [] Jac; Jac = NULL;
  delete [] dRr; dRr = NULL;
  delete [] Nns; Nns = NULL;
  delete [] Nnt; Nnt = NULL;
  delete [] Nnu; Nnu = NULL;
}


void FEAElement_NURBS_3D_der2::resize_container()
{
  clearBasisCache();
  R   = new double [10*rlength];
  Jac = new double [19*numQuapts];
  dRr = new double [9 * nLocBas];
  Nns = new double [3 * sdp1];
  Nnt = new double [3 * tdp1];
  Nnu = new double [3 * udp1];
}



void FEAElement_NURBS_3D_der2::reset_degree(
    const int &new_sdeg, const int &new_tdeg,
    const int &new_udeg )
{
  sdeg = new_sdeg;
  tdeg = new_tdeg;
  udeg = new_udeg;
  sdp1 = sdeg + 1;
  tdp1 = tdeg + 1;
  udp1 = udeg + 1;
  nLocBas = sdp1 * tdp1 * udp1;
  rlength = nLocBas * numQuapts;
  resize_container();
}


void FEAElement_NURBS_3D_der2::reset_numQua(
    const int &new_squa, const int &new_tqua,
    const int &new_uqua )
{
  num_qua_s = new_squa;
  num_qua_t = new_tqua;
  num_qua_u = new_uqua;
  numQuapts = num_qua_s * num_qua_t * num_qua_u;
  rlength = nLocBas * numQuapts;
  resize_container();
}


void FEAElement_NURBS_3D_der2::print() const
{
  SYS_T::commPrint("NURBS_3D_der2 : ");
  SYS_T::commPrint("Three-dimensional NURBS shape function with 2nd derivatives. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}


double FEAElement_NURBS_3D_der2::get_memory_usage() const
{
  double double_size = rlength * 10 + numQuapts * 19 + nLocBas * 9
    + 3 * (sdp1+tdp1+udp1);
  double int_size = 11;
  return double_size * 8.0 + int_size * 4.0;
}




void FEAElement_NURBS_3D_der2::buildBasis( const double &hx,
    const double &hy, const double &hz,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const BernsteinBasis_Array * const &bu,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z,
    const double * const &ctrl_w,
    const double * const &ext_x,
    const double * const &ext_y,
    const double * const &ext_z )
{
  const double invhx = 1.0 / hx;
  const double invhy = 1.0 / hy;
  const double invhz = 1.0 / hz;

  int ii, jj, ll;
  int counter = -1;

  for( int qua_z = 0; qua_z < num_qua_u; ++qua_z )
  {
    ll = -1;
    for( ii=0; ii<udp1; ++ii )
    {
      Nnu[ii] = 0.0; Nnu[udp1+ii] = 0.0; Nnu[2*udp1+ii] = 0.0;
      for(jj=0; jj<udp1; ++jj)
      {
        ll += 1;
        Nnu[ii]          += ext_z[ll] * bu->get_der0(jj, qua_z);
        Nnu[udp1 + ii]   += ext_z[ll] * bu->get_der1(jj, qua_z) * invhz;
        Nnu[2*udp1 + ii] += ext_z[ll] * bu->get_der2(jj, qua_z) * invhz * invhz;
      }
    }

    for( int qua_y = 0; qua_y < num_qua_t; ++qua_y )
    {
      ll = -1;
      for( ii=0; ii<tdp1; ++ii )
      {
        Nnt[ii] = 0.0; Nnt[tdp1 + ii] = 0.0; Nnt[2*tdp1 + ii] = 0.0;
        for( jj=0; jj<tdp1; ++jj )
        {
          ll += 1;
          Nnt[ii]          += ext_y[ll] * bt->get_der0(jj, qua_y);
          Nnt[tdp1 + ii]   += ext_y[ll] * bt->get_der1(jj, qua_y) * invhy;
          Nnt[2*tdp1 + ii] += ext_y[ll] * bt->get_der2(jj, qua_y) * invhy * invhy; 
        }
      }

      for( int qua_x = 0; qua_x < num_qua_s; ++qua_x )
      {
        ll = -1;
        for( ii=0; ii<sdp1; ++ii )
        {
          Nns[ii] = 0.0; Nns[sdp1+ii] = 0.0; Nns[2*sdp1+ii] = 0.0;
          for( jj=0; jj<sdp1; ++jj )
          {
            ll += 1;
            Nns[ii]          += ext_x[ll] * bs->get_der0(jj, qua_x);
            Nns[sdp1 + ii]   += ext_x[ll] * bs->get_der1(jj, qua_x) * invhx;
            Nns[2*sdp1 + ii] += ext_x[ll] * bs->get_der2(jj, qua_x) * invhx * invhx;
          }
        }

        counter += 1;
        BuildShape_atQua(counter, hx, hy, hz, ctrl_x, ctrl_y, ctrl_z, ctrl_w);
      }
    }
  }
}


void FEAElement_NURBS_3D_der2::BuildShape_atQua( const int &quaindex,
    const double &hx, const double &hy, const double &hz,
    const double * const &ctrl_x, const double * const &ctrl_y,
    const double * const &ctrl_z, const double * const &ctrl_w )
{
  const int offset = nLocBas * quaindex;

  double dx_ds = 0.0, dx_dt = 0.0, dx_du = 0.0;
  double dy_ds = 0.0, dy_dt = 0.0, dy_du = 0.0;
  double dz_ds = 0.0, dz_dt = 0.0, dz_du = 0.0;
  double d2x_dss = 0.0, d2x_dtt = 0.0, d2x_duu = 0.0, d2x_dst = 0.0, d2x_dsu = 0.0, d2x_dtu = 0.0;
  double d2y_dss = 0.0, d2y_dtt = 0.0, d2y_duu = 0.0, d2y_dst = 0.0, d2y_dsu = 0.0, d2y_dtu = 0.0;
  double d2z_dss = 0.0, d2z_dtt = 0.0, d2z_duu = 0.0, d2z_dst = 0.0, d2z_dsu = 0.0, d2z_dtu = 0.0;
  double w = 0.0, dw_ds = 0.0, dw_dt = 0.0, dw_du = 0.0;
  double d2w_dss = 0.0, d2w_dst = 0.0, d2w_dsu = 0.0, d2w_dtt = 0.0, d2w_dtu = 0.0, d2w_duu = 0.0;

  int ii, jj, kk, ll;
  ll = -1;

  for(kk=0; kk<udp1; ++kk)
  {
    for(jj=0; jj<tdp1; ++jj)
    {
      for(ii=0; ii<sdp1; ++ii)
      {
        ll += 1;
        R[offset + ll] = Nns[ii] * Nnt[jj] * Nnu[kk] * ctrl_w[ll];
        w += R[offset + ll];

        dRr[ll] = Nns[sdp1+ii] * Nnt[jj] * Nnu[kk] * ctrl_w[ll];
        dw_ds += dRr[ll];

        dRr[nLocBas + ll] = Nns[ii] * Nnt[tdp1+jj] * Nnu[kk] * ctrl_w[ll];
        dw_dt += dRr[nLocBas + ll];

        dRr[2*nLocBas + ll] = Nns[ii] * Nnt[jj] * Nnu[udp1+kk] * ctrl_w[ll];
        dw_du += dRr[2*nLocBas + ll];

        dRr[3*nLocBas + ll] = Nns[2*sdp1 + ii] * Nnt[jj] * Nnu[kk] * ctrl_w[ll];
        d2w_dss += dRr[3*nLocBas + ll];
        
        dRr[4*nLocBas + ll] = Nns[ii] * Nnt[2*tdp1 + jj] * Nnu[kk] * ctrl_w[ll];
        d2w_dtt += dRr[4*nLocBas + ll];
        
        dRr[5*nLocBas + ll] = Nns[ii] * Nnt[jj] * Nnu[2*udp1 + kk] * ctrl_w[ll];
        d2w_duu += dRr[5*nLocBas + ll];
        
        dRr[6*nLocBas + ll] = Nns[sdp1+ii] * Nnt[tdp1+jj] * Nnu[kk] * ctrl_w[ll];
        d2w_dst += dRr[6*nLocBas + ll];
        
        dRr[7*nLocBas + ll] = Nns[sdp1+ii] * Nnt[jj] * Nnu[udp1+kk] * ctrl_w[ll];
        d2w_dsu += dRr[7*nLocBas + ll];

        dRr[8*nLocBas + ll] = Nns[ii] * Nnt[tdp1+jj] * Nnu[udp1+kk] * ctrl_w[ll];
        d2w_dtu += dRr[8*nLocBas + ll];
      }
    }
  }

  const double inv_w = 1.0 / w;

  for(ii=0; ii<nLocBas; ++ii)
  {
    R[offset + ii] = R[offset + ii] * inv_w;
    dRr[ii]             = (dRr[ii]             - R[offset+ii] * dw_ds) * inv_w;
    dRr[nLocBas + ii]   = (dRr[nLocBas + ii]   - R[offset+ii] * dw_dt) * inv_w;
    dRr[2*nLocBas + ii] = (dRr[2*nLocBas + ii] - R[offset+ii] * dw_du) * inv_w;
    
    dRr[3*nLocBas + ii] = (dRr[3*nLocBas + ii] - 2.0 * dRr[ii] * dw_ds             - R[offset + ii] * d2w_dss ) * inv_w;
    dRr[4*nLocBas + ii] = (dRr[4*nLocBas + ii] - 2.0 * dRr[nLocBas + ii] * dw_dt   - R[offset + ii] * d2w_dtt ) * inv_w;
    dRr[5*nLocBas + ii] = (dRr[5*nLocBas + ii] - 2.0 * dRr[2*nLocBas + ii] * dw_du - R[offset + ii] * d2w_duu ) * inv_w;
    
    dRr[6*nLocBas + ii] = (dRr[6*nLocBas + ii] - dRr[ii] * dw_dt - dRr[nLocBas + ii] * dw_ds - R[offset + ii] * d2w_dst ) * inv_w;
    dRr[7*nLocBas + ii] = (dRr[7*nLocBas + ii] - dRr[ii] * dw_du - dRr[2*nLocBas + ii] * dw_ds - R[offset + ii] * d2w_dsu ) * inv_w;
    dRr[8*nLocBas + ii] = (dRr[8*nLocBas + ii] - dRr[nLocBas + ii] * dw_du - dRr[2*nLocBas + ii] * dw_dt - R[offset + ii] * d2w_dtu ) * inv_w;

    dx_ds += ctrl_x[ii] * dRr[ii];
    dy_ds += ctrl_y[ii] * dRr[ii];
    dz_ds += ctrl_z[ii] * dRr[ii];

    dx_dt += ctrl_x[ii] * dRr[nLocBas + ii];
    dy_dt += ctrl_y[ii] * dRr[nLocBas + ii];
    dz_dt += ctrl_z[ii] * dRr[nLocBas + ii];

    dx_du += ctrl_x[ii] * dRr[2*nLocBas + ii];
    dy_du += ctrl_y[ii] * dRr[2*nLocBas + ii];
    dz_du += ctrl_z[ii] * dRr[2*nLocBas + ii];

    d2x_dss += ctrl_x[ii] * dRr[3*nLocBas + ii];
    d2x_dtt += ctrl_x[ii] * dRr[4*nLocBas + ii]; 
    d2x_duu += ctrl_x[ii] * dRr[5*nLocBas + ii]; 
    d2x_dst += ctrl_x[ii] * dRr[6*nLocBas + ii]; 
    d2x_dsu += ctrl_x[ii] * dRr[7*nLocBas + ii]; 
    d2x_dtu += ctrl_x[ii] * dRr[8*nLocBas + ii]; 

    d2y_dss += ctrl_y[ii] * dRr[3*nLocBas + ii];
    d2y_dtt += ctrl_y[ii] * dRr[4*nLocBas + ii]; 
    d2y_duu += ctrl_y[ii] * dRr[5*nLocBas + ii]; 
    d2y_dst += ctrl_y[ii] * dRr[6*nLocBas + ii]; 
    d2y_dsu += ctrl_y[ii] * dRr[7*nLocBas + ii]; 
    d2y_dtu += ctrl_y[ii] * dRr[8*nLocBas + ii]; 

    d2z_dss += ctrl_z[ii] * dRr[3*nLocBas + ii];
    d2z_dtt += ctrl_z[ii] * dRr[4*nLocBas + ii]; 
    d2z_duu += ctrl_z[ii] * dRr[5*nLocBas + ii]; 
    d2z_dst += ctrl_z[ii] * dRr[6*nLocBas + ii]; 
    d2z_dsu += ctrl_z[ii] * dRr[7*nLocBas + ii]; 
    d2z_dtu += ctrl_z[ii] * dRr[8*nLocBas + ii]; 
  }

  const double detJac_temp = dx_ds * dy_dt * dz_du + dx_dt * dy_du * dz_ds
    + dx_du * dy_ds * dz_dt - dx_du * dy_dt * dz_ds - dx_dt * dy_ds * dz_du
    - dx_ds * dy_du * dz_dt;

  Jac[18*numQuapts + quaindex] = detJac_temp * hx * hy * hz;

  const int offset_dxds = 9 * quaindex;
  Jac[offset_dxds + 0] = dx_ds;
  Jac[offset_dxds + 1] = dx_dt;
  Jac[offset_dxds + 2] = dx_du;
  Jac[offset_dxds + 3] = dy_ds;
  Jac[offset_dxds + 4] = dy_dt;
  Jac[offset_dxds + 5] = dy_du;
  Jac[offset_dxds + 6] = dz_ds;
  Jac[offset_dxds + 7] = dz_dt;
  Jac[offset_dxds + 8] = dz_du;

  const double inv_detJac = 1.0 / detJac_temp;

  const int offset_dsdx = offset_dxds + 9 * numQuapts;

  double ds_dx = (dy_dt * dz_du - dy_du * dz_dt) * inv_detJac;
  double ds_dy = (dx_du * dz_dt - dx_dt * dz_du) * inv_detJac;
  double ds_dz = (dx_dt * dy_du - dx_du * dy_dt) * inv_detJac;
  double dt_dx = (dy_du * dz_ds - dy_ds * dz_du) * inv_detJac;
  double dt_dy = (dx_ds * dz_du - dx_du * dz_ds) * inv_detJac;
  double dt_dz = (dx_du * dy_ds - dx_ds * dy_du) * inv_detJac;
  double du_dx = (dy_ds * dz_dt - dy_dt * dz_ds) * inv_detJac;
  double du_dy = (dx_dt * dz_ds - dx_ds * dz_dt) * inv_detJac;
  double du_dz = (dx_ds * dy_dt - dx_dt * dy_ds) * inv_detJac;

  Jac[offset_dsdx + 0] = ds_dx;
  Jac[offset_dsdx + 1] = ds_dy;
  Jac[offset_dsdx + 2] = ds_dz;
  Jac[offset_dsdx + 3] = dt_dx;
  Jac[offset_dsdx + 4] = dt_dy;
  Jac[offset_dsdx + 5] = dt_dz;
  Jac[offset_dsdx + 6] = du_dx;
  Jac[offset_dsdx + 7] = du_dy;
  Jac[offset_dsdx + 8] = du_dz;

  for(ii=0; ii<nLocBas; ++ii)
  {
    R[offset + rlength + ii] = dRr[ii] * ds_dx
      + dRr[ii + nLocBas] * dt_dx + dRr[ii + 2*nLocBas] * du_dx;

    R[offset + 2*rlength + ii] = dRr[ii] * ds_dy
      + dRr[ii + nLocBas] * dt_dy + dRr[ii + 2*nLocBas] * du_dy;

    R[offset + 3*rlength + ii] = dRr[ii] * ds_dz
      + dRr[ii + nLocBas] * dt_dz + dRr[ii + 2*nLocBas] * du_dz;
  }

  Matrix_double_6by6_Array LHS(dx_ds, dx_dt, dx_du, dy_ds, dy_dt, dy_du,
      dz_ds, dz_dt, dz_du);

  LHS.LU_fac();
  double RHS [6];
  double sol [6];
  double Rx, Ry, Rz;

  for(ii=0; ii<nLocBas; ++ii)
  {
    Rx = R[offset + rlength + ii];
    Ry = R[offset + 2*rlength + ii];
    Rz = R[offset + 3*rlength + ii];

    RHS[0] = dRr[3*nLocBas + ii] - Rx * d2x_dss - Ry * d2y_dss - Rz * d2z_dss;
    RHS[1] = dRr[6*nLocBas + ii] - Rx * d2x_dst - Ry * d2y_dst - Rz * d2z_dst;
    RHS[2] = dRr[7*nLocBas + ii] - Rx * d2x_dsu - Ry * d2y_dsu - Rz * d2z_dsu;
    RHS[3] = dRr[4*nLocBas + ii] - Rx * d2x_dtt - Ry * d2y_dtt - Rz * d2z_dtt;
    RHS[4] = dRr[8*nLocBas + ii] - Rx * d2x_dtu - Ry * d2y_dtu - Rz * d2z_dtu;
    RHS[5] = dRr[5*nLocBas + ii] - Rx * d2x_duu - Ry * d2y_duu - Rz * d2z_duu;

    LHS.LU_solve(RHS, sol);

    R[offset + 4*rlength + ii] = sol[0];
    R[offset + 5*rlength + ii] = sol[1];
    R[offset + 6*rlength + ii] = sol[2];
    R[offset + 7*rlength + ii] = sol[3];
    R[offset + 8*rlength + ii] = sol[4];
    R[offset + 9*rlength + ii] = sol[5];
  }

}


void FEAElement_NURBS_3D_der2::get_R( const int &quaindex,
    double * const &basis ) const
{
  const int offset = quaindex * nLocBas;
  for(int ii=0; ii<nLocBas; ++ii) basis[ii] = R[offset + ii];
}


double FEAElement_NURBS_3D_der2::get_detJac(const int &quaindex) const
{
  return Jac[18*numQuapts + quaindex];
}



void FEAElement_NURBS_3D_der2::get_R_gradR( const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y, double * const &basis_z ) const
{
  const int offset1 = quaindex * nLocBas;
  const int offset2 = offset1 + rlength;
  const int offset3 = offset2 + rlength;
  const int offset4 = offset3 + rlength;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis[ii]   = R[offset1 + ii];
    basis_x[ii] = R[offset2 + ii];
    basis_y[ii] = R[offset3 + ii];
    basis_z[ii] = R[offset4 + ii];
  }
}


void FEAElement_NURBS_3D_der2::get_3D_R_dR_d2R( const int &quaindex, double * const &basis,
    double * const &basis_x, double * const &basis_y, double * const &basis_z,
    double * const &basis_xx, double * const &basis_yy, double * const &basis_zz,
    double * const &basis_xy, double * const &basis_xz, double * const &basis_yz ) const
{
  const int offset1 = quaindex * nLocBas;
  const int offset2 = offset1 + rlength;
  const int offset3 = offset2 + rlength;
  const int offset4 = offset3 + rlength;
  const int offset5 = offset4 + rlength;
  const int offset6 = offset5 + rlength;
  const int offset7 = offset6 + rlength;
  const int offset8 = offset7 + rlength;
  const int offset9 = offset8 + rlength;
  const int offset10 = offset9 + rlength;

  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis[ii]    = R[offset1 + ii];
    basis_x[ii]  = R[offset2 + ii];
    basis_y[ii]  = R[offset3 + ii];
    basis_z[ii]  = R[offset4 + ii];
    basis_xx[ii] = R[offset5 + ii];
    basis_yy[ii] = R[offset6 + ii];
    basis_zz[ii] = R[offset7 + ii];
    basis_xy[ii] = R[offset8 + ii];
    basis_xz[ii] = R[offset9 + ii];
    basis_yz[ii] = R[offset10 + ii];
  }
}


void FEAElement_NURBS_3D_der2::get_3D_R_gradR_LaplacianR( const int &quaindex,
    double * const &basis, double * const &basis_x, double * const &basis_y,
    double * const &basis_z, double * const &basis_xx,
    double * const &basis_yy, double * const &basis_zz ) const
{
  const int offset1 = quaindex * nLocBas;
  const int offset2 = offset1 + rlength;
  const int offset3 = offset2 + rlength;
  const int offset4 = offset3 + rlength;
  const int offset5 = offset4 + rlength;
  const int offset6 = offset5 + rlength;
  const int offset7 = offset6 + rlength;

  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis[ii]    = R[offset1 + ii];
    basis_x[ii]  = R[offset2 + ii];
    basis_y[ii]  = R[offset3 + ii];
    basis_z[ii]  = R[offset4 + ii];
    basis_xx[ii] = R[offset5 + ii];
    basis_yy[ii] = R[offset6 + ii];
    basis_zz[ii] = R[offset7 + ii];
  }
}

void FEAElement_NURBS_3D_der2::get_Jacobian(const int &quaindex, double * const &jac_value) const
{
  const int offset = 9 * quaindex;
  for(int ii=0; ii<9; ++ii) jac_value[ii] = Jac[offset + ii];
}


void FEAElement_NURBS_3D_der2::get_invJacobian(const int &quaindex, double * const &jac_value) const
{
  const int offset = 9 * quaindex + 9 * numQuapts;
  for(int ii=0; ii<9; ++ii) jac_value[ii] = Jac[offset + ii];
}



void FEAElement_NURBS_3D_der2::get_3d_normal_bottom( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const
{
  // r_t x r_s / |r_t x r_s|
  const int offset = 9 * quaindex;
  nx = Jac[offset+4] * Jac[offset+6] - Jac[offset+3] * Jac[offset+7];
  ny = Jac[offset+7] * Jac[offset+0] - Jac[offset+1] * Jac[offset+6];
  nz = Jac[offset+1] * Jac[offset+3] - Jac[offset+4] * Jac[offset+0];

  normalize_3d_vector(nx, ny, nz, surfaceArea);
}


void FEAElement_NURBS_3D_der2::get_3d_normal_top( const int &quaindex,
    double &nx, double &ny, double &nz, double &surfaceArea ) const
{
  const int offset = 9 * quaindex;
  nx = Jac[offset+3] * Jac[offset+7] - Jac[offset+4] * Jac[offset+6];
  ny = Jac[offset+1] * Jac[offset+6] - Jac[offset+7] * Jac[offset+0];
  nz = Jac[offset+4] * Jac[offset+0] - Jac[offset+1] * Jac[offset+3];

  normalize_3d_vector(nx, ny, nz, surfaceArea);
}



void FEAElement_NURBS_3D_der2::get_3d_normal_left( const int &quaindex,
    double &nx, double &ny, double &nz, double &surfaceArea ) const
{
  const int offset = 9 * quaindex;
  nx = Jac[offset+3] * Jac[offset+8] - Jac[offset+5] * Jac[offset+6];
  ny = Jac[offset+6] * Jac[offset+2] - Jac[offset+0] * Jac[offset+8];
  nz = Jac[offset+0] * Jac[offset+5] - Jac[offset+3] * Jac[offset+2];

  normalize_3d_vector(nx, ny, nz, surfaceArea);
}



void FEAElement_NURBS_3D_der2::get_3d_normal_right( const int &quaindex,
    double &nx, double &ny, double &nz, double &surfaceArea ) const
{
  const int offset = 9 * quaindex;
  nx = Jac[offset+5] * Jac[offset+6] - Jac[offset+3] * Jac[offset+8];
  ny = Jac[offset+0] * Jac[offset+8] - Jac[offset+6] * Jac[offset+2];
  nz = Jac[offset+3] * Jac[offset+2] - Jac[offset+0] * Jac[offset+5];

  normalize_3d_vector(nx, ny, nz, surfaceArea);
}



void FEAElement_NURBS_3D_der2::get_3d_normal_front( const int &quaindex,
    double &nx, double &ny, double &nz, double &surfaceArea ) const
{
  const int offset = 9 * quaindex;
  nx = Jac[offset+4] * Jac[offset+8] - Jac[offset+5] * Jac[offset+7];
  ny = Jac[offset+7] * Jac[offset+2] - Jac[offset+1] * Jac[offset+8];
  nz = Jac[offset+1] * Jac[offset+5] - Jac[offset+4] * Jac[offset+2];

  normalize_3d_vector(nx, ny, nz, surfaceArea);
}



void FEAElement_NURBS_3D_der2::get_3d_normal_back( const int &quaindex,
    double &nx, double &ny, double &nz, double &surfaceArea ) const
{
  const int offset = 9 * quaindex;
  nx = Jac[offset+5] * Jac[offset+7] - Jac[offset+4] * Jac[offset+8];
  ny = Jac[offset+1] * Jac[offset+8] - Jac[offset+7] * Jac[offset+2];
  nz = Jac[offset+4] * Jac[offset+2] - Jac[offset+1] * Jac[offset+5];

  normalize_3d_vector(nx, ny, nz, surfaceArea);
}


// EOF
