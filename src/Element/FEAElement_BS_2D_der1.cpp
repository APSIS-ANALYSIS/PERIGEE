#include "FEAElement_BS_2D_der1.hpp"

FEAElement_BS_2D_der1::FEAElement_BS_2D_der1( 
    const int &in_sdeg, const int &in_tdeg,
    const int &in_nquas, const int &in_nquat )
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

  R   = new double [3*rlength];
  Jac = new double [9*numQuapts];
  dRr = new double [2*nLocBas]; // dR_ds dR_dt
  Nns = new double [2*sdp1];    // Nns, dNn_ds
  Nnt = new double [2*tdp1];    // Nnt, dNn_dt
}



FEAElement_BS_2D_der1::~FEAElement_BS_2D_der1()
{
  clearBasisCache();
}



void FEAElement_BS_2D_der1::clearBasisCache()
{
  delete [] R;     R = NULL;
  delete [] Jac; Jac = NULL;
  delete [] dRr; dRr = NULL;
  delete [] Nns; Nns = NULL;
  delete [] Nnt; Nnt = NULL;
}


void FEAElement_BS_2D_der1::buildBasis( const double &hx, const double &hy,
    const BernsteinBasis_Array * const &bs,
    const BernsteinBasis_Array * const &bt,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ext_x,
    const double * const &ext_y )
{
  const double invhx = 1.0 / hx;
  const double invhy = 1.0 / hy;

  int ii, jj, ll;
  int counter = -1;

  for(int qua_y=0; qua_y < num_qua_t; ++qua_y)
  {
    ll = -1;
    for(ii=0; ii<tdp1; ++ii)
    {
      Nnt[ii] = 0.0; Nnt[tdp1+ii] = 0.0;
      for(jj=0; jj<tdp1; ++jj)
      {
        ll += 1;
        Nnt[ii]      += ext_y[ll] * bt->get_der0(jj,qua_y);
        Nnt[tdp1+ii] += ext_y[ll] * bt->get_der1(jj,qua_y) * invhy;
      }
    }

    for( int qua_x=0; qua_x < num_qua_s; ++qua_x )
    {
      ll = -1;
      for(ii=0; ii<sdp1; ++ii)
      {
        Nns[ii] = 0.0; Nns[sdp1+ii] = 0.0;
        for(jj=0; jj<sdp1; ++jj)
        {
          ll += 1;
          Nns[ii]      += ext_x[ll] * bs->get_der0(jj, qua_x);
          Nns[sdp1+ii] += ext_x[ll] * bs->get_der1(jj, qua_x) * invhx;
        }
      }
      counter += 1;
      BuildShape_atQua(counter, hx, hy, ctrl_x, ctrl_y);
    }
  } // End loop over all quadrature points
}


void FEAElement_BS_2D_der1::BuildShape_atQua( const int &quaindex,
    const double &hx, const double &hy,
    const double * const &ctrl_x,
    const double * const &ctrl_y )
{
  int ii, jj, ll;
  ll = -1;

  const int offset = nLocBas * quaindex;
  
  double dx_ds = 0.0, dx_dt = 0.0, dy_ds = 0.0, dy_dt = 0.0;

  for(jj=0; jj<tdp1; ++jj)
  {
    for(ii=0; ii<sdp1; ++ii)
    {
      ll += 1;
      R[offset + ll] = Nns[ii] * Nnt[jj];

      dRr[ll]   = Nns[sdp1+ii] * Nnt[jj]; // dN_ds

      dx_ds += ctrl_x[ll] * dRr[ll];      // dx_ds
      dy_ds += ctrl_y[ll] * dRr[ll];      // dy_ds

      dRr[nLocBas+ll] = Nns[ii] * Nnt[tdp1+jj]; // dN_dt

      dx_dt += ctrl_x[ll] * dRr[nLocBas+ll];    // dx_dt
      dy_dt += ctrl_y[ll] * dRr[nLocBas+ll];    // dy_dt
    }
  }

  // detJac at the quaindex point
  Jac[8*numQuapts+quaindex] = (dx_ds * dy_dt - dx_dt * dy_ds) * hx * hy;
 
  const int offset_dxds = 4 * quaindex;
  Jac[offset_dxds + 0] = dx_ds;
  Jac[offset_dxds + 1] = dx_dt;
  Jac[offset_dxds + 2] = dy_ds;
  Jac[offset_dxds + 3] = dy_dt;
  
  const double inv_detJac = 1.0 / (dx_ds * dy_dt - dx_dt * dy_ds); 
   
  const int offset_dsdx = offset_dxds + 4 * numQuapts;

  const double ds_dx = dy_dt * inv_detJac;           // ds_dx
  const double ds_dy = dx_dt * (-1.0) * inv_detJac;  // ds_dy
  const double dt_dx = dy_ds * (-1.0) * inv_detJac;  // dt_dx
  const double dt_dy = dx_ds * inv_detJac;           // dt_dy

  Jac[offset_dsdx + 0] = ds_dx;
  Jac[offset_dsdx + 1] = ds_dy;
  Jac[offset_dsdx + 2] = dt_dx;
  Jac[offset_dsdx + 3] = dt_dy;

  for(int ii=0; ii<nLocBas; ++ii)
  {
    R[offset +  rlength  + ii] = dRr[ii] * ds_dx + dRr[ii+nLocBas] * dt_dx;
    R[offset + 2*rlength + ii] = dRr[ii] * ds_dy + dRr[ii+nLocBas] * dt_dy; 
  }
}


inline void FEAElement_BS_2D_der1::get_R(const int &quaindex, double * const &basis) const
{
  const int offset = quaindex * nLocBas;
  for(int ii=0; ii<nLocBas; ++ii)
    basis[ii] = R[offset + ii];
}



inline double FEAElement_BS_2D_der1::get_detJac(const int &quaindex) const
{
  return Jac[8*numQuapts + quaindex];
}


inline void FEAElement_BS_2D_der1::get_R_gradR(const int &quaindex, double * const &basis,
            double * const &basis_x, double * const &basis_y) const
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


inline void FEAElement_BS_2D_der1::get_Jacobian(const int &quaindex, 
    double * const &jac_val) const
{
  const int offset = 4 * quaindex;
  jac_val[0] = Jac[offset];
  jac_val[1] = Jac[offset + 1];
  jac_val[2] = Jac[offset + 2];
  jac_val[3] = Jac[offset + 3];
}


inline void FEAElement_BS_2D_der1::get_invJacobian(const int &quaindex,
    double * const &jac_val) const
{
  const int offset = 4 * quaindex + 4 * numQuapts;
  jac_val[0] = Jac[offset];
  jac_val[1] = Jac[offset+1];
  jac_val[2] = Jac[offset+2];
  jac_val[3] = Jac[offset+3];
}


void FEAElement_BS_2D_der1::reset_degree(const int &new_sdeg, const int &new_tdeg)
{
  sdeg = new_sdeg;
  tdeg = new_tdeg;
  sdp1 = sdeg + 1;
  tdp1 = tdeg + 1;
  nLocBas = sdp1 * tdp1;
  rlength = nLocBas * numQuapts;
  resize_container();
}


void FEAElement_BS_2D_der1::resize_container()
{
  clearBasisCache();
  R   = new double [3*rlength];
  Jac = new double [9*numQuapts];
  dRr = new double [2*nLocBas];
  Nns = new double [2*sdp1];
  Nnt = new double [2*tdp1];
}


void FEAElement_BS_2D_der1::reset_numQua( const int &new_squa, const int &new_tqua )
{
  num_qua_s = new_squa;
  num_qua_t = new_tqua;
  numQuapts = num_qua_s * num_qua_t;
  rlength = nLocBas * numQuapts;
  resize_container();
}



inline void FEAElement_BS_2D_der1::print() const
{
  SYS_T::commPrint("BS_2D_der1 : ");
  SYS_T::commPrint("Two-dimensional B-Spline shape function with 1st derivatives. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d. \n", get_Type());
  SYS_T::commPrint("Note: Jacobian and inverse Jacobian are evaluated. \n");
}


inline double FEAElement_BS_2D_der1::get_memory_usage() const
{
  double double_size = rlength * 3 + numQuapts * 9 + 2*nLocBas + 2*sdp1 + 2*tdp1;
  double int_size = 9;
  return double_size * 8.0 + int_size * 4.0;
}


void FEAElement_BS_2D_der1::get_2d_normal_left( const int &quaindex,
    double &nx, double &ny, double &line ) const
{
  // n = -(r_t - (r_t dot r_s) r_s) / |its norm|
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


void FEAElement_BS_2D_der1::get_2d_normal_right( const int &quaindex,
    double &nx, double &ny, double &line ) const
{
  // n = r_t - (r_t dot r_s) r_s / |its norm|
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


void FEAElement_BS_2D_der1::get_2d_normal_back( const int &quaindex,
    double &nx, double &ny, double &line ) const
{
  // n = - r_s + (r_t dot r_s) r_t / |its norm|
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


void FEAElement_BS_2D_der1::get_2d_normal_front( const int &quaindex,
    double &nx, double &ny, double &line ) const
{
  // n =  r_s - (r_t dot r_s) r_t / |its norm|
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


// EOF
