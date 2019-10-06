#include "FEAElement_BS_2D_der0.hpp"

FEAElement_BS_2D_der0::FEAElement_BS_2D_der0( 
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

  R   = new double [rlength + numQuapts];
  dRr = new double [2*nLocBas];   // dRr_ds, dRr_dt
  Nns = new double [2*sdp1];      // Nns, dNns_ds
  Nnt = new double [2*tdp1];      // Nnt, dNnt_dt
}


FEAElement_BS_2D_der0::~FEAElement_BS_2D_der0()
{
  clearBasisCache();
}


inline void FEAElement_BS_2D_der0::clearBasisCache()
{
  delete [] R;      R = NULL;
  delete [] dRr;  dRr = NULL;
  delete [] Nns;  Nns = NULL;
  delete [] Nnt;  Nnt = NULL;
}



void FEAElement_BS_2D_der0::buildBasis( const double &hx, const double &hy,
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


inline void FEAElement_BS_2D_der0::get_R(const int &quaindex, double * const &basis) const
{
  const int offset = quaindex * nLocBas;
  for(int ii=0; ii<nLocBas; ++ii)
    basis[ii] = R[offset+ii];
}



inline double FEAElement_BS_2D_der0::get_detJac(const int &quaindex) const
{
  return R[rlength + quaindex];
}



inline void FEAElement_BS_2D_der0::print() const
{
  SYS_T::commPrint("BS_2D_der0 : ");
  SYS_T::commPrint("Two-dimensional B-Spline shape function with no derivatives. \n");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d. \n", get_Type());
  SYS_T::commPrint("Note: Jacobian matrix is NOT evaluated. \n");
}



inline double FEAElement_BS_2D_der0::get_memory_usage() const
{
  unsigned int double_size = rlength + numQuapts;
  unsigned int int_size = 9;
  double_size += 2*nLocBas + 2*sdp1 + 2*tdp1;

  return double(double_size) * 8.0 + double(int_size) * 4.0;
}


void FEAElement_BS_2D_der0::reset_degree(const int &new_sdeg, const int &new_tdeg)
{
  sdeg = new_sdeg;
  tdeg = new_tdeg;
  sdp1 = sdeg + 1;
  tdp1 = tdeg + 1;
  nLocBas = sdp1 * tdp1;
  rlength = nLocBas * numQuapts;
  resize_container();
}


void FEAElement_BS_2D_der0::reset_numQua( const int &new_squa, const int &new_tqua )
{
  num_qua_s = new_squa;
  num_qua_t = new_tqua;
  numQuapts = num_qua_s * num_qua_t;
  rlength = nLocBas * numQuapts;
  resize_container();
}


void FEAElement_BS_2D_der0::resize_container()
{
  clearBasisCache();
  R   = new double [rlength + numQuapts];
  dRr = new double [2*nLocBas];
  Nns = new double [2*sdp1];  
  Nnt = new double [2*tdp1]; 
}



void FEAElement_BS_2D_der0::BuildShape_atQua( const int &quaindex,
    const double &hx, const double &hy,
    const double * const &ctrl_x,
    const double * const &ctrl_y )
{
  int ii, jj, ll;
  ll = -1;

  double m0 = 0.0, m1 = 0.0, m2 = 0.0, m3 = 0.0;
  const int offset = nLocBas * quaindex;
  for(jj=0; jj<tdp1; ++jj)
  {
    for(ii=0; ii<sdp1; ++ii)
    {
      ll += 1;
      R[offset + ll] = Nns[ii] * Nnt[jj];

      dRr[ll]   = Nns[sdp1+ii] * Nnt[jj]; // dN_ds

      m0 += ctrl_x[ll] * dRr[ll]; // dx_ds
      m2 += ctrl_y[ll] * dRr[ll]; // dy_ds

      dRr[nLocBas+ll] = Nns[ii] * Nnt[tdp1+jj]; // dN_dt

      m1 += ctrl_x[ll] * dRr[nLocBas+ll]; // dx_dt
      m3 += ctrl_y[ll] * dRr[nLocBas+ll]; // dy_dt
    }
  }

  // detJac at the quaindex point
  R[rlength+quaindex] = (m0*m3 - m1*m2) * hx * hy;
}


// EOF
