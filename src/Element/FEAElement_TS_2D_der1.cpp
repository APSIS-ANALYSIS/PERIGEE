#include "FEAElement_TS_2D_der1.hpp"

FEAElement_TS_2D_der1::FEAElement_TS_2D_der1( const int &in_eIndex,
    const IBernsteinBasis * const &bbasis,
    const IALocal_IEN * const &locIEN )
: elem_index(in_eIndex), nLocBas(locIEN->get_nLocBas(elem_index)),
   numQuapts(bbasis->get_nQuaPts_s() * bbasis->get_nQuaPts_t())
{
  R = new double [3*nLocBas * numQuapts + numQuapts];
}

FEAElement_TS_2D_der1::FEAElement_TS_2D_der1( const int &in_eIndex, const int &in_nLocBas,
    const int &in_nQuaPts )
: elem_index(in_eIndex), nLocBas(in_nLocBas), numQuapts(in_nQuaPts)
{
  R = new double [3*nLocBas * numQuapts + numQuapts];
}


FEAElement_TS_2D_der1::~FEAElement_TS_2D_der1()
{
  clearBasisCache();
}

void FEAElement_TS_2D_der1::resize_container()
{
  delete [] R;
  R = new double [3*nLocBas*numQuapts + numQuapts];
}


void FEAElement_TS_2D_der1::buildBasis( const IBernsteinBasis * const &bbasis,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z,
    const double * const &ctrl_w,
    const IAExtractor * const &ext_full )
{
  const int num_qua_s = bbasis->get_nQuaPts_s();
  const int num_qua_t = bbasis->get_nQuaPts_t();
  const int sdeg = ext_full->get_sdegree();
  const int tdeg = ext_full->get_tdegree();
  const int sp1tp1 = (sdeg + 1) * (tdeg + 1);
  const int vseg = nLocBas * numQuapts;

  std::vector<double> vecb, vecbs, vecbt;
  std::vector<double> ext;

  // N, dN_ds, dN_dt in the same vector
  double * N = new double [nLocBas*3];

  // iterator
  int ii, qua_t, qua_s;

  double w , dw_ds, dw_dt, dR_ds, dR_dt;
   
  int counter = -1; 
  for(qua_t=0; qua_t < num_qua_t; ++qua_t)
  {
    for(qua_s=0; qua_s < num_qua_s; ++qua_s)
    {
      counter += 1;

      bbasis->get_B( qua_s, qua_t, vecb );
      bbasis->get_dB_ds( qua_s, qua_t, vecbs );
      bbasis->get_dB_dt( qua_s, qua_t, vecbt ); 
      
      w = 0.0; dw_ds = 0.0; dw_dt = 0.0;

      for(ii=0; ii<nLocBas; ++ii)
      {
        ext_full->get_EXT(elem_index, ii, ext); // ii-th row of the extor mapping

        N[ii]             = VEC_T::vDot(vecb,  ext, sp1tp1) * ctrl_w[ii];
        N[nLocBas + ii]   = VEC_T::vDot(vecbs, ext, sp1tp1) * ctrl_w[ii];
        N[2*nLocBas + ii] = VEC_T::vDot(vecbt, ext, sp1tp1) * ctrl_w[ii]; 
      
        w += N[ii]; 
        dw_ds += N[nLocBas+ii];
        dw_dt += N[2*nLocBas+ii];
      }

      const double inv_w = 1.0 / w;
      const int sindex = nLocBas * counter; 

      double m0 = 0.0; double m1 = 0.0; double m2 = 0.0; double m3 = 0.0;

      for(ii=0; ii<nLocBas; ++ii)
      {
        R[sindex+ii] = N[ii] * inv_w;

        dR_ds = (N[nLocBas+ii] - R[sindex+ii] * dw_ds) * inv_w;
        dR_dt = (N[2*nLocBas+ii] - R[sindex+ii] * dw_dt) * inv_w;

        m0 += ctrl_x[ii] * dR_ds;
        m1 += ctrl_x[ii] * dR_dt;
        m2 += ctrl_y[ii] * dR_ds;
        m3 += ctrl_y[ii] * dR_dt;
      }

      R[3*vseg + counter] = m0 * m3 - m1 * m2;

      const double inv_detJac = 1.0 / R[3*vseg+counter];
      
      const double ds_dx = m3 * inv_detJac;
      const double ds_dy = -1.0 * m1 * inv_detJac;
      const double dt_dx = -1.0 * m2 * inv_detJac;
      const double dt_dy = m0 * inv_detJac; 

      for(ii=0; ii<nLocBas; ++ii)
      {
        dR_ds = (N[nLocBas+ii] - R[sindex+ii] * dw_ds) * inv_w;
        dR_dt = (N[2*nLocBas+ii] - R[sindex+ii] * dw_dt) * inv_w;

        R[vseg + sindex + ii] = dR_ds * ds_dx + dR_dt * dt_dx;
        R[2*vseg+ sindex + ii] = dR_ds * ds_dy + dR_dt * dt_dy; 
      }
    }
  }
  delete [] N;
}


void FEAElement_TS_2D_der1::clearBasisCache()
{
  if(R != NULL)
  {
    delete [] R;
    R = NULL;
  }
}

double FEAElement_TS_2D_der1::get_memory_usage() const
{
  const unsigned int double_size = nLocBas * numQuapts * 3 + numQuapts;
  const unsigned int int_size = 3;
  return double(double_size) * 8.0 + double(int_size) * 4.0;
}

void FEAElement_TS_2D_der1::print() const
{
  SYS_T::commPrint("\n  TS_2D_der1: \n");
  SYS_T::commPrint("2-Dimensional T-Spline basis function with 1st order derivatives. \n");
  SYS_T::commPrint("  -- elemType = 521 \n");
  SYS_T::cmdPrint("  -- elem_index = ", elem_index);
  SYS_T::cmdPrint("  -- nLocBas = ", nLocBas);
  SYS_T::cmdPrint("  -- numQuapts = ", numQuapts);
  SYS_T::commPrint("  -- get_R, get_R_gradR, and get_detJac enabled. \n");
}

void FEAElement_TS_2D_der1::get_R( const int &quaindex, 
    double * const &basis) const
{
  const unsigned int starter = quaindex * nLocBas;
  for(int ii=0; ii<nLocBas; ++ii)
    basis[ii] = R[starter + ii];
}


void FEAElement_TS_2D_der1::get_R_gradR( const int &quaindex,
    double * const &basis, double * const &basis_x,
    double * const &basis_y ) const
{
  const unsigned int starter = quaindex * nLocBas;
  const unsigned int vseg    = numQuapts * nLocBas;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    basis[ii] = R[starter + ii];
    basis_x[ii] = R[vseg + starter + ii];
    basis_y[ii] = R[2*vseg + starter + ii];
  }
}

double FEAElement_TS_2D_der1::get_detJac(const int &quaindex) const
{
  return R[3*nLocBas * numQuapts + quaindex];
}


// EOF
