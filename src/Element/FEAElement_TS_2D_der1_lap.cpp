#include "FEAElement_TS_2D_der1_lap.hpp"

FEAElement_TS_2D_der1_lap::FEAElement_TS_2D_der1_lap(const int &in_eIndex,
    const IBernsteinBasis * const &bbasis,
    const IALocal_IEN * const &locIEN )
: elem_index(in_eIndex), nLocBas(locIEN->get_nLocBas(elem_index)),
  numQuapts(bbasis->get_nQuaPts_s() * bbasis->get_nQuaPts_t())
{
  R = new double [5*nLocBas*numQuapts + numQuapts];
}



FEAElement_TS_2D_der1_lap::FEAElement_TS_2D_der1_lap(const int &in_eIndex,
    const int &in_nLocBas, const int &in_nQuaPts)
: elem_index(in_eIndex), nLocBas(in_nLocBas),
  numQuapts(in_nQuaPts)
{
  R = new double [5*nLocBas*numQuapts + numQuapts];
}



FEAElement_TS_2D_der1_lap::~FEAElement_TS_2D_der1_lap()
{
  clearBasisCache();
}



void FEAElement_TS_2D_der1_lap::clearBasisCache()
{
  if(R!=NULL)
  {
    delete [] R; R = NULL;
  }
}


void FEAElement_TS_2D_der1_lap::print() const
{
  SYS_T::commPrint("\n  TS_2D_der1_lap: \n");
  SYS_T::commPrint("2-Dimensional T-Spline basis function with 1st order derivatives ");
  SYS_T::commPrint("and homogeneous 2nd order derivatives (Laplacian). \n ");
  SYS_T::commPrint("  -- elemType = 523 \n");
  SYS_T::cmdPrint("  -- elem_index = ", elem_index);
  SYS_T::cmdPrint("  -- nLocBas = ", nLocBas);
  SYS_T::cmdPrint("  -- numQuapts = ", numQuapts);
  SYS_T::commPrint("  -- get_R, get_R_gradR_LaplacianR, and get_detJac enabled. \n");
}


double FEAElement_TS_2D_der1_lap::get_memory_usage() const
{
  unsigned int int_size = 3;
  unsigned int double_size = nLocBas * numQuapts * 5 + numQuapts;
  return (double) double_size * 8.0 + (double) int_size * 4.0;
}


void FEAElement_TS_2D_der1_lap::get_R(const int &quaindex, std::vector<double> &basis) const
{
  const unsigned int starter = quaindex * nLocBas;
  basis.resize(nLocBas);
  for(int ii=0; ii<nLocBas; ++ii)
    basis[ii] = R[starter + ii];
}



void FEAElement_TS_2D_der1_lap::get_R_gradR(const int &quaindex,
    std::vector<double> &basis) const
{
  const unsigned int starter = quaindex * nLocBas;
  const unsigned int vseg    = numQuapts * nLocBas;
  basis.resize(3*nLocBas);
  int index_l, index_r;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    index_l = ii; index_r = starter + ii;
    basis[index_l] = R[index_r];
    index_l += nLocBas; index_r += vseg;
    basis[index_l] = R[index_r];
    index_l += nLocBas; index_r += vseg;
    basis[index_l] = R[index_r];
  }
}


void FEAElement_TS_2D_der1_lap::get_2D_R_gradR_LaplacianR(const int &quaindex, 
    std::vector<double> &basis) const
{
  const unsigned int starter = quaindex * nLocBas;
  const unsigned int vseg    = numQuapts * nLocBas;
  basis.resize(5*nLocBas);
  int index_l, index_r;
  for(int ii=0; ii<nLocBas; ++ii)
  {
    index_l = ii; index_r = starter + ii;
    basis[ii] = R[starter + ii];
    index_l += nLocBas; index_r += vseg;
    basis[index_l] = R[index_r];
    index_l += nLocBas; index_r += vseg;
    basis[index_l] = R[index_r];
    index_l += nLocBas; index_r += vseg;
    basis[index_l] = R[index_r];
    index_l += nLocBas; index_r += vseg;
    basis[index_l] = R[index_r];
  }
}


void FEAElement_TS_2D_der1_lap::buildBasis( const IBernsteinBasis * const &bbasis,
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

  std::vector<double> vecb, vecbs, vecbt, vecbss, vecbst, vecbtt;
  std::vector<double> ext;

  // N, dN_ds, dN_dt in the same vector
  double * N  = new double [nLocBas*6];
  double * dR = new double [nLocBas*5];
  double dR_ds, dR_dt, d2R_dss, d2R_dst, d2R_dtt;

  // iterator
  int ii, qua_t, qua_s;

  double w , dw_ds, dw_dt, d2w_dss, d2w_dst, d2w_dtt;
  double temp; 
  double m0, m1, m2, m3, d2x_dss, d2x_dst, d2x_dtt, d2y_dss, d2y_dst, d2y_dtt; 
  double rhs1, rhs2, rhs3;

  int counter = -1; 
  for(qua_t=0; qua_t < num_qua_t; ++qua_t)
  {
    for(qua_s=0; qua_s < num_qua_s; ++qua_s)
    {
      counter += 1;

      bbasis->get_B( qua_s, qua_t, vecb );
      bbasis->get_dB_ds( qua_s, qua_t, vecbs );
      bbasis->get_dB_dt( qua_s, qua_t, vecbt ); 
      bbasis->get_d2B_dss( qua_s, qua_t, vecbss );
      bbasis->get_d2B_dst( qua_s, qua_t, vecbst );
      bbasis->get_d2B_dtt( qua_s, qua_t, vecbtt ); 

      w = 0.0; dw_ds = 0.0; dw_dt = 0.0; 
      d2w_dss = 0.0; d2w_dst = 0.0; d2w_dtt = 0.0;

      for(ii=0; ii<nLocBas; ++ii)
      {
        ext_full->get_EXT(elem_index, ii, ext); // ii-th row of the extor mapping

        N[ii]             = VEC_T::vDot(vecb,  ext, sp1tp1) * ctrl_w[ii];
        N[nLocBas + ii]   = VEC_T::vDot(vecbs, ext, sp1tp1) * ctrl_w[ii];
        N[2*nLocBas + ii] = VEC_T::vDot(vecbt, ext, sp1tp1) * ctrl_w[ii]; 
        N[3*nLocBas + ii] = VEC_T::vDot(vecbss, ext, sp1tp1) * ctrl_w[ii]; 
        N[4*nLocBas + ii] = VEC_T::vDot(vecbst, ext, sp1tp1) * ctrl_w[ii]; 
        N[5*nLocBas + ii] = VEC_T::vDot(vecbtt, ext, sp1tp1) * ctrl_w[ii]; 

        w += N[ii]; 
        dw_ds += N[nLocBas+ii];
        dw_dt += N[2*nLocBas+ii];
        d2w_dss += N[3*nLocBas+ii];
        d2w_dst += N[4*nLocBas+ii];
        d2w_dtt += N[5*nLocBas+ii];
      }

      const double inv_w = 1.0 / w;
      const int sindex = nLocBas * counter; 

      m0 = 0.0; m1 = 0.0; m2 = 0.0; m3 = 0.0;
      d2x_dss = 0.0; d2x_dst = 0.0; d2x_dtt = 0.0;
      d2y_dss = 0.0; d2y_dst = 0.0; d2y_dtt = 0.0;

      for(ii=0; ii<nLocBas; ++ii)
      {
        R[sindex+ii] = N[ii] * inv_w;

        temp = R[sindex+ii];

        // dR_ds = dR[ii]
        dR_ds = (N[nLocBas+ii]   - temp * dw_ds) * inv_w;
        // dR_dt = dR[nLocBas+ii]
        dR_dt = (N[2*nLocBas+ii] - temp * dw_dt) * inv_w;
        // d2R_dss = dR[2*nLocBas+ii]
        d2R_dss = (N[3*nLocBas+ii] - temp * d2w_dss - 2.0*dw_ds*dR_ds) * inv_w;
        // d2R_dst = dR[3*nLocBas+ii]
        d2R_dst = (N[4*nLocBas+ii] - temp * d2w_dst - dw_dt * dR_ds - dw_ds * dR_dt )*inv_w;
        // d2R_dtt = dR[4*nLocBas+ii]
        d2R_dtt = (N[5*nLocBas+ii] - temp * d2w_dtt - 2.0*dw_dt*dR_dt)*inv_w;

        m0 += ctrl_x[ii] * dR_ds;
        m1 += ctrl_x[ii] * dR_dt;
        m2 += ctrl_y[ii] * dR_ds;
        m3 += ctrl_y[ii] * dR_dt;

        d2x_dss += ctrl_x[ii] * d2R_dss;
        d2x_dst += ctrl_x[ii] * d2R_dst;
        d2x_dtt += ctrl_x[ii] * d2R_dtt;
        d2y_dss += ctrl_y[ii] * d2R_dss;
        d2y_dst += ctrl_y[ii] * d2R_dst;
        d2y_dtt += ctrl_y[ii] * d2R_dtt;

        dR[ii] = dR_ds; 
        dR[ii+nLocBas] = dR_dt;
        dR[ii+2*nLocBas] = d2R_dss;
        dR[ii+3*nLocBas] = d2R_dst;
        dR[ii+4*nLocBas] = d2R_dtt;
      }

      R[5*vseg + counter] = m0 * m3 - m1 * m2;

      const double inv_detJac = 1.0 / R[5*vseg+counter];

      const double ds_dx = m3 * inv_detJac;
      const double ds_dy = -1.0 * m1 * inv_detJac;
      const double dt_dx = -1.0 * m2 * inv_detJac;
      const double dt_dy = m0 * inv_detJac; 

      for(ii=0; ii<nLocBas; ++ii)
      {
        dR_ds = dR[ii];
        dR_dt = dR[ii+nLocBas];

        R[vseg + sindex + ii] = dR_ds * ds_dx + dR_dt * dt_dx;
        R[2*vseg+ sindex + ii] = dR_ds * ds_dy + dR_dt * dt_dy; 
      }

      const double inv_bottom = inv_detJac * inv_detJac;
      for(int ii=0; ii<nLocBas; ++ii)
      {
        rhs1 = dR[ii+2*nLocBas] - dR[ii] * d2x_dss - dR[ii+nLocBas] * d2y_dss;
        rhs2 = dR[ii+3*nLocBas] - dR[ii] * d2x_dst - dR[ii+nLocBas] * d2y_dst;
        rhs3 = dR[ii+4*nLocBas] - dR[ii] * d2x_dtt - dR[ii+nLocBas] * d2y_dtt;

        R[3*vseg+sindex+ii] = (m3*m3*rhs1 - 2.0*m2*m3*rhs2 + m2*m2*rhs3) * inv_bottom;
        R[4*vseg+sindex+ii] = (m1*m1*rhs1 - 2.0*m0*m1*rhs2 + m0*m0*rhs3) * inv_bottom;
      }
    }
  }

  delete [] dR;
  delete [] N;
}

// EOF
