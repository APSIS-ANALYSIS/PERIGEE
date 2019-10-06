#include "FEAElement_TS_2D_der0.hpp"


FEAElement_TS_2D_der0::FEAElement_TS_2D_der0( const int &in_eIndex,
    const IBernsteinBasis * const &bbasis,
    const IALocal_IEN * const &locIEN )
: elem_index(in_eIndex), nLocBas(locIEN->get_nLocBas(elem_index)),
 numQuapts(bbasis->get_nQuaPts_s() * bbasis->get_nQuaPts_t()) 
{
  // allocate space for holding R and detJac
  R = new double [nLocBas * numQuapts]; 
  detJac = new double [numQuapts];
}


FEAElement_TS_2D_der0::FEAElement_TS_2D_der0( const int &in_eIndex,
            const int &in_nLocBas, const int &in_nQuaPts )
: elem_index(in_eIndex), nLocBas(in_nLocBas), numQuapts(in_nQuaPts)
{
  R = new double [nLocBas * numQuapts];
  detJac = new double [numQuapts];
}



FEAElement_TS_2D_der0::~FEAElement_TS_2D_der0()
{
  clearBasisCache();
}



void FEAElement_TS_2D_der0::clearBasisCache()
{
  if(R != NULL)
  {
    delete [] R;
    R = NULL;
  } 
  if(detJac != NULL)
  {
    delete [] detJac;
    detJac = NULL;
  }
}



double FEAElement_TS_2D_der0::get_memory_usage() const
{
  unsigned int double_size = nLocBas * numQuapts + numQuapts;
  unsigned int int_size = 3;
  return double(double_size) * 8.0 + double(int_size) * 4.0;
}



void FEAElement_TS_2D_der0::print() const
{
  SYS_T::commPrint("\n  TS_2D_der0: \n");
  SYS_T::commPrint("2-Dimensional T-Spline basis function with no derivatives. \n");
  SYS_T::commPrint("  -- elemType = 520 \n");
  SYS_T::cmdPrint("  -- elem_index = ", elem_index);
  SYS_T::cmdPrint("  -- nLocBas = ", nLocBas);
  SYS_T::cmdPrint("  -- numQuapts = ", numQuapts);
  SYS_T::commPrint("  -- get_R and get_detJac enabled. \n");
}



void FEAElement_TS_2D_der0::get_R(const int &quaindex, 
    double * const &basis) const
{
  const unsigned int starter = quaindex * nLocBas;
  for(int ii=0; ii<nLocBas; ++ii)
    basis[ii] = R[starter + ii];
}


void FEAElement_TS_2D_der0::resize_container()
{
  delete [] R; 
  R = new double [nLocBas * numQuapts];
  delete [] detJac; 
  detJac = new double [numQuapts];
}


void FEAElement_TS_2D_der0::buildBasis(const IBernsteinBasis * const &bbasis,
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

  std::vector<double> vecb, vecbs, vecbt;
  std::vector<double> ext;

  // N, dN_ds, dN_dt in the same vector
  double * N = new double [nLocBas*3];

  // iterator
  int ii, qua_t, qua_s;

  double w , dw_ds, dw_dt;

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

        const double dR_ds = (N[nLocBas+ii] - R[sindex+ii] * dw_ds) * inv_w;
        const double dR_dt = (N[2*nLocBas+ii] - R[sindex+ii] * dw_dt) * inv_w;

        m0 += ctrl_x[ii] * dR_ds;
        m1 += ctrl_x[ii] * dR_dt;
        m2 += ctrl_y[ii] * dR_ds;
        m3 += ctrl_y[ii] * dR_dt;
      }

      detJac[counter] = m0 * m3 - m1 * m2;
    }
  }
  delete [] N;
}


// EOF
