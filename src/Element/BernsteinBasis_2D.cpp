#include "BernsteinBasis_2D.hpp"

BernsteinBasis_2D::BernsteinBasis_2D( const int &in_sdeg, const int &in_tdeg,
    const IQuadPts * const &in_quaPt_s,
    const IQuadPts * const &in_quaPt_t )
: sdeg(in_sdeg), tdeg(in_tdeg), sdeg_p1(in_sdeg + 1), tdeg_p1(in_tdeg+1),
  sp1tp1(sdeg_p1*tdeg_p1), 
  nqua_s(in_quaPt_s->get_num_quadPts()),
  nqua_t(in_quaPt_t->get_num_quadPts())
{
  bs = new BernsteinBasis_1D(sdeg, in_quaPt_s);
  bt = new BernsteinBasis_1D(tdeg, in_quaPt_t);
}


BernsteinBasis_2D::~BernsteinBasis_2D()
{
  delete bs; delete bt;
  bs = NULL; bt = NULL;
}


double BernsteinBasis_2D::get_B(const int &ii, const int &qua) const
{
  int ii_s, ii_t;
  SYS_T::get_xy_index(ii, sdeg_p1, ii_s, ii_t); // ii = ii_t * (sdeg+1) + ii_s
  int qua_s, qua_t;
  SYS_T::get_xy_index(qua, nqua_s, qua_s, qua_t); // qua = qua_t * nqua_s + qua_s

  return bs->get_B(ii_s, qua_s) * bt->get_B(ii_t, qua_t);
}



double BernsteinBasis_2D::get_dB_ds(const int &ii, const int &qua) const
{
  int ii_s, ii_t;
  SYS_T::get_xy_index(ii, sdeg_p1, ii_s, ii_t); // ii = ii_t * (sdeg+1) + ii_s
  int qua_s, qua_t;
  SYS_T::get_xy_index(qua, nqua_s, qua_s, qua_t); // qua = qua_t * nqua_s + qua_s

  return bs->get_dB_ds(ii_s, qua_s) * bt->get_B(ii_t, qua_t);
}



double BernsteinBasis_2D::get_dB_dt(const int &ii, const int &qua) const
{
  int ii_s, ii_t;
  SYS_T::get_xy_index(ii, sdeg_p1, ii_s, ii_t); // ii = ii_t * (sdeg+1) + ii_s
  int qua_s, qua_t;
  SYS_T::get_xy_index(qua, nqua_s, qua_s, qua_t); // qua = qua_t * nqua_s + qua_s

  return bs->get_B(ii_s, qua_s) * bt->get_dB_ds(ii_t, qua_t);
}



double BernsteinBasis_2D::get_d2B_dss(const int &ii, const int &qua) const
{
  int ii_s, ii_t;
  SYS_T::get_xy_index(ii, sdeg_p1, ii_s, ii_t); // ii = ii_t * (sdeg+1) + ii_s
  int qua_s, qua_t;
  SYS_T::get_xy_index(qua, nqua_s, qua_s, qua_t); // qua = qua_t * nqua_s + qua_s

  return bs->get_d2B_dss(ii_s, qua_s) * bt->get_B(ii_t, qua_t);
}



double BernsteinBasis_2D::get_d2B_dst(const int &ii, const int &qua) const
{
  int ii_s, ii_t;
  SYS_T::get_xy_index(ii, sdeg_p1, ii_s, ii_t); // ii = ii_t * (sdeg+1) + ii_s
  int qua_s, qua_t;
  SYS_T::get_xy_index(qua, nqua_s, qua_s, qua_t); // qua = qua_t * nqua_s + qua_s

  return bs->get_dB_ds(ii_s, qua_s) * bt->get_dB_ds(ii_t, qua_t);
}


double BernsteinBasis_2D::get_d2B_dtt(const int &ii, const int &qua) const
{
  int ii_s, ii_t;
  SYS_T::get_xy_index(ii, sdeg_p1, ii_s, ii_t); // ii = ii_t * (sdeg+1) + ii_s
  int qua_s, qua_t;
  SYS_T::get_xy_index(qua, nqua_s, qua_s, qua_t); // qua = qua_t * nqua_s + qua_s

  return bs->get_B(ii_s, qua_s) * bt->get_d2B_dss(ii_t, qua_t);
}



void BernsteinBasis_2D::get_B(const int &qua_s, const int &qua_t, std::vector<double> &vec) const
{
  vec.resize(sp1tp1);
  std::vector<double> temp_bs, temp_bt;
  
  bs->get_B(qua_s, temp_bs);
  bt->get_B(qua_t, temp_bt);

  int counter = -1;
  for(int jj=0; jj<tdeg_p1; ++jj)
  {
    for(int ii=0; ii<sdeg_p1; ++ii)
    {
      counter += 1;
      vec[counter] = temp_bs[ii] * temp_bt[jj];
    }
  }
  VEC_T::shrink2fit(vec);
}



void BernsteinBasis_2D::get_dB_ds(const int &qua_s, const int &qua_t, std::vector<double> &vec) const
{
  vec.resize(sp1tp1);
  std::vector<double> temp_bs, temp_bt;

  bs->get_dB_ds(qua_s, temp_bs);
  bt->get_B(qua_t, temp_bt);

  int counter = -1;
  for(int jj=0; jj<tdeg_p1; ++jj)
  {
    for(int ii=0; ii<sdeg_p1; ++ii)
    {
      counter += 1;
      vec[counter] = temp_bs[ii] * temp_bt[jj];
    }
  }
  VEC_T::shrink2fit(vec);
}



void BernsteinBasis_2D::get_dB_dt(const int &qua_s, const int &qua_t, std::vector<double> &vec) const
{
  vec.resize(sp1tp1);
  std::vector<double> temp_bs, temp_bt;

  bs->get_B(qua_s, temp_bs);
  bt->get_dB_ds(qua_t, temp_bt);

  int counter = -1;
  for(int jj=0; jj<tdeg_p1; ++jj)
  {
    for(int ii=0; ii<sdeg_p1; ++ii)
    {
      counter += 1;
      vec[counter] = temp_bs[ii] * temp_bt[jj];
    }
  }
  VEC_T::shrink2fit(vec);
}



void BernsteinBasis_2D::get_d2B_dss(const int &qua_s, const int &qua_t, std::vector<double> &vec) const
{
  vec.resize(sp1tp1);
  std::vector<double> temp_bs, temp_bt;

  bs->get_d2B_dss(qua_s, temp_bs);
  bt->get_B(qua_t, temp_bt);

  int counter = -1;
  for(int jj=0; jj<tdeg_p1; ++jj)
  {
    for(int ii=0; ii<sdeg_p1; ++ii)
    {
      counter += 1;
      vec[counter] = temp_bs[ii] * temp_bt[jj];
    }
  }
  VEC_T::shrink2fit(vec);
}



void BernsteinBasis_2D::get_d2B_dst(const int &qua_s, const int &qua_t, std::vector<double> &vec) const
{
  vec.resize(sp1tp1);
  std::vector<double> temp_bs, temp_bt;

  bs->get_dB_ds(qua_s, temp_bs);
  bt->get_dB_ds(qua_t, temp_bt);

  int counter = -1;
  for(int jj=0; jj<tdeg_p1; ++jj)
  {
    for(int ii=0; ii<sdeg_p1; ++ii)
    {
      counter += 1;
      vec[counter] = temp_bs[ii] * temp_bt[jj];
    }
  }
  VEC_T::shrink2fit(vec);
}



void BernsteinBasis_2D::get_d2B_dtt(const int &qua_s, const int &qua_t, std::vector<double> &vec) const
{
  vec.resize(sp1tp1);
  std::vector<double> temp_bs, temp_bt;

  bs->get_B(qua_s, temp_bs);
  bt->get_d2B_dss(qua_t, temp_bt);

  int counter = -1;
  for(int jj=0; jj<tdeg_p1; ++jj)
  {
    for(int ii=0; ii<sdeg_p1; ++ii)
    {
      counter += 1;
      vec[counter] = temp_bs[ii] * temp_bt[jj];
    }
  }
  VEC_T::shrink2fit(vec);
}





void BernsteinBasis_2D::print_info() const
{
  std::cout<<"B: \n";
  for(int ii=0; ii<(sdeg+1)*(tdeg+1); ++ii)
  {
    std::cout<<"basis "<<ii<<'\t';
    for(int jj=0; jj<nqua_s*nqua_t; ++jj)
      std::cout<<get_B(ii,jj)<<'\t';
    std::cout<<std::endl;
  }

  std::cout<<std::endl;

  std::cout<<"dB_ds: \n";
  for(int ii=0; ii<(sdeg+1)*(tdeg+1); ++ii)
  {
    std::cout<<"basis "<<ii<<'\t';
    for(int jj=0; jj<nqua_s*nqua_t; ++jj)
      std::cout<<get_dB_ds(ii,jj)<<'\t';
    std::cout<<std::endl;
  }

  std::cout<<std::endl;

  std::cout<<"dB_dt: \n";
  for(int ii=0; ii<(sdeg+1)*(tdeg+1); ++ii)
  {
    std::cout<<"basis "<<ii<<'\t';
    for(int jj=0; jj<nqua_s*nqua_t; ++jj)
      std::cout<<get_dB_dt(ii,jj)<<'\t';
    std::cout<<std::endl;
  }

  std::cout<<std::endl;
  std::cout<<"d2B_dss: \n";
  for(int ii=0; ii<(sdeg+1)*(tdeg+1); ++ii)
  {
    std::cout<<"basis "<<ii<<'\t';
    for(int jj=0; jj<nqua_s*nqua_t; ++jj)
      std::cout<<get_d2B_dss(ii,jj)<<'\t';
    std::cout<<std::endl;
  }

  std::cout<<std::endl;

  std::cout<<"d2B_dst: \n";
  for(int ii=0; ii<(sdeg+1)*(tdeg+1); ++ii)
  {
    std::cout<<"basis "<<ii<<'\t';
    for(int jj=0; jj<nqua_s*nqua_t; ++jj)
      std::cout<<get_d2B_dst(ii,jj)<<'\t';
    std::cout<<std::endl;
  }
  std::cout<<std::endl;

  std::cout<<"d2B_dtt: \n";
  for(int ii=0; ii<(sdeg+1)*(tdeg+1); ++ii)
  {
    std::cout<<"basis "<<ii<<'\t';
    for(int jj=0; jj<nqua_s*nqua_t; ++jj)
      std::cout<<get_d2B_dtt(ii,jj)<<'\t';
    std::cout<<std::endl;
  }
}

// EOF
