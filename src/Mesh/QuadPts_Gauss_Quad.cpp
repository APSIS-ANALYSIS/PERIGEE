#include "QuadPts_Gauss_Quad.hpp"

QuadPts_Gauss_Quad::QuadPts_Gauss_Quad( const int &in_num_pts_1d )
: num_pts( in_num_pts_1d * in_num_pts_1d )
{
  qp.clear(); qw.clear();
  
  // Use QuadPts_Gauss to generate a rule in 1D with in_num_pts_1d points
  const QuadPts_Gauss qpg1d( in_num_pts_1d );

  for(int ii=0; ii<in_num_pts_1d; ++ii)
  {
    for(int jj=0; jj<in_num_pts_1d; ++jj)
    {
      qp.push_back( qpg1d.get_qp(ii) );
      qp.push_back( qpg1d.get_qp(jj) );
      qw.push_back( qpg1d.get_qw(ii) * qpg1d.get_qw(jj) );
    }
  }

  VEC_T::shrink2fit(qp); VEC_T::shrink2fit(qw);
}

QuadPts_Gauss_Quad::~QuadPts_Gauss_Quad()
{
  VEC_T::clean(qp);
  VEC_T::clean(qw);
}

void QuadPts_Gauss_Quad::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Gauss Points for Quad ======="<<std::endl;
  std::cout<<"Num of pt = "<<num_pts<<std::endl;
  std::cout<<"qp.size() = "<<qp.size()<<std::endl;
  std::cout<<"qw.size() = "<<qw.size()<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<std::setprecision(16)<<qp[2*ii]<<'\t'<<qp[2*ii+1]<<'\t'<<qw[ii]<<'\n';
  std::cout<<"========================================="<<std::endl;
}

// EOF
