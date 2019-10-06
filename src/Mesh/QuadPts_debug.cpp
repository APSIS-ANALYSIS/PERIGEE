#include "QuadPts_debug.hpp"


QuadPts_debug::QuadPts_debug(const int &len, const std::vector<double> &in_qp,
    const std::vector<double> &in_qw )
: num_pts( len )
{
  qp = in_qp;
  qw = in_qw;
}


QuadPts_debug::~QuadPts_debug()
{
  VEC_T::clean(qw);
  VEC_T::clean(qp);
}


void QuadPts_debug::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Debug Points ======"<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<qp[ii]<<'\t';
  std::cout<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<qw[ii]<<'\t';
  std::cout<<std::endl;
  std::cout<<"========================="<<std::endl;
}

// EOF
