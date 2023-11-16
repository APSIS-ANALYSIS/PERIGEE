#include "QuadPts_bc1.hpp"

void QuadPts_bc1::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== BC 1 Points ======"<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<qp[ii]<<'\t';
  std::cout<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<qw[ii]<<'\t';
  std::cout<<std::endl;
  std::cout<<"========================="<<std::endl;
}

// EOF
