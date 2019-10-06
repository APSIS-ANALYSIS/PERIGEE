#include "QuadPts_bc0.hpp"

QuadPts_bc0::QuadPts_bc0()
{
  num_pts = 1;

  qp.push_back(0.0);

  qw.push_back(1.0);
}

QuadPts_bc0::~QuadPts_bc0()
{}

void QuadPts_bc0::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== BC 0 Points ======"<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<qp[ii]<<'\t';
  std::cout<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<qw[ii]<<'\t';
  std::cout<<std::endl;
  std::cout<<"========================="<<std::endl;
}

// EOF
