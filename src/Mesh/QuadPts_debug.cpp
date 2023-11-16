#include "QuadPts_debug.hpp"

QuadPts_debug::QuadPts_debug( const std::vector<double> &in_qp, 
    const std::vector<double> &in_qw, const int &in_dim )
: num_pts( VEC_T::get_size(in_qp) ), dim( in_dim ), qp( in_qp ), qw( in_qw )
{
  SYS_T::print_fatal_if( qp.size() != dim*qw.size(), 
      "Error: QuadPts_debug, input vector does not have the same length.\n");
}

void QuadPts_debug::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"num_pts = "<<num_pts<<std::endl;
  std::cout<<"dim = "<<dim<<std::endl;
  std::cout<<"====== Debug Points ======"<<std::endl;
  for(int ii=0; ii<num_pts*dim; ++ii)
    std::cout<<qp[ii]<<'\t';
  std::cout<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<qw[ii]<<'\t';
  std::cout<<std::endl;
  std::cout<<"========================="<<std::endl;
}

// EOF
