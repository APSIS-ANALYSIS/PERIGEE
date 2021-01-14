#include "QuadPts_debug.hpp"

QuadPts_debug::QuadPts_debug(const int &len, const std::vector<double> &in_qp,
    const std::vector<double> &in_qw )
: num_pts( len ), dim(1)
{
  qp = in_qp;
  qw = in_qw;

  SYS_T::print_fatal_if( qp.size() != qw.size(), 
      "Error: QuadPts_debug, input vector does not have the same length.\n");
  SYS_T::print_fatal_if( qp.size() != static_cast<unsigned int>(num_pts), 
      "Error: QuadPts_debug, input vector does not match the input parameter.\n");
}


QuadPts_debug::QuadPts_debug( const int &in_dim, const int &in_numpt, 
    const std::vector<double> &in_qp, const std::vector<double> &in_qw )
: num_pts( in_numpt ), dim( in_dim )
{
  qp = in_qp;
  qw = in_qw;

  SYS_T::print_fatal_if( qp.size() != dim*qw.size(), 
      "Error: QuadPts_debug, input vector does not have the same length.\n");
  SYS_T::print_fatal_if( qp.size() != static_cast<unsigned int>(num_pts*dim), 
      "Error: QuadPts_debug, input vector does not match the input parameters.\n");
}


QuadPts_debug::~QuadPts_debug()
{
  VEC_T::clean(qw);
  VEC_T::clean(qp);
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
