#include "QuadPts_Gauss_Hex.hpp"

QuadPts_Gauss_Hex::QuadPts_Gauss_Hex( const int &in_num_pts_x,
    const int &in_num_pts_y, const int &in_num_pts_z,
    const double &x_min, const double &x_max,
    const double &y_min, const double &y_max,
    const double &z_min, const double &z_max )
: num_pts( in_num_pts_x * in_num_pts_y * in_num_pts_z ),
  num_pts_x( in_num_pts_x ), num_pts_y( in_num_pts_y ), num_pts_z( in_num_pts_z )
{
  qp.clear(); qw.clear();

  // Use QuadPts_Gauss_1D to generate a rule in 1D with in_num_pts_1d points
  const QuadPts_Gauss_1D qpg1d_x( in_num_pts_x, x_min, x_max );
  const QuadPts_Gauss_1D qpg1d_y( in_num_pts_y, y_min, y_max );
  const QuadPts_Gauss_1D qpg1d_z( in_num_pts_z, z_min, z_max );
  for(int kk = 0; kk<in_num_pts_z; ++kk)
  {
    for(int jj = 0; jj<in_num_pts_y; ++jj)
    {
      for(int ii = 0; ii<in_num_pts_x; ++ii)
      {
        qp.push_back( qpg1d_x.get_qp(ii) );
        qp.push_back( qpg1d_y.get_qp(jj) );
        qp.push_back( qpg1d_z.get_qp(kk) );
        qw.push_back( qpg1d_x.get_qw(ii) * qpg1d_y.get_qw(jj) * qpg1d_z.get_qw(kk) );
      }
    }
  }
  VEC_T::shrink2fit(qp); VEC_T::shrink2fit(qw);
}

QuadPts_Gauss_Hex::QuadPts_Gauss_Hex( const int &in_num_pts_1d, 
    const double &x_min, const double &x_max, 
    const double &y_min, const double &y_max,
    const double &z_min, const double &z_max )
: QuadPts_Gauss_Hex(in_num_pts_1d, in_num_pts_1d, in_num_pts_1d,
    x_min, x_max, y_min, y_max, z_min, z_max)
{}

QuadPts_Gauss_Hex::~QuadPts_Gauss_Hex()
{
  VEC_T::clean(qp); VEC_T::clean(qw);
}

void QuadPts_Gauss_Hex::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Gauss Points for Hexagon ======="<<std::endl;
  std::cout<<"Num of pt = "<<num_pts<<std::endl;
  std::cout<<"qp.size() = "<<qp.size()<<std::endl;
  std::cout<<"qw.size() = "<<qw.size()<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<std::setprecision(16)
      <<'\t'<<qp[3*ii  ]<<'\t'<<qp[3*ii+1]
      <<'\t'<<qp[3*ii+2]<<'\t'<<qw[ii]<<'\n';
  std::cout<<"==========================================="<<std::endl;
}

// EOF
