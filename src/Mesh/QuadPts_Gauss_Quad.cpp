#include "QuadPts_Gauss_Quad.hpp"

QuadPts_Gauss_Quad::QuadPts_Gauss_Quad( const int &in_num_pts_x,
    const int &in_num_pts_y, 
    const double &x_min, const double &x_max, 
    const double &y_min, const double &y_max )
: num_pts( in_num_pts_x * in_num_pts_y ), num_pts_x( in_num_pts_x ), num_pts_y( in_num_pts_y )
{
  qp.clear(); qw.clear();

  // Use QuadPts_Gauss_1D to generate a rule in 1D with in_num_pts_1d points
  const QuadPts_Gauss_1D qpg1d_x( in_num_pts_x, x_min, x_max );
  const QuadPts_Gauss_1D qpg1d_y( in_num_pts_y, y_min, y_max );

  for(int jj=0; jj<in_num_pts_y; ++jj)
  {
    for(int ii=0; ii<in_num_pts_x; ++ii)
    {
      qp.push_back( qpg1d_x.get_qp(ii) );
      qp.push_back( qpg1d_y.get_qp(jj) );
      qw.push_back( qpg1d_x.get_qw(ii) * qpg1d_y.get_qw(jj) );
    }
  }

  qp.shrink_to_fit();
  qw.shrink_to_fit();
}

QuadPts_Gauss_Quad::QuadPts_Gauss_Quad( const int &in_num_pts_1d, 
    const double &x_min, const double &x_max, 
    const double &y_min, const double &y_max )
: QuadPts_Gauss_Quad(in_num_pts_1d, in_num_pts_1d, x_min, x_max, y_min, y_max)
{}

void QuadPts_Gauss_Quad::print_info() const
{
  SYS_T::commPrint("====== Gauss Points for Quad =======\n");
  SYS_T::commPrint("Number of points = %d\n", num_pts);
  SYS_T::commPrint("qp.size() = %d\n", qp.size());
  SYS_T::commPrint("qw.size() = %d\n", qw.size());
  for(int ii=0; ii<num_pts; ++ii)
    SYS_T::commPrint("  %.15f %.15f %.15f\n", 
        qp[2*ii], qp[2*ii+1], qw[ii]);
  SYS_T::commPrint("====================================\n");
}

// EOF
