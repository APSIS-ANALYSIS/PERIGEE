#include "QuadPts_Gauss_Hex.hpp"

QuadPts_Gauss_Hex::QuadPts_Gauss_Hex( int in_num_pts_x,
    int in_num_pts_y, int in_num_pts_z,
    double x_min, double x_max,
    double y_min, double y_max,
    double z_min, double z_max )
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
  qp.shrink_to_fit();
  qw.shrink_to_fit();
}

QuadPts_Gauss_Hex::QuadPts_Gauss_Hex( int in_num_pts_1d,
    double x_min, double x_max,
    double y_min, double y_max,
    double z_min, double z_max )
: QuadPts_Gauss_Hex(in_num_pts_1d, in_num_pts_1d, in_num_pts_1d,
    x_min, x_max, y_min, y_max, z_min, z_max)
{}

void QuadPts_Gauss_Hex::print_info() const
{
  SYS_T::commPrint("====== Gauss Points for Hexagon =======\n");
  SYS_T::commPrint("Number of points = %d\n", num_pts);
  SYS_T::commPrint("qp.size() = %d\n", qp.size());
  SYS_T::commPrint("qw.size() = %d\n", qw.size());
  for(int ii=0; ii<num_pts; ++ii)
    SYS_T::commPrint("  %.15f %.15f %.15f %.15f\n",
        qp[3*ii], qp[3*ii+1], qp[3*ii+2], qw[ii]);
  SYS_T::commPrint("=======================================\n");
}

// EOF
