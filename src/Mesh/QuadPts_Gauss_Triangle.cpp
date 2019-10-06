#include "QuadPts_Gauss_Triangle.hpp"

QuadPts_Gauss_Triangle::QuadPts_Gauss_Triangle( const int &in_num_pts )
: num_pts( in_num_pts )
{
  qp.resize( 3 * num_pts );
  qw.resize( num_pts );

  //int idx = 0;
  double a = 0.0, b = 0.0, c = 0.0, w = 0.0;

  switch( num_pts )
  {
    case 3:
      qp[0] = 0.5;
      qp[1] = 0.5;
      qp[2] = 0.0;
      qw[0] = 1.0 / 3.0;

      qp[3] = 0.5;
      qp[4] = 0.0;
      qp[5] = 0.5;
      qw[1] = 1.0 / 3.0;

      qp[6] = 0.0;
      qp[7] = 0.5;
      qp[8] = 0.5;
      qw[2] = 1.0 / 3.0;
      break;

    case 4:
      qp[0] = 1.0 / 3.0;
      qp[1] = 1.0 / 3.0;
      qp[2] = 1.0 / 3.0;
      qw[0] = -0.5625;

      w = 0.520833333333333;
      qp[3] = 0.6;
      qp[4] = 0.2;
      qp[5] = 0.2;
      qw[1] = w;

      qp[6] = 0.2;
      qp[7] = 0.6;
      qp[8] = 0.2;
      qw[2] = w;

      qp[9] = 0.2;
      qp[10] = 0.2;
      qp[11] = 0.6;
      qw[3] = w;
      break;

    case 6:
      a = 0.816847572980459;
      b = 0.091576213509771;
      w = 0.109951743655322;

      qp[0] = a; qp[1] = b; qp[2] = b; qw[0] = w;
      qp[3] = b; qp[4] = a; qp[5] = b; qw[1] = w;
      qp[6] = b; qp[7] = b; qp[8] = a; qw[2] = w;

      a = 0.108103018168070;
      b = 0.445948490915965;
      w = 0.223381589678011;
      qp[9] = a; qp[10] = b; qp[11] = b; qw[3] = w;
      qp[12] = b; qp[13] = a; qp[14] = b; qw[4] = w;
      qp[15] = b; qp[16] = b; qp[17] = a; qw[5] = w;
      break;

    case 13:
      qp[0] = 1.0 / 3.0;
      qp[1] = 1.0 / 3.0;
      qp[2] = 1.0 / 3.0;
      qw[0] = -0.149570044467670;

      a = 0.479308067841923;
      b = 0.260345966079038;
      w = 0.175615257433204;
      qp[3] = a; qp[4] = b; qp[5] = b; qw[1] = w;
      qp[6] = b; qp[7] = a; qp[8] = b; qw[2] = w;
      qp[9] = b; qp[10] = b; qp[11] = a; qw[3] = w;
      
      a = 0.869739794195568;
      b = 0.065130102902216;
      w = 0.053347235608839;
      qp[12] = a; qp[13] = b; qp[14] = b; qw[4] = w;
      qp[15] = b; qp[16] = a; qp[17] = b; qw[5] = w;
      qp[18] = b; qp[19] = b; qp[20] = a; qw[6] = w;

      a = 0.638444188569809;
      b = 0.312865496004875;
      c = 1.0 - a - b;
      w = 0.077113760890257;
      qp[21] = a; qp[22] = b; qp[23] = c; qw[7] = w;
      qp[24] = a; qp[25] = c; qp[26] = b; qw[8] = w;
      qp[27] = b; qp[28] = a; qp[29] = c; qw[9] = w;
      qp[30] = c; qp[31] = a; qp[32] = b; qw[10] = w;
      qp[33] = b; qp[34] = c; qp[35] = a; qw[11] = w;
      qp[36] = c; qp[37] = b; qp[38] = a; qw[12] = w;

      break;

    default:
      SYS_T::print_fatal("Error: QuadPts_Gauss_Triangle: input number of quadrature points is not implemeneted. \n");
      break;
  }

  VEC_T::shrink2fit(qp);
  VEC_T::shrink2fit(qw);

  // Correct the formula by timeing the triangle area 0.5
  for(int ii=0; ii<num_pts; ++ii) qw[ii] *= 0.5;
}

QuadPts_Gauss_Triangle::~QuadPts_Gauss_Triangle()
{
  VEC_T::clean(qp);
  VEC_T::clean(qw);
}

void QuadPts_Gauss_Triangle::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Gauss Points for Triangles ======="<<std::endl;
  std::cout<<"Num of pt = "<<num_pts<<std::endl;
  std::cout<<"qp.size() = "<<qp.size()<<std::endl;
  std::cout<<"qw.size() = "<<qw.size()<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<std::setprecision(16)<<qp[3*ii]
      <<'\t'<<qp[3*ii+1]<<'\t'<<qp[3*ii+2]<<'\t'<<qw[ii]<<'\n';
  std::cout<<"========================================="<<std::endl;
}

// EOF
