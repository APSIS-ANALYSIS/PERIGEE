#include "QuadPts_Gauss_Tet.hpp"

QuadPts_Gauss_Tet::QuadPts_Gauss_Tet( const int &in_num_pts ) : num_pts( in_num_pts )
{
  qp.resize( 4 * num_pts );
  qw.resize( num_pts );
  
  int offset;
  double a,b,c,w;
  std::vector<double> temp;

  switch( num_pts )
  {
    case 1:
      qp[0] = 0.25; qp[1] = 0.25; qp[2] = 0.25; qp[3] = 0.25;
      qw[0] = 1.0;
      break;

    case 4:
      a = 0.5854101966249685;
      b = (1-a)/3.0;

      qp[0]  = a; qp[1]  = b; qp[2] =  b; qp[3] =  b; qw[0] = 0.25;
      qp[4]  = b; qp[5]  = a; qp[6] =  b; qp[7] =  b; qw[1] = 0.25;
      qp[8]  = b; qp[9]  = b; qp[10] = a; qp[11] = b; qw[2] = 0.25;
      qp[12] = b; qp[13] = b; qp[14] = b; qp[15] = a; qw[3] = 0.25;
      break;

    case 5:
      qp[0] = 0.25; qp[1] = 0.25; qp[2] = 0.25; qp[3] = 0.25; qw[0] = -0.8;
      a = 0.5; b = (1-a) / 3.0;
      qp[4]  = a; qp[5]  = b; qp[6]  = b; qp[7]  = b; qw[1] = 0.45;
      qp[8]  = b; qp[9]  = a; qp[10] = b; qp[11] = b; qw[2] = 0.45;
      qp[12] = b; qp[13] = b; qp[14] = a; qp[15] = b; qw[3] = 0.45;
      qp[16] = b; qp[17] = b; qp[18] = b; qp[19] = a; qw[4] = 0.45;
      break;

    case 17:
      a = 0.1325810999384657;
      b = 0.02454003792903;
      c = (1.0 - a - b) / 2.0;
      temp = gen_permutations(a,b,c);
      for(int ii=0; ii<48; ++ii) qp[ii] = temp[ii];
      
      for(int ii=0; ii<12; ++ii) qw[ii] = 0.04528559236327399;

      a = 0.7316369079576180;
      b = (1.0 - a) / 3.0;
      w = 0.06703858372604275;

      offset = 48;
      qp[offset+0]  = a; qp[offset+1]  = b; qp[offset+2]  = b; qp[offset+3]  = b;
      qw[12] = w;

      qp[offset+4]  = b; qp[offset+5]  = a; qp[offset+6]  = b; qp[offset+7]  = b;
      qw[13] = w;

      qp[offset+8]  = b; qp[offset+9]  = b; qp[offset+10] = a; qp[offset+11] = b;
      qw[14] = w;

      qp[offset+12] = b; qp[offset+13] = b; qp[offset+14] = b; qp[offset+15] = a;
      qw[15] = w;

      qp[offset+16] = 0.25; qp[offset+17] = 0.25; 
      qp[offset+18] = 0.25; qp[offset+19] = 0.25;
      qw[16] = 0.1884185567365411;

      break;

    case 29:
      a = 0.05135188412556341;
      b = 0.4860510285706072;
      c = (1.0 - a - b) * 0.5;
      w = 0.04361493840666568;
      temp = gen_permutations(a,b,c);
      for(int ii=0; ii<48; ++ii) qp[ii] = temp[ii];

      for(int ii=0; ii<12; ++ii) qw[ii] = w;

      a = 0.2967538129690260;
      b = 0.6081079894015281;
      c = (1.0 - a - b) * 0.5;
      w = 0.02581167596199161;
      temp = gen_permutations(a,b,c);
      offset = 48;
      for(int ii=0; ii<48; ++ii) qp[offset+ii] = temp[ii];
      
      for(int ii=0; ii<12; ++ii) qw[12+ii] = w;

      a = 0.8277192480479295;
      b = (1.0 -a) / 3.0;
      w = 0.01911983427899124;

      offset = 96;
      qp[offset+0] = a; qp[offset+1] = b; qp[offset+2] = b; qp[offset+3] = b;
      qw[24] = w;

      qp[offset+4] = b; qp[offset+5] = a; qp[offset+6] = b; qp[offset+7] = b;
      qw[25] = w;

      qp[offset+8] = b; qp[offset+9] = b; qp[offset+10] = a; qp[offset+11] = b;
      qw[26] = w;

      qp[offset+12] = b; qp[offset+13] = b; qp[offset+14] = b; qp[offset+15] = a;
      qw[27] = w;

      qp[offset+16] = 0.25; qp[offset+17] = 0.25; 
      qp[offset+18] = 0.25; qp[offset+19] = 0.25;
      qw[28] = 0.09040129046014750;

      break;

    default:
      SYS_T::print_fatal("Error: QuadPts_Gauss_Tet: input number of quadrature points is not implemented. \n");
      break;
  }

  qp.shrink_to_fit(); qw.shrink_to_fit();

  // Correct the formula by time c = 1.0 / 6.0, the volume of ref. tet.
  for(int ii=0; ii<num_pts; ++ii) qw[ii] /= 6.0;
}

void QuadPts_Gauss_Tet::print_info() const
{
  SYS_T::commPrint("====== Gauss Points for Tetrahedron =======\n");
  SYS_T::commPrint("Number of points = %d\n", num_pts);
  SYS_T::commPrint("qp.size() = %d\n", qp.size());
  SYS_T::commPrint("qw.size() = %d\n", qw.size());
  for(int ii=0; ii<num_pts; ++ii)
    SYS_T::commPrint("  %.15f %.15f %.15f %.15f %.15f\n", 
        qp[4*ii], qp[4*ii+1], qp[4*ii+2], qp[4*ii+3], qw[ii]);
  SYS_T::commPrint("===========================================\n");
}

std::vector<double> QuadPts_Gauss_Tet::gen_permutations(const double &a,
    const double &b, const double &c ) const
{
  return {a,b,c,c,
    a,c,b,c,
    a,c,c,b,
    b,a,c,c,
    b,c,a,c,
    b,c,c,a,
    c,a,b,c,
    c,b,a,c,
    c,a,c,b,
    c,b,c,a,
    c,c,a,b,
    c,c,b,a };
}

// EOF
