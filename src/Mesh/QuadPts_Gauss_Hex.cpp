#include "QuadPts_Gauss_Hex.hpp"

QuadPts_Gauss_Hex::QuadPts_Gauss_Hex( const int &in_num_pts )
: num_pts( in_num_pts )
{
  qp.resize( 3 * num_pts );
  qw.resize( num_pts );
  std::vector<double> xx, ww;
  
  switch( num_pts )
  {
    case 1:
      xx.push_back(0.5);
      ww.push_back(1.0);
      break;
    case 8:
      xx.push_back(0.788675134594813);
      xx.push_back(0.211324865405187);
      ww.push_back(0.5);
      ww.push_back(0.5);
      break;
    case 27:
      xx.push_back(0.887298334620742);
      xx.push_back(0.5);
      xx.push_back(0.112701665379258);
      ww.push_back(0.277777777777777);
      ww.push_back(0.444444444444444);
      ww.push_back(0.277777777777777);
      break;
    case 64:
      xx.push_back(0.930568155797026);
      xx.push_back(0.669990521792428);
      xx.push_back(0.330009478207572);
      xx.push_back(0.069431844202974);
      ww.push_back(0.173927422568727);
      ww.push_back(0.326072577431273);
      ww.push_back(0.326072577431273);
      ww.push_back(0.173927422568727);
      break;
    case 125:
      xx.push_back(0.953089922969332);
      xx.push_back(0.769234655052841);
      xx.push_back(0.500000000000000);
      xx.push_back(0.230765344947158);
      xx.push_back(0.046910077030668);
      ww.push_back(0.118463442528095);
      ww.push_back(0.239314335249683);
      ww.push_back(0.284444444444444);
      ww.push_back(0.239314335249683);
      ww.push_back(0.118463442528095);
      break;
    default:
      SYS_T::print_fatal("Error: QuadPts_Gauss_Hex: input number of quadrature points is not implemented. \n");
      break;
  }
  int nn = xx.size();
  int nn2 = nn * nn; 
  for(int i = 0; i < nn; ++i){
    for(int j = 0; j < nn; ++j){
      for(int k = 0; k < nn; ++k){
        qp[nn2*k+nn*j+i          ] = xx[i];
        qp[nn2*k+nn*j+i+num_pts  ] = xx[j];
        qp[nn2*k+nn*j+i+num_pts*2] = xx[k];

        qw[nn2*k+nn*j+i] = ww[i] * ww[j] * ww[k]; 
      }
    }
  }
  VEC_T::shrink2fit(qp); VEC_T::shrink2fit(qw);
}


QuadPts_Gauss_Hex::~QuadPts_Gauss_Hex()
{
  VEC_T::clean(qp); VEC_T::clean(qw);
}


void QuadPts_Gauss_Hex::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Gauss Points for Tetrahedron ======="<<std::endl;
  std::cout<<"Num of pt = "<<num_pts<<std::endl;
  std::cout<<"qp.size() = "<<qp.size()<<std::endl;
  std::cout<<"qw.size() = "<<qw.size()<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<std::setprecision(16)<<qp[4*ii]
      <<'\t'<<qp[3*ii  ]<<'\t'<<qp[3*ii+1]
      <<'\t'<<qp[3*ii+2]<<'\t'<<qw[ii]<<'\n';
  std::cout<<"==========================================="<<std::endl;
}


void QuadPts_Gauss_Hex::gen_permutations(const double &a,
    const double &b, const double &c, std::vector<double> &out ) const
{
  out.clear();
  out.push_back(a); out.push_back(b); out.push_back(c); out.push_back(c);
  out.push_back(a); out.push_back(c); out.push_back(b); out.push_back(c);
  out.push_back(a); out.push_back(c); out.push_back(c); out.push_back(b);
  
  out.push_back(b); out.push_back(a); out.push_back(c); out.push_back(c);
  out.push_back(b); out.push_back(c); out.push_back(a); out.push_back(c);
  out.push_back(b); out.push_back(c); out.push_back(c); out.push_back(a);
  
  out.push_back(c); out.push_back(a); out.push_back(b); out.push_back(c);
  out.push_back(c); out.push_back(b); out.push_back(a); out.push_back(c);
  
  out.push_back(c); out.push_back(a); out.push_back(c); out.push_back(b);
  out.push_back(c); out.push_back(b); out.push_back(c); out.push_back(a);
  
  out.push_back(c); out.push_back(c); out.push_back(a); out.push_back(b);
  out.push_back(c); out.push_back(c); out.push_back(b); out.push_back(a);
}


// EOF
