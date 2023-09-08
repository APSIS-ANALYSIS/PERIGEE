#include "QuadPts_Gauss_Quad.hpp"

QuadPts_Gauss_Quad::QuadPts_Gauss_Quad( const int &in_num_pts )
: num_pts( in_num_pts )
{
  qp.resize( 2 * num_pts );
  qw.resize( num_pts );
  std::vector<double> xx, ww;
  
  switch( num_pts )
  {
    case 1:
      xx.push_back(0.5);
      ww.push_back(1.0);
      break;
    case 4:
      xx.push_back(0.788675134594813);
      xx.push_back(0.211324865405187);
      ww.push_back(0.5);
      ww.push_back(0.5);
      break;
    case 9:
      xx.push_back(0.887298334620742);
      xx.push_back(0.5);
      xx.push_back(0.112701665379258);
      ww.push_back(0.277777777777777);
      ww.push_back(0.444444444444444);
      ww.push_back(0.277777777777777);
      break;
    case 16:
      xx.push_back(0.930568155797026);
      xx.push_back(0.669990521792428);
      xx.push_back(0.330009478207572);
      xx.push_back(0.069431844202974);
      ww.push_back(0.173927422568727);
      ww.push_back(0.326072577431273);
      ww.push_back(0.326072577431273);
      ww.push_back(0.173927422568727);
      break;
    case 25:
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
  for(int i = 0; i < nn; ++i){
    for(int j = 0; j < nn; ++j){
      qp[2*(nn*j+i)  ] = xx[i];
      qp[2*(nn*j+i)+1] = xx[j];

      qw[nn*j+i] = ww[i] * ww[j]; 
    }
  }
  VEC_T::shrink2fit(qp); VEC_T::shrink2fit(qw);
}


QuadPts_Gauss_Quad::~QuadPts_Gauss_Quad()
{
  VEC_T::clean(qp);
  VEC_T::clean(qw);
}


void QuadPts_Gauss_Quad::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"====== Gauss Points for Quad ======="<<std::endl;
  std::cout<<"Num of pt = "<<num_pts<<std::endl;
  std::cout<<"qp.size() = "<<qp.size()<<std::endl;
  std::cout<<"qw.size() = "<<qw.size()<<std::endl;
  for(int ii=0; ii<num_pts; ++ii)
    std::cout<<std::setprecision(16)<<qp[2*ii]
      <<'\t'<<qp[2*ii+1]<<'\t'<<qw[ii]<<'\n';
  std::cout<<"========================================="<<std::endl;
}


// EOF
