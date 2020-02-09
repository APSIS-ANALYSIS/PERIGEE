#include "AInt_Weight.hpp"

AInt_Weight::AInt_Weight( const IQuadPts * const &qua_s,
    const IQuadPts * const &qua_t,
    const IQuadPts * const &qua_u )
{
  const int n_s = qua_s->get_num_quadPts();
  const int n_t = qua_t->get_num_quadPts();
  const int n_u = qua_u->get_num_quadPts();

  num = n_s * n_t * n_u;

  Weight = new double [num];

  int index;
  double wei_s, wei_t, wei_u;
  for(int ii=0; ii<n_s; ++ii)
  {
    wei_s = qua_s->get_qw(ii);
    for(int jj=0; jj<n_t; ++jj)
    {
      wei_t = qua_t->get_qw(jj);
      for(int kk=0; kk<n_u; ++kk)
      {
        wei_u = qua_u->get_qw(kk);
        index = ii + jj * n_s + kk * n_s * n_t;
        Weight[index] = wei_s * wei_t * wei_u; 
      }
    }
  }
}


AInt_Weight::AInt_Weight( const IQuadPts * const &qua_s,
    const IQuadPts * const &qua_t  )
{
  const int n_s = qua_s->get_num_quadPts();
  const int n_t = qua_t->get_num_quadPts();

  num = n_s * n_t;

  Weight = new double [num];

  int index;
  double wei_s, wei_t;
  for(int ii=0; ii<n_s; ++ii)
  {
    wei_s = qua_s->get_qw(ii);
    for(int jj=0; jj<n_t; ++jj)
    {
      wei_t = qua_t->get_qw(jj);
      index = ii + jj * n_s;
      Weight[index] = wei_s * wei_t; 
    }
  }
}


AInt_Weight::AInt_Weight( const IQuadPts * const &qua )
{
  num = qua->get_num_quadPts();
  Weight = new double [num];
  for(int ii=0; ii<num; ++ii) Weight[ii] = qua -> get_qw(ii);
}


AInt_Weight::~AInt_Weight()
{
  delete [] Weight; Weight = nullptr;
}


void AInt_Weight::print_info() const
{
  std::cout<<'\n';
  std::cout<<"===== AInt_Weight ===== \n";
  std::cout<<"num = "<<num<<std::endl;
  for(int ii=0; ii<num; ++ii) std::cout<<Weight[ii]<<'\n';
  std::cout<<"======================= \n";
}

// EOF
