#include "IEN_Tetra_P2.hpp"

IEN_Tetra_P2::IEN_Tetra_P2( const int &in_nelem,
    const std::vector<int> &in_ien ) : nElem( in_nelem )
{
  IEN = new int [nElem * 10];

  for(unsigned int ii=0; ii<in_ien.size(); ++ii) IEN[ii] = in_ien[ii];
}

IEN_Tetra_P2::~IEN_Tetra_P2()
{
  delete [] IEN; IEN = nullptr;
}

int IEN_Tetra_P2::get_IEN( const int &ee, const int &l_node ) const
{
  return IEN[ee*10+l_node];
}

void IEN_Tetra_P2::print_IEN() const
{
  std::cout<<std::endl;
  std::cout<<"====== IEN ====== \n";
  for(int ii=0; ii<nElem; ++ii)
  {
    for(int jj=0; jj<10; ++jj)
      std::cout<<get_IEN(ii, jj)<<'\t';
    std::cout<<std::endl;
  }
  std::cout<<"================= \n";
}

// EOF
