#include "IEN_Tetra_P1.hpp"

IEN_Tetra_P1::IEN_Tetra_P1( const int &in_nelem,
    const std::vector<int> &in_ien )
{
  nElem = in_nelem;

  // Linear tets has four local nodes
  IEN = new int [nElem * 4];

  for(unsigned int ii=0; ii<in_ien.size(); ++ii) IEN[ii] = in_ien[ii];
}


IEN_Tetra_P1::~IEN_Tetra_P1()
{
  delete [] IEN; IEN = NULL;
}


int IEN_Tetra_P1::get_IEN( const int &ee, const int &l_node ) const
{
  return IEN[ee*4+l_node];
}


void IEN_Tetra_P1::print_IEN() const
{
  std::cout<<std::endl;
  std::cout<<"====== IEN ====== \n";
  for(int ii=0; ii<nElem; ++ii)
  {
    for(int jj=0; jj<4; ++jj)
      std::cout<<get_IEN(ii, jj)<<'\t';
    std::cout<<std::endl;
  }
  std::cout<<"================= \n";
}

// EOF
