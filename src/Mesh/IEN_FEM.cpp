#include "IEN_FEM.hpp"

IEN_FEM::IEN_FEM( const int &in_nelem, const std::vector<int> &in_ien ) 
: nElem( in_nelem ), nLocBas( VEC_T::get_size( in_ien )/nElem )
{
  SYS_T::print_exit_if_not( VEC_T::get_size( in_ien ) % nElem == 0,
     "Error: IEN_FEM input IEN vector size cannot do a perfect division by input number of elements. This suggests the mesh main contain different types of elements. \n" );

  IEN = new int [nElem * nLocBas];

  for(unsigned int ii=0; ii<in_ien.size(); ++ii) IEN[ii] = in_ien[ii];
}

IEN_FEM::~IEN_FEM()
{
  delete [] IEN; IEN = nullptr;
}

int IEN_FEM::get_IEN( const int &ee, const int &ii ) const
{
  return IEN[ee*nLocBas + ii];
}

void IEN_FEM::print_IEN() const
{
  std::cout<<std::endl;
  std::cout<<"====== IEN ====== \n";
  for(int ii=0; ii<nElem; ++ii)
  {
    for(int jj=0; jj<nLocBas; ++jj)
      std::cout<<get_IEN(ii, jj)<<'\t';
    std::cout<<std::endl;
  }
  std::cout<<"================= \n";
}

// EOF
