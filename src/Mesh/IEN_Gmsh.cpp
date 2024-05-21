#include "IEN_Gmsh.hpp"

IEN_Gmsh::IEN_Gmsh( const int &in_nelem, const int &in_nlocbas,
    const std::vector<int> &in_ien ) : nElem(in_nelem), nLocBas(in_nlocbas)
{
  IEN = new int [nElem * nLocBas];

  for(unsigned int ii=0; ii<in_ien.size(); ++ii) IEN[ii] = in_ien[ii];
}


IEN_Gmsh::~IEN_Gmsh()
{
  delete [] IEN; IEN = nullptr;
}

int IEN_Gmsh::get_IEN( const int &ee, const int &l_node ) const
{
  return IEN[ee*nLocBas + l_node];
}

void IEN_Gmsh::print_info() const
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
