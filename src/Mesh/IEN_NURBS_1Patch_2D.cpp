#include "IEN_NURBS_1Patch_2D.hpp"

IEN_NURBS_1Patch_2D::IEN_NURBS_1Patch_2D( const IMesh * const &mesh )
{
  nElem = mesh->get_nElem();
  nLocBas = mesh->get_nLocBas();

  IEN = new int * [nElem];

  for(int ii=0; ii<nElem; ++ii)
    IEN[ii] = new int [nLocBas];

  int nFunc_x = mesh->get_nFunc_x();
  int nFunc_y = mesh->get_nFunc_y();

  int sdegree = mesh->get_s_degree();
  int tdegree = mesh->get_t_degree();

  int A = -1;
  int e = -1;
  int B = 0;
  int b = 0;

  for(int jj=0; jj<nFunc_y; ++jj)
  {
    for(int ii=0; ii<nFunc_x; ++ii)
    {
      A = A + 1;
      if( (ii>=sdegree) && (jj>=tdegree) )
      {
        e = e+1;
        for(int jloc=0; jloc<=tdegree; ++jloc)
        {
          for(int iloc=0; iloc<=sdegree; ++iloc)
          {
            B = A - jloc * nFunc_x - iloc;
            b = jloc * (sdegree+1) + iloc;
            IEN[e][nLocBas - b - 1] = B;
          }
        }
      }
    }
  }

  std::cout<<"\n=== IEN generated. Memory usage: ";
  SYS_T::print_mem_size(nElem * nLocBas * sizeof(int));
  std::cout<<std::endl;
}


IEN_NURBS_1Patch_2D::~IEN_NURBS_1Patch_2D()
{
  for(int ii=0; ii<nElem; ++ii)
    delete [] IEN[ii];
  delete [] IEN;
  std::cout<<"-- IEN array deleted. \n";
}


void IEN_NURBS_1Patch_2D::print_IEN() const
{
  std::cout<<std::endl;
  std::cout<<"====== IEN ====== \n";
  for(int ii=0; ii<nElem; ++ii)
  {
    for(int jj=0; jj<nLocBas; ++jj)
      std::cout<<IEN[ii][jj]<<'\t';
    std::cout<<std::endl;
  }
  std::cout<<"================= \n";
}

// EOF
