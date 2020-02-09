#include "IEN_NURBS_1Patch_3D.hpp"

IEN_NURBS_1Patch_3D::IEN_NURBS_1Patch_3D( const IMesh * const &mesh )
{
  nElem = mesh->get_nElem();
  nLocBas = mesh->get_nLocBas();

  IEN = new int * [nElem];
  for(int ii=0; ii<nElem; ++ii)
    IEN[ii] = new int [nLocBas]; 

  const int nFunc_x = mesh->get_nFunc_x();
  const int nFunc_y = mesh->get_nFunc_y();
  const int nFunc_z = mesh->get_nFunc_z();

  const int sdegree = mesh->get_s_degree();
  const int tdegree = mesh->get_t_degree();
  const int udegree = mesh->get_u_degree();

  int A = -1;
  int e = -1;
  int B = 0;
  int b = 0;
  for(int kk = 0; kk<nFunc_z; ++kk)
  {
    for(int jj=0; jj<nFunc_y; ++jj)
    {
      for(int ii=0; ii<nFunc_x; ++ii)
      {
        A = A + 1; // Global function index increment
        if( (ii>=sdegree) && (jj>=tdegree) && (kk>=udegree) )
        {
          e = e + 1;
          for(int kloc=0; kloc<=udegree; ++kloc)
          {
            for(int jloc=0; jloc<=tdegree; ++jloc)
            {
              for(int iloc=0; iloc<=sdegree; ++iloc)
              {
                B = A - kloc * nFunc_x * nFunc_y - jloc * nFunc_x - iloc;
                b = kloc*(sdegree+1)*(tdegree+1) + jloc * (sdegree+1) + iloc;
                IEN[e][nLocBas - b - 1] = B;     
              }
            }
          }
        }
      }
    }
  }
  std::cout<<"\n=== IEN generated. Memory usage: ";
  SYS_T::print_mem_size(double(nElem) * double(nLocBas) * sizeof(int));
  std::cout<<std::endl;
}

IEN_NURBS_1Patch_3D::~IEN_NURBS_1Patch_3D()
{
  for(int ii=0; ii<nElem; ++ii)
    delete [] IEN[ii];
  delete [] IEN;
  std::cout<<"-- IEN array deleted. \n";
}

void IEN_NURBS_1Patch_3D::print_IEN() const
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
