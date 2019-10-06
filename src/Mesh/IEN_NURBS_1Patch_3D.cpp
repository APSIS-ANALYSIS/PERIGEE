#include "IEN_NURBS_1Patch_3D.hpp"

IEN_NURBS_1Patch_3D::IEN_NURBS_1Patch_3D( const IMesh * const &mesh )
{
  nElem = mesh->get_nElem();
  nLocBas = mesh->get_nLocBas();

  IEN = new s_int * [nElem];
  for(s_int ii=0; ii<nElem; ++ii)
    IEN[ii] = new s_int [nLocBas]; 

  const s_int nFunc_x = mesh->get_nFunc_x();
  const s_int nFunc_y = mesh->get_nFunc_y();
  const s_int nFunc_z = mesh->get_nFunc_z();

  const int sdegree = mesh->get_s_degree();
  const int tdegree = mesh->get_t_degree();
  const int udegree = mesh->get_u_degree();

  s_int A = -1;
  s_int e = -1;
  s_int B = 0;
  s_int b = 0;
  for(s_int kk = 0; kk<nFunc_z; ++kk)
  {
    for(s_int jj=0; jj<nFunc_y; ++jj)
    {
      for(s_int ii=0; ii<nFunc_x; ++ii)
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
  SYS_T::print_mem_size(double(nElem) * double(nLocBas) * sizeof(s_int));
  std::cout<<std::endl;
}

IEN_NURBS_1Patch_3D::~IEN_NURBS_1Patch_3D()
{
  for(s_int ii=0; ii<nElem; ++ii)
    delete [] IEN[ii];
  delete [] IEN;
  std::cout<<"-- IEN array deleted. \n";
}

void IEN_NURBS_1Patch_3D::print_IEN() const
{
  std::cout<<std::endl;
  std::cout<<"====== IEN ====== \n";
  for(s_int ii=0; ii<nElem; ++ii)
  {
    for(s_int jj=0; jj<nLocBas; ++jj)
      std::cout<<IEN[ii][jj]<<'\t';
    std::cout<<std::endl;
  }
  std::cout<<"================= \n";
}

// EOF
