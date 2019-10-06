#include "IEN_NURBS_1Patch_2D_wPtr.hpp"

IEN_NURBS_1Patch_2D_wPtr::IEN_NURBS_1Patch_2D_wPtr( const IMesh * const &mesh )
{
  nElem = mesh->get_nElem();
  nLocBas = mesh->get_nLocBas();

  IEN = new int [nElem * nLocBas];
  
  const int nFunc_x = mesh->get_nFunc_x();
  const int nFunc_y = mesh->get_nFunc_y();

  const int sdegree = mesh->get_s_degree();
  const int tdegree = mesh->get_t_degree();

  int A = -1, e = -1, B = 0, b = 0;
  
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
            IEN[e*nLocBas + nLocBas - b - 1] = B;
          }
        }
      }
    }
  }

  std::cout<<"\n=== IEN generated. Memory usage : ";
  SYS_T::print_mem_size(nElem * nLocBas * sizeof(int));
  std::cout<<std::endl;
}


IEN_NURBS_1Patch_2D_wPtr::IEN_NURBS_1Patch_2D_wPtr( IIEN const * const &ien_wz,
            IMesh const * const &mesh )
{
  const int nElem_x_nz = mesh->get_nElem_x_nz();
  const int nElem_y_nz = mesh->get_nElem_y_nz();
  const int nElem_x = mesh->get_nElem_x();
  nElem   = mesh->get_nElem_nz();
  nLocBas = mesh->get_nLocBas();
  
  IEN = new int [nElem * nLocBas];

  int * eindex_x = new int [nElem_x_nz];
  int * eindex_y = new int [nElem_y_nz];

  int counter = 0;
  for(int ii=0; ii<mesh->get_nElem_x(); ++ii)
  {
    if(mesh->get_hx(ii) > 0.0)
    {
      eindex_x[counter] = ii;
      counter += 1;
    }
  }

  counter = 0;
  for(int ii=0; ii<mesh->get_nElem_y(); ++ii)
  {
    if(mesh->get_hy(ii*nElem_x) > 0.0)
    {
      eindex_y[counter] = ii;
      counter += 1;
    }
  }

  int e = 0, enz = 0;
  for(int jj=0; jj<nElem_y_nz; ++jj)
  {
    for(int ii=0; ii<nElem_x_nz; ++ii)
    {
      e   = eindex_x[ii] + eindex_y[jj] * nElem_x;
      enz = ii + jj * nElem_x_nz;
      for(int kk=0; kk<nLocBas; ++kk)
        IEN[enz*nLocBas + kk] = ien_wz->get_IEN(e, kk);
    }
  }

  delete [] eindex_x; delete [] eindex_y;
}


IEN_NURBS_1Patch_2D_wPtr::~IEN_NURBS_1Patch_2D_wPtr()
{
  delete [] IEN; IEN = NULL;
}



int IEN_NURBS_1Patch_2D_wPtr::get_IEN( const int &e, const int &node ) const
{
  return IEN[e*nLocBas + node];
}



void IEN_NURBS_1Patch_2D_wPtr::print_IEN() const
{
  std::cout<<std::endl;
  std::cout<<"====== IEN ====== \n";
  for(int ii=0; ii<nElem; ++ii)
  {
    for(int jj=0; jj<nLocBas; ++jj)
    {
      std::cout<<get_IEN(ii,jj)<<'\t';
    }
    std::cout<<'\n';
  }
  std::cout<<"================="<<std::endl;
}


int * IEN_NURBS_1Patch_2D_wPtr::get_IENptr() const
{
  return IEN;
}

// EOF
