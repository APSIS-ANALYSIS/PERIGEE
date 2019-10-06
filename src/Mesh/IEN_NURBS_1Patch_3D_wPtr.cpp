#include "IEN_NURBS_1Patch_3D_wPtr.hpp"

IEN_NURBS_1Patch_3D_wPtr::IEN_NURBS_1Patch_3D_wPtr( const IMesh * const &mesh )
{
  nElem = mesh->get_nElem();
  nLocBas = mesh->get_nLocBas();

  IEN = new s_int [nElem*nLocBas];

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
                IEN[e*nLocBas+nLocBas - b - 1] = B;     
              }
            }
          }
        }
      }
    }
  }

}


IEN_NURBS_1Patch_3D_wPtr::IEN_NURBS_1Patch_3D_wPtr( const IIEN * const &ien_wz,
           const IMesh * const &mesh )
{
  const int nElem_x_nz = mesh->get_nElem_x_nz();
  const int nElem_y_nz = mesh->get_nElem_y_nz();
  const int nElem_z_nz = mesh->get_nElem_z_nz();

  const int nElem_x = mesh->get_nElem_x();
  const int nElem_y = mesh->get_nElem_y();
  const int nElem_z = mesh->get_nElem_z();

  nElem = mesh->get_nElem_nz();
  nLocBas = mesh->get_nLocBas();

  IEN = new int [nElem * nLocBas];

  int * eindex_x = new int [nElem_x_nz];
  int * eindex_y = new int [nElem_y_nz];
  int * eindex_z = new int [nElem_z_nz];

  int counter = 0;
  for(int ii=0; ii<nElem_x; ++ii)
  {
    if(mesh->get_hx(ii) > 0.0)
    {
      eindex_x[counter] = ii;
      counter += 1;
    }
  }

  counter = 0;
  for(int ii=0; ii<nElem_y; ++ii)
  {
    if(mesh->get_hy(ii*nElem_x) > 0.0)
    {
      eindex_y[counter] = ii;
      counter += 1;
    }
  }
 
  counter = 0;
  for(int ii=0; ii<nElem_z; ++ii)
  {
    if(mesh->get_hz(ii * nElem_x * nElem_y) > 0.0)
    {
      eindex_z[counter] = ii;
      counter += 1;
    }
  } 

  int e = 0, enz = 0;
  for(int kk=0; kk<nElem_z_nz; ++kk)
  {
    for(int jj=0; jj<nElem_y_nz; ++jj)
    {
      for(int ii=0; ii<nElem_x_nz; ++ii)
      {
        e = eindex_x[ii] + eindex_y[jj]*nElem_x + eindex_z[kk] * nElem_x * nElem_y;

        enz = ii + jj * nElem_x_nz + kk * nElem_x_nz * nElem_y_nz;
        
        for(int mm=0; mm<nLocBas; ++mm) IEN[enz*nLocBas + mm] = ien_wz->get_IEN(e, mm);
      }
    }
  }

  delete [] eindex_x; delete [] eindex_y; delete [] eindex_z;
}




IEN_NURBS_1Patch_3D_wPtr::~IEN_NURBS_1Patch_3D_wPtr()
{
  delete [] IEN; IEN = NULL;
}

void IEN_NURBS_1Patch_3D_wPtr::print_IEN() const
{
  std::cout<<std::endl;
  std::cout<<"====== IEN ====== \n";
  for(s_int ii=0; ii<nElem; ++ii)
  {
    for(s_int jj=0; jj<nLocBas; ++jj)
      std::cout<<get_IEN(ii, jj)<<'\t';
    std::cout<<std::endl;
  }
  std::cout<<"================= \n";
}

inline s_int IEN_NURBS_1Patch_3D_wPtr::get_IEN( const s_int &e, const s_int &l_node) const
{
  return IEN[e*nLocBas + l_node];
}

// EOF
