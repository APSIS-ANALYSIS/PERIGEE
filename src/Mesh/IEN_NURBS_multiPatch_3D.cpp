#include "IEN_NURBS_multiPatch_3D.hpp"

IEN_NURBS_multiPatch_3D::IEN_NURBS_multiPatch_3D( const IMesh * const &mmesh )
: nElem(mmesh->get_nElem()), nLocBas(mmesh->get_nLocBas())
{
  IEN = new int [nElem * nLocBas];

  const int num_patch = mmesh->get_num_patch();

  IIEN * pien_ptr = NULL;

  for(int ii=0; ii<num_patch; ++ii)
  {
    pien_ptr = new IEN_NURBS_1Patch_3D_wPtr(mmesh->get_patch_mesh(ii));

    const int e_start = mmesh->get_patch_mesh(ii)->get_nElem_start();
    const int ien_start = e_start * nLocBas;
    const int pnelem  = mmesh->get_patch_mesh(ii)->get_nElem();
    const int f_start = mmesh->get_patch_mesh(ii)->get_nFunc_start();

    int counter = 0;
    for(int jj=0; jj<pnelem; ++jj)
    {
      for(int kk=0; kk<nLocBas; ++kk)
      {
        IEN[ien_start + counter] = pien_ptr->get_IEN(jj, kk) + f_start;
        counter += 1; 
      }
    }
    delete pien_ptr;
  }

  std::cout<<"\n=== IEN for multi-patch NURBS strongly match mesh generated. Memory usage: ";
  SYS_T::print_mem_size(double(nElem) * double(nLocBas) * sizeof(int));
  std::cout<<std::endl;
}



IEN_NURBS_multiPatch_3D::~IEN_NURBS_multiPatch_3D()
{
  if( IEN != NULL )
    delete [] IEN;
  IEN = NULL;
  std::cout<<"-- IEN array deleted. \n";
}



void IEN_NURBS_multiPatch_3D::print_IEN() const
{
  std::cout<<std::endl;
  std::cout<<"====== IEN ====== \n";
  for(int ii=0; ii<nElem; ++ii)
  {
    std::cout<<ii<<'\t';
    for(int jj=0; jj<nLocBas; ++jj)
      std::cout<<get_IEN(ii, jj)<<'\t';
    std::cout<<std::endl;
  }
  std::cout<<"================= \n";
}



// EOF
