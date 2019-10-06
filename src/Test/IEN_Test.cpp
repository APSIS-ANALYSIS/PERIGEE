#include "IEN_Test.hpp"

void TEST_T::IEN_NURBS_1Patch_3D_Test( const IIEN * const &ien,
   const IMesh * const &mesh )
{
  int sdegree = mesh->get_s_degree();
  int tdegree = mesh->get_t_degree();
  int udegree = mesh->get_u_degree();

  s_int nFunc_x = mesh->get_nFunc_x();
  s_int nFunc_y = mesh->get_nFunc_y();

  s_int nElem_x = mesh->get_nElem_x();
  s_int nElem_y = mesh->get_nElem_y();
  s_int nElem_z = mesh->get_nElem_z();

  s_int nElem = nElem_x * nElem_y * nElem_z;
  int nLocBas = (sdegree + 1) * (tdegree + 1) * (udegree + 1);

  std::vector< std::vector<s_int> > myIEN;
  myIEN.resize(nElem);
  for(int ii =0; ii<nElem; ++ii)
    myIEN[ii].resize(nLocBas);


  s_int e = 0;
  s_int lnode = 0;
  for(s_int ex=0; ex<nElem_x; ++ex)
  {
    for(s_int ey=0; ey<nElem_y; ++ey)
    {
      for(s_int ez=0; ez<nElem_z; ++ez)
      {
        e = ez * nElem_x * nElem_y + ey * nElem_x + ex;
        for(s_int kk=0; kk<=udegree; ++kk)
        {
          for(s_int jj=0; jj<=tdegree; ++jj)
          {
            for(s_int ii=0; ii<=sdegree; ++ii)
            {
              lnode = kk *(sdegree +1) * (tdegree + 1) + jj * (sdegree + 1) + ii;
              myIEN[e][lnode] = (ex + ii) + (ey + jj) * nFunc_x + (ez + kk) * nFunc_x * nFunc_y; 
            }
          }
        }
      }
    }
  }

  for(s_int e = 0; e < nElem; ++e)
  {
    for(s_int ii=0; ii<nLocBas; ++ii)
    {
      if(myIEN[e][ii] != ien->get_IEN(e, ii))
      {
        std::cerr<<e<<'\t'<<ii<<'\t'<<myIEN[e][ii]<<'\n';
        exit(1);
      }
    }
  }

  std::cout<<"IEN Test Passed!..."<<std::endl;

}


void TEST_T::IEN_NURBS_multiPatch_3D_Test( const IMesh * const mmesh,
         const IIEN * const &in_ien )
{
  const int numPat = mmesh->get_num_patch();

  s_int test_a, test_b;
  // loop over all patches
  for(int ii=0; ii<numPat; ++ii)
  {
    // Create a IEN array for the local patch. 
    IIEN * patch_ien = new IEN_NURBS_1Patch_3D(mmesh->get_patch_mesh(ii));

    const int e_start = mmesh->get_patch_mesh(ii)->get_nElem_start();
    const int f_start = mmesh->get_patch_mesh(ii)->get_nFunc_start();
    const int nLocBas = mmesh->get_patch_mesh(ii)->get_nLocBas();
    const int nElem   = mmesh->get_patch_mesh(ii)->get_nElem();

    for(int ee=0; ee<nElem; ++ee)
    {
      for(int kk=0; kk<nLocBas; ++kk)
      {
        test_a = patch_ien->get_IEN(ee, kk);
        test_b = in_ien->get_IEN(ee+e_start, kk);
        if( test_a + f_start != test_b )
        {
          std::cerr<<test_a<<'\t'<<test_b<<'\t'<<f_start<<std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    delete patch_ien;
    std::cout<<"Patch "<<ii<<" IEN passed test. \n";
  }
}



void TEST_T::IEN_NURBS_compare( IIEN const * const &ien1, IIEN const * const &ien2,
   const int &nElem, const int &nLocBas )
{
  int const * const ienptr = ien2->get_IENptr();

  for(int ee=0; ee<nElem; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
    {
      if(ien1->get_IEN(ee, ii) != ien2->get_IEN(ee, ii))
      {
        std::cerr<<"Error : element "<<ee<<" node "<<ii<<'\t'<<"does not match. ";
        std::cerr<<ien1->get_IEN(ee, ii)<<" != "<<ien2->get_IEN(ee,ii)<<std::endl;
        exit(EXIT_FAILURE);
      }
      if(ien1->get_IEN(ee, ii) != *(ienptr+ee*nLocBas+ii) )
      {
        std::cerr<<"Error : element "<<ee<<" node "<<ii<<'\t'<<"does not match. ";
        std::cerr<<ien1->get_IEN(ee, ii)<<" != "<<ien2->get_IEN(ee,ii)<<std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  std::cout<<"IENs are identical. \n";
}


void TEST_T::IEN_NURBS_2D_nz_check( IIEN const * const &ien_wz, IIEN const * const &ien_nz,
          IMesh const * const &mesh )
{
  const int nElem = mesh->get_nElem();
  const int nLocBas = mesh->get_nLocBas();

  int enz = -1;
  for(int ee=0; ee<nElem; ++ee)
  {
    if(mesh->get_hx(ee) * mesh->get_hy(ee) > 0.0)
    {
      enz += 1;
      for(int kk=0; kk<nLocBas; ++kk)
      {
        if(ien_wz->get_IEN(ee, kk) != ien_nz->get_IEN(enz, kk))
        {
          std::cerr<<"Error : element "<<ee<<'\t'<<enz<<'\t'<<kk<<std::endl;
        }
      }
    }
  }

  if(enz != mesh->get_nElem_nz()-1)
  {
    std::cerr<<"Error enz !="<<mesh->get_nElem_nz()-1<<std::endl;
  }

  std::cout<<"Nonzero IEN is compatible with the original IEN. \n";
}


// EOF
