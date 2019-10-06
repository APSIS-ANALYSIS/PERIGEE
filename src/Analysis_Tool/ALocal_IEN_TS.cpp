#include "ALocal_IEN_TS.hpp"

ALocal_IEN_TS::ALocal_IEN_TS(const HDF5_PartReader * const &h5reader)
{
  h5reader->get_1D_LIEN(LIEN);
  
  h5reader->get_GMI_nLocBas(nLocBas);

  nLocBas_max = *max_element(nLocBas.begin(), nLocBas.end());

  nLocBas_min = *min_element(nLocBas.begin(), nLocBas.end());

  
  elem_ptr.resize(nLocBas.size()+1);

  elem_ptr[0] = 0;

  for(unsigned int ii=1; ii<elem_ptr.size(); ++ii)
    elem_ptr[ii] = elem_ptr[ii-1] + nLocBas[ii-1];

  VEC_T::shrink2fit(elem_ptr);
  VEC_T::shrink2fit(LIEN);
  VEC_T::shrink2fit(nLocBas);

  if(LIEN.size() != elem_ptr[elem_ptr.size()-1])
  {
    SYS_T::commPrint("Error: LIEN length does not match with nLocBas. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}


ALocal_IEN_TS::~ALocal_IEN_TS()
{
  VEC_T::clean(LIEN);
  VEC_T::clean(nLocBas);
}



void ALocal_IEN_TS::get_LIEN(const int &ee, std::vector<int> &elem_ien) const
{
  elem_ien.resize(nLocBas[ee]);
  const int start_pt = elem_ptr[ee];
  for(int ii=0; ii<nLocBas[ee]; ++ii)
    elem_ien[ii] = LIEN[start_pt + ii];
}



void ALocal_IEN_TS::get_LIEN(const int &ee, int * &elem_ien) const
{
  const int start_pt = elem_ptr[ee];
  for(int ii=0; ii<nLocBas[ee]; ++ii)
    elem_ien[ii] = LIEN[start_pt + ii];
}



bool ALocal_IEN_TS::isNode_in_Elem(const int &elem, const int &node) const
{
  std::vector<int> temp;
  get_LIEN(elem, temp);
  std::vector<int>::const_iterator it = find(temp.begin(), temp.end(), node);
  return (it!=temp.end());
}



void ALocal_IEN_TS::print_info() const
{
  std::cout<<"ALocal_IEN: \n";
  VEC_T::print(LIEN);

  std::cout<<"nLocBas: "<<nLocBas_max<<'\t'<<nLocBas_min<<'\n';
  VEC_T::print(nLocBas);
}



// EOF
