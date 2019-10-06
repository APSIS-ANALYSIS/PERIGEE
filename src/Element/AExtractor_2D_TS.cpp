#include "AExtractor_2D_TS.hpp"

AExtractor_2D_TS::AExtractor_2D_TS( const HDF5_PartReader * const &h5reader )
{
  h5reader->get_EXT_TS_full( extractor );
  h5reader->get_GMI_degree(sdegree, tdegree);
  
  sp1tp1 = (sdegree+1)*(tdegree+1);

  std::vector<int> nlocbas;

  h5reader->get_GMI_nLocBas(nlocbas);

  elem_ptr.resize(nlocbas.size() + 1);
  
  elem_ptr[0] = 0;

  for(unsigned int ii=1; ii<elem_ptr.size(); ++ii)
    elem_ptr[ii] = elem_ptr[ii-1] + nlocbas[ii-1];

  VEC_T::shrink2fit(elem_ptr);
  
  if( extractor.size() != sp1tp1 * elem_ptr[elem_ptr.size()-1] )
  {
    SYS_T::commPrint("Error: the extractor size does not match with nLocBas, sDegree, tDegree. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

}


AExtractor_2D_TS::AExtractor_2D_TS( const int &in_sdeg, const int &in_tdeg,
            const std::vector<double> &in_ext, const std::vector<unsigned int> in_elem_ptr )
{
  sdegree = in_sdeg;
  tdegree = in_tdeg;

  sp1tp1 = (sdegree + 1) * (tdegree + 1);

  extractor = in_ext;
  VEC_T::shrink2fit(extractor);

  elem_ptr = in_elem_ptr;
  VEC_T::shrink2fit(elem_ptr);

  if( extractor.size() != sp1tp1 * elem_ptr[elem_ptr.size()-1] )
  {
    SYS_T::commPrint("Error: the extractor size does not match with nLocBas, sDegree, tDegree. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

}





AExtractor_2D_TS::~AExtractor_2D_TS()
{
  VEC_T::clean(extractor);
  VEC_T::clean(elem_ptr);
}


void AExtractor_2D_TS::get_EXT( const int &ee, const int &ii, std::vector<double> &ext ) const
{
  ext.resize(sp1tp1);
  const int start = ( elem_ptr[ee] + ii ) * sp1tp1;

  for(int ii=0; ii<sp1tp1; ++ii)
    ext[ii] = extractor[start + ii];
}



void AExtractor_2D_TS::get_EXT( const int &ee, const int &ii, double * &ext ) const
{
  ext = new double [sp1tp1];
  const int start = ( elem_ptr[ee] + ii )* sp1tp1;

  for(int ii=0; ii<sp1tp1; ++ii)
    ext[ii] = extractor[start + ii];
}



void AExtractor_2D_TS::print_info() const
{
  std::cout<<"sdegree = "<<sdegree<<std::endl;
  std::cout<<"tdegree = "<<tdegree<<std::endl;
  std::cout<<"extractor: \n";
  VEC_T::print(extractor);
  std::cout<<"\n elem_ptr : \n";
  VEC_T::print(elem_ptr);
}










//EOF
