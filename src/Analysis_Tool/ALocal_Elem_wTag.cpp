#include "ALocal_Elem_wTag.hpp"

ALocal_Elem_wTag::ALocal_Elem_wTag( const std::string &fileBaseName, 
    const int &cpu_rank )
: ALocal_Elem( fileBaseName, cpu_rank )
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  elem_tag = h5r->read_intVector("/Local_Elem", "elem_phy_tag");

  delete h5r; H5Fclose( file_id );
}

ALocal_Elem_wTag::~ALocal_Elem_wTag()
{
  VEC_T::clean(elem_tag);
}

void ALocal_Elem_wTag::print_info() const
{
  std::cout<<"nlocalelem: "<<get_nlocalele()<<std::endl;
  std::cout<<"elem_loc(elem_tag): \n";
  for(int ii=0; ii<get_nlocalele(); ++ii)
    std::cout<<get_elem_loc(ii)<<"("<<get_elem_tag(ii)<<")"<<'\t';
  std::cout<<std::endl;
}

// EOF
