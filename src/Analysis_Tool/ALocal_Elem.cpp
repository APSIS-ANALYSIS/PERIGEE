#include "ALocal_Elem.hpp"

ALocal_Elem::ALocal_Elem(const std::string &fileBaseName, const int &cpu_rank)
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  elem_loc = h5r->read_intVector( "/Local_Elem", "elem_loc" );

  nlocalele = h5r->read_intScalar( "/Local_Elem", "nlocalele" );

  delete h5r; H5Fclose( file_id );
}

ALocal_Elem::~ALocal_Elem()
{
  VEC_T::clean(elem_loc);
}

void ALocal_Elem::print_info() const
{
  std::cout<<"nlocalelem: "<<nlocalele<<std::endl;
  std::cout<<"elem_loc: \n";
  VEC_T::print(elem_loc);
}

// EOF
