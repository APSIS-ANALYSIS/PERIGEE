#include "ALocal_WeakBC.hpp"

ALocal_WeakBC::ALocal_WeakBC( const std::string &fileBaseName, const int &cpu_rank)
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/weak");

  weakbc_type = h5r -> read_intScalar( gname.c_str(), "weak_bc_type" );

  if (weakbc_type > 0)
  {
    
  }

  if(weakbc_type == 2)
  {

  }

}