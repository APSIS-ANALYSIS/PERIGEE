#include "ALocal_WeakBC.hpp"

ALocal_WeakBC::ALocal_WeakBC( const std::string &fileBaseName, const int &cpu_rank, const double &in_C_bI)
: weakbc_type {0}, C_bI {in_C_bI}
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("weak");

  weakbc_type = h5r -> read_intScalar( gname.c_str(), "weakBC_type" );

  if (weakbc_type > 0)
  {
    num_sur_ele = h5r -> read_intScalar( gname.c_str(), "num_local_cell" );

    part_vol_ele_id = h5r -> read_intVector( gname.c_str(), "part_volume_cell_id" );

    ele_face_id = h5r -> read_intVector( gname.c_str(), "cell_face_id" );
  }

  delete h5r; H5Fclose( file_id );
}

ALocal_WeakBC::~ALocal_WeakBC()
{
  if(weakbc_type > 0)
  {
    VEC_T::clean( part_vol_ele_id );
    VEC_T::clean( ele_face_id );
  }
}