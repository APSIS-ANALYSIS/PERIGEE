#include "ALocal_WeakBC.hpp"

ALocal_WeakBC::ALocal_WeakBC( const std::string &fileBaseName, 
    const int &cpu_rank )
: wall_model_type {0}
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("weak");

  wall_model_type = h5r -> read_intScalar( gname.c_str(), "wall_model_type" );

  if (wall_model_type > 0)
  {
    num_sur_ele = h5r -> read_intScalar( gname.c_str(), "num_local_cell" );

    part_vol_ele_id = h5r -> read_intVector( gname.c_str(), "part_volume_cell_id" );

    ele_face_id = h5r -> read_intVector( gname.c_str(), "cell_face_id" );
  }
  else
  {
    num_sur_ele = 0;
    part_vol_ele_id = {};
    ele_face_id = {};
  }

  delete h5r; H5Fclose( file_id );
}

ALocal_WeakBC::~ALocal_WeakBC()
{
  if(wall_model_type > 0)
  {
    VEC_T::clean( part_vol_ele_id );
    VEC_T::clean( ele_face_id );
  }
}

void ALocal_WeakBC::print_info() const
{
  SYS_T::commPrint("Weakly enforced wall boundary condition:\n");
  SYS_T::commPrint("-wall_model_type: %d\n", wall_model_type);
}

// EOF
