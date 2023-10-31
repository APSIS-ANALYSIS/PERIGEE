#include "EBC_Partition_wall_turbulence.hpp"

EBC_Partition_wall_turbulence::EBC_Partition_wall_turbulence(const IPart * const &part,
    const Map_Node_Index * const &mnindex, const ElemBC * const &ebc)
: EBC_Partition(part, mnindex, ebc),
weak_bc_type {ebc->get_weak_bc_type()}, C_bI {ebc->get_C_bI()}
{
  if(weak_bc_type == 0)
    ;   // do nothing
  else if(weak_bc_type == 1 || weak_bc_type == 2)
  {
    for(int ee{0}; ee < get_num_local_cell(0); ++ee)
    {
      const int global_vol_ele_id = ebc -> get_global_cell(0, ee);
      const int loc_id = part -> get_elemLocIndex(global_vol_ele_id);
      if( loc_id != -1 )
      {
        part_vol_ele_id.push_back( loc_id );
        ele_face_id.push_back( ebc -> get_faceID(ee) );
      }
    }
  }
  else
    SYS_T::print_fatal("Error: EBC_Partition_wall_turbulence, unknown weak bc type.\n");
}

EBC_Partition_wall_turbulence::~EBC_Partition_wall_turbulence()
{
  part_vol_ele_id.clear();
  ele_face_id.clear();
}

void EBC_Partition_wall_turbulence::write_hdf5(const std::string &FileName) const
{
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, "weak", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_intScalar( g_id, "weakBC_type", weak_bc_type );

  if(weak_bc_type > 0)
  {
    h5w -> write_doubleScalar( g_id, "C_bI", C_bI );

    h5w -> write_intScalar( g_id, "num_local_cell", num_local_cell[0] );

    h5w -> write_intVector( g_id, "part_volume_cell_id", part_vol_ele_id );

    h5w -> write_intVector( g_id, "cell_face_id", ele_face_id );

  }
  else
    ;   // stop writing if weak_bc_type = 0

  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}