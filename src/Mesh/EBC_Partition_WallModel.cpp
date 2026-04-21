#include "EBC_Partition_WallModel.hpp"
#include "HDF5_Group.hpp"
#include "HDF5_Writer.hpp"

EBC_Partition_WallModel::EBC_Partition_WallModel(const IPart * const &part,
    const Map_Node_Index * const &mnindex, const ElemBC * const &ebc)
: EBC_Partition(part, mnindex, ebc), wall_model_type {ebc->get_wall_model_type()}
{
  if(wall_model_type == 0)
    ;   // do nothing
  else if(wall_model_type == 1 || wall_model_type == 2)
  {
    for(int ee=0; ee < ebc->get_num_cell(0); ++ee)
    {
      const int global_vol_ele_id = ebc -> get_global_cell(0, ee);
      const int loc_id = part -> get_elemLocIndex(global_vol_ele_id);
      if( loc_id != -1 )
      {
        part_vol_ele_id.push_back( loc_id );
        ele_face_id.push_back( ebc -> get_faceID(ee) );
      }
    }

    SYS_T::print_fatal_if( num_local_cell[0] != VEC_T::get_size(part_vol_ele_id), "Error: EBC_Partition_WallModel the part_vol_ele_id vector size does not match with the number of local cell given in EBC_Partition.\n" ); 
  }
  else
    SYS_T::print_fatal("Error: EBC_Partition_WallModel, unknown wall model type.\n");
}

void EBC_Partition_WallModel::write_hdf5(const std::string &FileName) const
{
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  auto h5w = SYS_T::make_unique<HDF5_Writer>( fName, H5F_ACC_RDWR );
  const hid_t file_id = h5w->get_file_id();

  auto weak_group = HDF5_Group::create( file_id, "/weak" );

  h5w -> write_intScalar( weak_group.id(), "wall_model_type", wall_model_type );

  if(wall_model_type > 0)
  {
    h5w -> write_intScalar( weak_group.id(), "num_local_cell", num_local_cell[0] );

    h5w -> write_intVector( weak_group.id(), "part_volume_cell_id", part_vol_ele_id );

    h5w -> write_intVector( weak_group.id(), "cell_face_id", ele_face_id );
  }
  else
    ;   // stop writing if wall_model_type = 0

}

// EOF
