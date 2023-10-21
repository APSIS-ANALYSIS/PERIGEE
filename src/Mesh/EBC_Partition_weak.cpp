#include "EBC_Partition_weak.hpp"

EBC_Partition_weak::EBC_Partition_weak(const IPart * const &part,
    const Map_Node_Index * const &mnindex, const ElemBC * const &ebc)
: EBC_Partition(part, mnindex, ebc),
weak_bc_type {ebc->get_weak_bc_type()}, C_bI {ebc->get_C_bI()}
{
  if(weak_bc_type == 0)
    ;   // do nothing
  else if(weak_bc_type == 1 || weak_bc_type == 2)
  {

    vol_ele_id.resize(num_ebc);
    ele_face_id.resize(num_ebc);
    if(weak_bc_type == 2)
      rot_mat.resize(num_ebc);

    for(int ii{0}; ii < num_ebc; ++ii)
    {   
      vol_ele_id[ii].resize(get_num_local_cell(ii));
      for(int ee{0}; ee < get_num_local_cell(ii); ++ee)
        vol_ele_id[ii][ee] = ebc -> get_global_cell(ii, ee);

      ele_face_id[ii] = ebc -> get_faceID(ii);

      if(weak_bc_type == 2)
        rot_mat[ii] = ebc -> get_rotation_matrix(ii);
    }
  }
  else
    SYS_T::print_fatal("Error: EBC_Partition_weak, unknown weak bc type.\n");

  if(weak_bc_type == 2)
    ;   // Rotation matrix, unimplemented
}

EBC_Partition_weak::~EBC_Partition_weak()
{
  vol_ele_id.clear();
  ele_face_id.clear();
}

void EBC_Partition_weak::write_hdf5(const std::string &FileName) const
{
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, "weak", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_intScalar( g_id, "weakBC_type", weak_bc_type );

  h5w -> write_intScalar( g_id, "num_weakBC", num_ebc );

  if(weak_bc_type > 0)
  {
    h5w -> write_doubleScalar( g_id, "C_bI", C_bI );

    h5w -> write_intVector( g_id, "num_local_cell", num_local_cell );

    if(weak_bc_type == 2)
      h5w -> write_intVector( g_id, "num_local_cell_node", num_local_cell_node );

    const std::string groupbase("weakBCid_");

    for(int ii{0}; ii < num_ebc; ++ii)
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( std::to_string(ii) );

      hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(), 
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      h5w -> write_intVector( group_id, "volume_cell_id", vol_ele_id[ii] );

      h5w -> write_intVector( group_id, "cell_face_id", ele_face_id[ii] );

      if(weak_bc_type == 2)
      {
        h5w -> write_intVector( group_id, "global_node_id", local_cell_node_vol_id[ii] );

        h5w -> write_doubleVector( group_id, "node_rotation_matrix", rot_mat[ii] );
      }
      
      H5Gclose( group_id );
    }
  }
  else
    ;   // stop writing if weak_bc_type = 0

  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}