#include "EBC_Partition_vtp_wall.hpp"

EBC_Partition_vtp_wall::EBC_Partition_vtp_wall( 
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const ElemBC * const &ebc )
: EBC_Partition_vtp(part, mnindex, ebc)
{
  const int ebc_id = 0;
  if( num_local_cell[ebc_id] > 0 )
  {
    ebc -> get_wall_thickness( part_thickness );
    ebc -> get_wall_youngsmod( part_youngsmod );
  }
  else
  {
    part_thickness.clear();
    part_youngsmod.clear();
  }
}


EBC_Partition_vtp_wall::~EBC_Partition_vtp_wall()
{
  VEC_T::clean( part_thickness );
  VEC_T::clean( part_youngsmod );
}


void EBC_Partition_vtp_wall::write_hdf5( const char * FileName ) const
{
  // Call base class writer to write base class data
  EBC_Partition_vtp::write_hdf5( FileName );

  const std::string input_fName(FileName);
  const std::string fName = SYS_T::gen_partfile_name( input_fName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t g_id = H5Gopen( file_id, "/ebc", H5P_DEFAULT );

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  const std::string groupbase("ebcid_");

  // num_ebc = 1 for wall elem bc
  const int ebc_id = 0;
  if( num_local_cell[ebc_id] > 0 )
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( SYS_T::to_string(ebc_id) );
    hid_t subgroup_id = H5Gopen(g_id, subgroup_name.c_str(), H5P_DEFAULT );

    h5w->write_doubleVector( subgroup_id, "thickness", part_thickness );
    h5w->write_doubleVector( subgroup_id, "youngsmod", part_youngsmod );

    H5Gclose( subgroup_id );
  }

  delete h5w;
  H5Gclose( g_id );
  H5Fclose( file_id );
}


void EBC_Partition_vtp_wall::write_hdf5( const char * FileName, 
    const char * GroupName ) const
{
  // Call the base class writer to write the base class data
  EBC_Partition_vtp::write_hdf5( FileName, GroupName );

  const std::string input_fName(FileName);
  const std::string fName = SYS_T::gen_partfile_name( input_fName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t g_id = H5Gopen( file_id, GroupName, H5P_DEFAULT );

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  const std::string groupbase("ebcid_");

  // num_ebc = 1 for wall elem bc
  const int ebc_id = 0;
  if( num_local_cell[ebc_id] > 0 )
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( SYS_T::to_string(ebc_id) );
    hid_t subgroup_id = H5Gopen(g_id, subgroup_name.c_str(), H5P_DEFAULT );

    h5w->write_doubleVector( subgroup_id, "thickness", part_thickness );
    h5w->write_doubleVector( subgroup_id, "youngsmod", part_youngsmod );

    H5Gclose( subgroup_id );
  }

  delete h5w;
  H5Gclose( g_id );
  H5Fclose( file_id );
}

// EOF
