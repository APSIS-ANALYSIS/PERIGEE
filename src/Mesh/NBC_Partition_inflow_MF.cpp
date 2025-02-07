#include "NBC_Partition_inflow_MF.hpp"

NBC_Partition_inflow_MF::NBC_Partition_inflow_MF( const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const INodalBC * const &nbc,
    const std::vector< std::vector<int> > &grid2id )
: NBC_Partition_inflow(part, mnindex, nbc)
{
  SYS_T::print_fatal_if( VEC_T::get_size(grid2id) != 3, "Error: NBC_Partition_MF, the grid2id array size should be 3 in NBC_Partition_inflow_MF.\n" );

  LDN_MF.resize(num_nbc);

  // The inflow has 3 dofs
  for(int ii=0; ii<num_nbc; ++ii) LDN_MF[ii].resize(3*Num_LD[ii]);

  for(int ii=0; ii<num_nbc; ++ii)
  {
    for(int jj=0; jj<Num_LD[ii]; ++jj)
    {
      LDN_MF[ii][3*jj+0] = grid2id[0][ LDN[ii][jj] ];
      LDN_MF[ii][3*jj+1] = grid2id[1][ LDN[ii][jj] ];
      LDN_MF[ii][3*jj+2] = grid2id[2][ LDN[ii][jj] ];
    }
  }
}

void NBC_Partition_inflow_MF::write_hdf5( const std::string &FileName ) const
{
  // --------------------------------------------------------------------------
  // Call the base class writer to write the base class data
  NBC_Partition_inflow::write_hdf5( FileName );
  // --------------------------------------------------------------------------

  std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gopen(file_id, "/inflow", H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  for(int ii=0; ii<num_nbc; ++ii)
  {
    std::string subgroup_name( "nbcid_" );
    subgroup_name.append( std::to_string(ii) );

    hid_t group_id = H5Gopen(g_id, subgroup_name.c_str(), H5P_DEFAULT);

    h5w->write_intVector( group_id, "LDN_MF", LDN_MF[ii] );
  }

  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}

// EOF
