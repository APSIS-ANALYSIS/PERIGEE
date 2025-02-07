#include "ALocal_EBC.hpp"

ALocal_EBC::ALocal_EBC( const std::string &fileBaseName, 
    const int &cpu_rank, const std::string &gname )
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

  num_ebc = h5r -> read_intScalar(gname.c_str(), "num_ebc");

  if( num_ebc > 0)
  {
    num_local_cell_node = h5r -> read_intVector(gname.c_str(), "num_local_cell_node" );

    num_local_cell = h5r -> read_intVector(gname.c_str(), "num_local_cell" );

    cell_nLocBas = h5r -> read_intVector(gname.c_str(), "cell_nLocBas" );
  }

  std::string groupbase(gname);
  groupbase.append("/ebcid_");

  local_cell_node_xyz.resize(num_ebc);
  local_cell_ien.resize(num_ebc);
  local_cell_node_vol_id.resize(num_ebc);
  local_cell_node_pos.resize(num_ebc);
  local_cell_vol_id.resize(num_ebc);

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( std::to_string(ii) );

      local_cell_node_xyz[ii] = h5r -> read_doubleVector( subgroup_name.c_str(), "local_cell_node_xyz" );

      local_cell_ien[ii] = h5r -> read_intVector( subgroup_name.c_str(), "local_cell_ien" );

      local_cell_node_vol_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "local_cell_node_vol_id" );

      local_cell_node_pos[ii] = h5r -> read_intVector( subgroup_name.c_str(), "local_cell_node_pos" );

      local_cell_vol_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "local_cell_vol_id" );
    }
    else
    {
      local_cell_node_xyz[ii].clear();
      local_cell_ien[ii].clear();
      local_cell_node_vol_id[ii].clear();
      local_cell_node_pos[ii].clear();
      local_cell_vol_id[ii].clear();
    }
  }

  H5Fclose( file_id );
}

ALocal_EBC::ALocal_EBC( const HDF5_Reader * const &h5r,
    const std::string &gname )
{
  num_ebc = h5r -> read_intScalar(gname.c_str(), "num_ebc");

  if( num_ebc > 0)
  {
    num_local_cell_node = h5r -> read_intVector(gname.c_str(), "num_local_cell_node" );

    num_local_cell = h5r -> read_intVector(gname.c_str(), "num_local_cell" );

    cell_nLocBas = h5r -> read_intVector(gname.c_str(), "cell_nLocBas" );
  }

  std::string groupbase(gname);
  groupbase.append("/ebcid_");

  local_cell_node_xyz.resize(num_ebc);
  local_cell_ien.resize(num_ebc);
  local_cell_node_vol_id.resize(num_ebc);
  local_cell_node_pos.resize(num_ebc);
  local_cell_vol_id.resize(num_ebc);

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( std::to_string(ii) );

      local_cell_node_xyz[ii] = h5r -> read_doubleVector( subgroup_name.c_str(), "local_cell_node_xyz" );

      local_cell_ien[ii] = h5r -> read_intVector( subgroup_name.c_str(), "local_cell_ien" );

      local_cell_node_vol_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "local_cell_node_vol_id" );

      local_cell_node_pos[ii] = h5r -> read_intVector( subgroup_name.c_str(), "local_cell_node_pos" );

      local_cell_vol_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "local_cell_vol_id" );
    }
    else
    {
      local_cell_node_xyz[ii].clear();
      local_cell_ien[ii].clear();
      local_cell_node_vol_id[ii].clear();
      local_cell_node_pos[ii].clear();
      local_cell_vol_id[ii].clear();
    }
  }
}

void ALocal_EBC::get_ctrlPts_xyz(const int &ii,
    const int &eindex, double * const &ctrl_x,
    double * const &ctrl_y, double * const &ctrl_z ) const
{
  const int len = cell_nLocBas[ii];
  for(int jj=0; jj<len; ++jj)
  {
    const int pos = local_cell_ien[ii][len*eindex+jj];
    ctrl_x[jj] = local_cell_node_xyz[ii][3*pos];
    ctrl_y[jj] = local_cell_node_xyz[ii][3*pos+1];
    ctrl_z[jj] = local_cell_node_xyz[ii][3*pos+2];
  }
}

void ALocal_EBC::get_SIEN( const int &ii,
    const int &eindex, int * const &sien ) const
{
  const int len = cell_nLocBas[ii];
  for(int jj=0; jj<len; ++jj)
  {
    const int pos = local_cell_ien[ii][len*eindex+jj];
    sien[jj] = local_cell_node_pos[ii][pos];
  }
}

std::vector<int> ALocal_EBC::get_SIEN( const int &ii, const int &eindex ) const
{
  const int len = cell_nLocBas[ii];
  std::vector<int> out (len, 0);
  for(int jj=0; jj<len; ++jj)
  {
    const int pos = local_cell_ien[ii][len*eindex+jj];
    out[jj] = local_cell_node_pos[ii][pos];
  }
  return out;
}

// EOF
