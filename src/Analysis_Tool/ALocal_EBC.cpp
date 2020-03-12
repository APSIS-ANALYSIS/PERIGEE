#include "ALocal_EBC.hpp"

ALocal_EBC::ALocal_EBC( const std::string &fileBaseName, 
    const int &cpu_rank, const std::string &gname )
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string ebcgroup_name(gname);

  num_ebc = h5r -> read_intScalar(ebcgroup_name.c_str(), "num_ebc");

  if( num_ebc > 0)
  {
    h5r -> read_intVector( ebcgroup_name.c_str(), "num_local_node",
        num_local_node );

    h5r -> read_intVector( ebcgroup_name.c_str(), "num_local_cell",
        num_local_cell );

    h5r -> read_intVector( ebcgroup_name.c_str(), "cell_nLocBas",
        cell_nLocBas );
  }

  std::string groupbase(gname);
  groupbase.append("/ebcid_");

  local_pt_xyz.resize(num_ebc);
  local_tri_ien.resize(num_ebc);
  local_global_node.resize(num_ebc);
  local_node_pos.resize(num_ebc);
  local_global_cell.resize(num_ebc);

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( SYS_T::to_string(ii) );

      h5r -> read_doubleVector( subgroup_name.c_str(), "local_pt_xyz",
          local_pt_xyz[ii] );

      h5r -> read_intVector( subgroup_name.c_str(), "local_tri_ien",
          local_tri_ien[ii] );

      h5r -> read_intVector( subgroup_name.c_str(), "local_global_node",
          local_global_node[ii] );

      h5r -> read_intVector( subgroup_name.c_str(), "local_node_pos",
          local_node_pos[ii] );

      h5r -> read_intVector( subgroup_name.c_str(), "local_global_cell", 
          local_global_cell[ii] );
    }
    else
    {
      local_pt_xyz[ii].clear();
      local_tri_ien[ii].clear();
      local_global_node[ii].clear();
      local_node_pos[ii].clear();
      local_global_cell[ii].clear();
    }
  }

  delete h5r; H5Fclose( file_id );
}


ALocal_EBC::~ALocal_EBC()
{
  VEC_T::clean( num_local_node );
  VEC_T::clean( num_local_cell );
  VEC_T::clean( cell_nLocBas );
  VEC_T::clean( local_pt_xyz );
  VEC_T::clean( local_tri_ien );
  VEC_T::clean( local_global_node );
  VEC_T::clean( local_node_pos );
  VEC_T::clean( local_global_cell );
}


void ALocal_EBC::get_ctrlPts_xyz(const int &ii,
    const int &eindex, double * const &ctrl_x,
    double * const &ctrl_y, double * const &ctrl_z ) const
{
  const int len = cell_nLocBas[ii];
  for(int jj=0; jj<len; ++jj)
  {
    const int pos = local_tri_ien[ii][len*eindex+jj];
    ctrl_x[jj] = local_pt_xyz[ii][3*pos];
    ctrl_y[jj] = local_pt_xyz[ii][3*pos+1];
    ctrl_z[jj] = local_pt_xyz[ii][3*pos+2];
  }
}


void ALocal_EBC::get_SIEN( const int &ii,
    const int &eindex, int * const &sien ) const
{
  const int len = cell_nLocBas[ii];
  for(int jj=0; jj<len; ++jj)
  {
    const int pos = local_tri_ien[ii][len*eindex+jj];
    sien[jj] = local_node_pos[ii][pos];
  }
}


// EOF
