#include "ALocal_MovingBC.hpp"

ALocal_MovingBC::ALocal_MovingBC( 
    const std::string &fileBaseName, const int &cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/moving");

  num_nbc = h5r -> read_intScalar( gname.c_str(), "num_nbc" );

  // Allocate the size of the member data
  Num_LD.resize(num_nbc); 
  LDN.resize(num_nbc);
  LDN_pt_xyz.resize(num_nbc); 
  num_local_node.resize(num_nbc); 
  num_local_cell.resize(num_nbc); 
  cell_nLocBas.resize(num_nbc);
  local_pt_xyz.resize(num_nbc); 
  local_cell_ien.resize(num_nbc); 
  local_node_pos.resize(num_nbc);  

  std::string groupbase(gname);
  groupbase.append("/nbcid_");

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( std::to_string(nbc_id) );

    Num_LD[nbc_id] = h5r -> read_intScalar( subgroup_name.c_str(), "Num_LD" );

    // If this sub-domain of this CPU contains local moving bc points, load the LDN array.
    if( Num_LD[nbc_id] > 0 )
    {
      LDN[nbc_id]            = h5r -> read_intVector(    subgroup_name.c_str(), "LDN" );

      const std::vector<double> temp_ldn_xyz = h5r->read_doubleVector( subgroup_name.c_str(), "LDN_pt_xyz" );

      ASSERT( VEC_T::get_size(temp_ldn_xyz) == Num_LD[nbc_id]*3, "Error: ALocal_MovingBC LDN_pt_xyz format is wrong.\n");

      LDN_pt_xyz[nbc_id] = std::vector<Vector_3> (Num_LD[nbc_id], Vector_3{ 0, 0, 0 });
      
      for(int ii {0}; ii < Num_LD[nbc_id]; ++ii)
        LDN_pt_xyz[nbc_id][ii] = Vector_3{ temp_ldn_xyz[3 * ii], temp_ldn_xyz[3 * ii + 1], temp_ldn_xyz[3 * ii + 2] };
    }
    else
    {
      LDN[nbc_id].clear();
      LDN_pt_xyz[nbc_id].clear();
    }

    SYS_T::print_fatal_if( Num_LD[nbc_id] != static_cast<int>( LDN[nbc_id].size() ), "Error: the LDN vector size does not match with the value of Num_LD.\n" );

    // Cell-related data
    cell_nLocBas[nbc_id]   = h5r->read_intScalar( subgroup_name.c_str(), "cell_nLocBas" );
    num_local_cell[nbc_id] = h5r->read_intScalar( subgroup_name.c_str(), "num_local_cell" );
    num_local_node[nbc_id] = h5r->read_intScalar( subgroup_name.c_str(), "num_local_node" );

    // If this partitioned sub-domain contains moving surface element,
    // load its geometrical info 
    if(num_local_cell[nbc_id] > 0)
    {      
      local_cell_ien[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "local_cell_ien" );
    }
    else
    {
      local_cell_ien[nbc_id].clear();
    }

    if(num_local_node[nbc_id] > 0)
    {  
      const std::vector<double> temp_xyz = h5r->read_doubleVector( subgroup_name.c_str(), "local_pt_xyz" );

      ASSERT( VEC_T::get_size(temp_xyz) == num_local_node[nbc_id]*3, "Error: ALocal_MovingBC local_pt_xyz format is wrong.\n");

      local_pt_xyz[nbc_id] = std::vector<Vector_3> (num_local_node[nbc_id], Vector_3{ 0, 0, 0 });
        
      for(int ii {0}; ii < num_local_node[nbc_id]; ++ii)
        local_pt_xyz[nbc_id][ii] = Vector_3{ temp_xyz[3 * ii], temp_xyz[3 * ii + 1], temp_xyz[3 * ii + 2] };

      local_node_pos[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "local_node_pos" );
    }
    else
    {
      local_pt_xyz.clear();
      local_node_pos.clear();    
    }

  }// end nbc_id-loop

  delete h5r; H5Fclose( file_id );    
}

void ALocal_MovingBC::get_ctrlPts_xyz( const int &nbc_id,
    const int &eindex, double * const &ctrl_x, double * const &ctrl_y, 
    double * const &ctrl_z ) const
{
  for(int jj=0; jj<cell_nLocBas[nbc_id]; ++jj)
  {
    const int pos = local_cell_ien[nbc_id][ cell_nLocBas[nbc_id]*eindex+jj ];
    ctrl_x[jj] = local_pt_xyz[nbc_id][pos].x();
    ctrl_y[jj] = local_pt_xyz[nbc_id][pos].y();
    ctrl_z[jj] = local_pt_xyz[nbc_id][pos].z();
  }
}

void ALocal_MovingBC::get_SIEN( const int &nbc_id,
    const int &eindex, int * const &sien ) const
{
  for(int jj=0; jj<cell_nLocBas[nbc_id]; ++jj)
  {
    const int pos = local_cell_ien[nbc_id][ cell_nLocBas[nbc_id]*eindex+jj ];
    sien[jj] = local_node_pos[nbc_id][pos];
  }
}

std::vector<int> ALocal_MovingBC::get_SIEN( const int &nbc_id,
    const int &eindex ) const
{
  std::vector<int> out( cell_nLocBas[nbc_id], 0 );
  for(int jj=0; jj<cell_nLocBas[nbc_id]; ++jj)
  {
    const int pos = local_cell_ien[nbc_id][ cell_nLocBas[nbc_id]*eindex+jj ];
    out[jj] = local_node_pos[nbc_id][pos];
  }
  return out;
}

// EOF
