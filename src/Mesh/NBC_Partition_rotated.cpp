#include "NBC_Partition_rotated.hpp"

NBC_Partition_rotated::NBC_Partition_rotated(
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const INodalBC * const &nbc ) 
: cpu_rank(part->get_cpu_rank())
{  
  // Collect the Dirichlet nodes on the rotated boundary surface
  LDN.clear();
  Num_LD = 0;
  for(unsigned int jj=0; jj<nbc->get_num_dir_nodes_on_rotated_surface(); ++jj)
  {
    int node_index = static_cast<int>( nbc -> get_dir_nodes_on_rotated_surface(jj) );
    node_index = mnindex -> get_old2new(node_index);
    if(part->isNodeInPart(node_index))
    {
      LDN.push_back(node_index);
      Num_LD += 1;
    }
  }

  // Record the geometrical info of the rotated surface in this CPU
  cell_nLocBas = nbc -> get_nLocBas();

  std::vector<int> local_node {}, local_elem {};

  local_global_cell.clear();

  // Loop over all cells on rotated boundary surface
  for(int jj=0; jj < nbc->get_num_cell(); ++jj)
  {
    const int elem_index = nbc -> get_global_cell(jj); // cell vol id

    if( part -> get_elemLocIndex( elem_index ) != -1  )
    {
      local_elem.push_back( jj );

      local_global_cell.push_back( elem_index );

      for( int kk=0; kk<cell_nLocBas; ++kk )
        local_node.push_back( nbc->get_ien(jj, kk) );
    }
  }

  VEC_T::sort_unique_resize( local_node );

  num_local_node = static_cast<int>( local_node.size() );
  num_local_cell = static_cast<int>( local_global_cell.size() );

  local_pt_xyz.resize(      num_local_node * 3 );
  local_global_node.resize( num_local_node );
  local_node_pos.resize(    num_local_node );

  for(int jj=0; jj<num_local_node; ++jj)
  {
    local_pt_xyz[3*jj+0] = nbc -> get_pt_xyz( local_node[jj], 0 );
    local_pt_xyz[3*jj+1] = nbc -> get_pt_xyz( local_node[jj], 1 );
    local_pt_xyz[3*jj+2] = nbc -> get_pt_xyz( local_node[jj], 2 );

    local_global_node[jj] = nbc -> get_global_node( local_node[jj] );

    local_node_pos[jj] = part->get_nodeLocGhoIndex( mnindex->get_old2new(local_global_node[jj]) );
  }

  LDN_pt_xyz.resize( Num_LD * 3 );

  for(int jj=0; jj<Num_LD; ++jj)
  {
    const int LDN_old_index = mnindex->get_new2old(LDN[jj]);
    
    // LDN_old_pos: the position of old_LDN_index in local_global_node
    int LDN_old_pos = VEC_T::get_pos(local_global_node, LDN_old_index);

    LDN_pt_xyz[3*jj+0] = local_pt_xyz[ 3 * LDN_old_pos + 0 ]; 
    LDN_pt_xyz[3*jj+1] = local_pt_xyz[ 3 * LDN_old_pos + 1 ]; 
    LDN_pt_xyz[3*jj+2] = local_pt_xyz[ 3 * LDN_old_pos + 2 ]; 
  }

  // create new IEN
  local_cell_ien.resize( num_local_cell * cell_nLocBas );

  for(int jj=0; jj<num_local_cell; ++jj)
  {
    for(int kk=0; kk<cell_nLocBas; ++kk)
    {
      const int temp_node = nbc -> get_ien( local_elem[jj], kk );
      const int temp_npos = VEC_T::get_pos( local_node, temp_node );
      local_cell_ien[jj*cell_nLocBas + kk] = temp_npos;
    }
  }
}

void NBC_Partition_rotated::write_hdf5( const std::string &FileName ) const
{
  std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, "/rotated_nbc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  h5w->write_intScalar( g_id, "Num_LD", Num_LD );

  h5w->write_intVector( g_id, "LDN", LDN );

  h5w->write_intScalar( g_id, "num_local_node", num_local_node );

  h5w->write_intScalar( g_id, "num_local_cell", num_local_cell );

  h5w->write_intScalar( g_id, "cell_nLocBas", cell_nLocBas );

  h5w->write_doubleVector( g_id, "local_pt_xyz", local_pt_xyz );

  h5w->write_doubleVector( g_id, "LDN_pt_xyz", LDN_pt_xyz );

  h5w->write_intVector( g_id, "local_cell_ien", local_cell_ien );

  h5w->write_intVector( g_id, "local_global_node", local_global_node );

  h5w->write_intVector( g_id, "local_node_pos", local_node_pos );

  h5w->write_intVector( g_id, "local_global_cell", local_global_cell );

  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}

// EOF
