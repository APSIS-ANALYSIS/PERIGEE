#include "EBC_Partition.hpp"

EBC_Partition::EBC_Partition( const IPart * const &part,
    const Map_Node_Index * const &mnindex, const ElemBC * const &ebc )
: cpu_rank( part->get_cpu_rank() ), num_ebc( ebc->get_num_ebc() )
{
  // Clean up all the vectors
  num_local_cell_node.clear();
  num_local_cell.clear();
  cell_nLocBas.clear();
  local_cell_node_xyz.clear();
  local_cell_ien.clear();
  local_cell_node_vol_id.clear();
  local_cell_node.clear();
  local_cell_node_pos.clear();
  local_cell_vol_id.clear();

  // resize the vectors with length equaling to the number of ebc's
  num_local_cell_node.resize( num_ebc );
  num_local_cell.resize( num_ebc );
  cell_nLocBas.resize( num_ebc );
  local_cell_node_xyz.resize( num_ebc );
  local_cell_ien.resize( num_ebc );
  local_cell_node_vol_id.resize( num_ebc );
  local_cell_node.resize( num_ebc );
  local_cell_node_pos.resize( num_ebc );
  local_cell_vol_id.resize( num_ebc );

  // fill the cell_nLocBas array
  for(int ii=0; ii<num_ebc; ++ii) 
    cell_nLocBas[ii] = ebc->get_cell_nLocBas(ii);

  // Loop over each element bc surface
  for(int ii=0; ii<num_ebc; ++ii)
  {
    // empty the ii-th surface's local_cell_node array
    local_cell_node[ii].clear(); 
    
    // obtain the node belonging to this partition
    std::vector<int> local_elem;
    local_elem.clear();

    const int num_global_bccell = ebc->get_num_cell( ii );

    local_cell_vol_id[ii].clear();

    PERIGEE_OMP_PARALLEL
    {
      std::vector<int> temp_local_elem {};
      std::vector<int> temp_local_cell_vol_id {};
      std::vector<int> temp_local_cell_node {};

      PERIGEE_OMP_FOR
      for(int jj=0; jj<num_global_bccell; ++jj)
      {
        const int elem_index = ebc -> get_global_cell(ii, jj);

        // If the element belongs to the partitioned subdomain, record it
        if( part -> get_elemLocIndex( elem_index ) != -1 )
        {
          temp_local_elem.push_back( jj );

          // record the cell's global (volumetric) element index
          temp_local_cell_vol_id.push_back( elem_index );

          // now put this cell's nodes into local_cell_node list
          for(int kk=0; kk<cell_nLocBas[ii]; ++kk)
            temp_local_cell_node.push_back( ebc->get_ien(ii, jj, kk) );
        }
      }

      PERIGEE_OMP_CRITICAL
      {
        VEC_T::insert_end(local_elem, temp_local_elem);
        VEC_T::insert_end(local_cell_vol_id[ii], temp_local_cell_vol_id);
        VEC_T::insert_end(local_cell_node[ii], temp_local_cell_node);
      }
    }

    VEC_T::sort_unique_resize( local_cell_node[ii] );

    // local_cell_node[ii] now stores indices (in surface mesh ii) for all cell nodes
    // in this partition of the surface; 
    // local_cell_vol_id[ii] stores all cells in this partition of the surface.
    num_local_cell_node[ii] = static_cast<int>( local_cell_node[ii].size() );
    num_local_cell[ii] = static_cast<int>( local_cell_vol_id[ii].size() );

    // For each local cell node, extract nodal xyz-coordinates and the global
    // volumetric mesh index.
    local_cell_node_xyz[ii].resize(num_local_cell_node[ii] * 3);
    local_cell_node_vol_id[ii].resize( num_local_cell_node[ii] );
    local_cell_node_pos[ii].resize( num_local_cell_node[ii] );

    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_local_cell_node[ii]; ++jj)
    {
      local_cell_node_xyz[ii][3*jj]   = ebc->get_pt_xyz( ii, local_cell_node[ii][jj], 0 );
      local_cell_node_xyz[ii][3*jj+1] = ebc->get_pt_xyz( ii, local_cell_node[ii][jj], 1 );
      local_cell_node_xyz[ii][3*jj+2] = ebc->get_pt_xyz( ii, local_cell_node[ii][jj], 2 );
      local_cell_node_vol_id[ii][jj] = ebc->get_global_node( ii, local_cell_node[ii][jj] );
      local_cell_node_pos[ii][jj] = part->get_nodeLocGhoIndex( mnindex->get_old2new( local_cell_node_vol_id[ii][jj] ) );
      
      SYS_T::print_fatal_if( local_cell_node_pos[ii][jj] < 0, "ERROR: there are nodes on the ebc surface not found in the partition's localtogloal array.\n" );
    }

    // now create the new IEN
    local_cell_ien[ii].resize( num_local_cell[ii] * cell_nLocBas[ii] );
    
    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_local_cell[ii]; ++jj)
    {
      for(int kk=0; kk<cell_nLocBas[ii]; ++kk)
      {
        const int temp_node = ebc->get_ien(ii, local_elem[jj], kk);
        const int temp_npos = VEC_T::get_pos( local_cell_node[ii], temp_node );
        SYS_T::print_fatal_if( temp_npos < 0, "Error: EBC_Partition, local_cell_node is incomplete. \n" );
        local_cell_ien[ii][jj*cell_nLocBas[ii] + kk] = temp_npos;
      }
    }
  } // end loop over num_ebc
}

void EBC_Partition::write_hdf5( const std::string &FileName, 
    const std::string &GroupName ) const
{
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, GroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_intScalar( g_id, "num_ebc", num_ebc );

  h5w -> write_intVector( g_id, "num_local_cell_node", num_local_cell_node );

  h5w -> write_intVector( g_id, "num_local_cell", num_local_cell );

  h5w -> write_intVector( g_id, "cell_nLocBas", cell_nLocBas );

  const std::string groupbase("ebcid_");

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( std::to_string(ii) );

      hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(), 
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      h5w->write_doubleVector( group_id, "local_cell_node_xyz", local_cell_node_xyz[ii] );

      h5w->write_intVector( group_id, "local_cell_ien", local_cell_ien[ii] );

      h5w->write_intVector( group_id, "local_cell_node_vol_id", local_cell_node_vol_id[ii] );

      h5w->write_intVector( group_id, "local_cell_node_pos", local_cell_node_pos[ii] );

      h5w->write_intVector( group_id, "local_cell_vol_id", local_cell_vol_id[ii] );

      H5Gclose( group_id );
    }
  }

  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}

void EBC_Partition::print_info() const
{
  std::cout<<"=========================================== \n";
  std::cout<<"EBC_Partition : \n";
  std::cout<<"-- num_ebc = "<<num_ebc<<std::endl;
  std::cout<<"-- num_local_cell_node : ";
  VEC_T::print(num_local_cell_node);
  std::cout<<"-- num_local_cell : ";
  VEC_T::print(num_local_cell);
  std::cout<<"-- cell_nLocBas : ";
  VEC_T::print(cell_nLocBas);
  for(int ii=0; ii<num_ebc; ++ii)
  {
    std::cout<<"-- ebc_id = "<<ii<<std::endl;
    std::cout<<"   local_cell_node_vol_id : \n";
    VEC_T::print( local_cell_node_vol_id[ii] );
    std::cout<<"   local_cell_node_xyz : \n";
    VEC_T::print( local_cell_node_xyz[ii] );
    std::cout<<"   local_cell_vol_id : \n";
    VEC_T::print( local_cell_vol_id[ii] );
    std::cout<<"   local_cell_ien : \n";
    VEC_T::print( local_cell_ien[ii] );
  } 
  std::cout<<"=========================================== \n";
}

// EOF
