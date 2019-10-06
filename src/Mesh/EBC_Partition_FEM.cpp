#include "EBC_Partition_FEM.hpp"

EBC_Partition_FEM::EBC_Partition_FEM( const IPart * const &part,
    const Map_Node_Index * const &mnindex,     
    const ElemBC * const &ebc )
: cpu_rank( part->get_cpu_rank() ), num_ebc( ebc->get_num_ebc() )
{
  // Clean up all the vectors
  num_local_node.clear();
  num_local_cell.clear();
  cell_nLocBas.clear();
  local_pt_xyz.clear();
  local_tri_ien.clear();
  local_global_node.clear();
  local_node_pos.clear();
  local_global_cell.clear();
  local_intpt.clear();

  // resize the vectors with length equaling to the number of ebc's
  num_local_node.resize( num_ebc );
  num_local_cell.resize( num_ebc );
  cell_nLocBas.resize( num_ebc );
  local_pt_xyz.resize( num_ebc );
  local_tri_ien.resize( num_ebc );
  local_global_node.resize( num_ebc );
  local_node_pos.resize( num_ebc );
  local_global_cell.resize( num_ebc );
  local_intpt.resize( num_ebc );

  // fill the cell_nLocBas array
  for(int ii=0; ii<num_ebc; ++ii) 
    cell_nLocBas[ii] = ebc->get_cell_nLocBas(ii);

  // obtain the node belonging to this partition
  std::vector<int> local_node, local_elem;

  int elem_index, elem_pos; 
  for(int ii=0; ii<num_ebc; ++ii)
  {
    local_node.clear(); local_elem.clear();

    const int num_global_bccell = ebc->get_num_cell( ii );

    local_global_cell[ii].clear();

    for(int jj=0; jj<num_global_bccell; ++jj)
    {
      elem_index = ebc -> get_global_cell(ii, jj);

      // Determine if the element is in the partitioned subdomain
      elem_pos = part -> get_elemLocIndex( elem_index );

      if( elem_pos != -1 )
      {
        local_elem.push_back( jj );

        // If this element is in this partition, add it to the local
        // cell's global cell index list
        local_global_cell[ii].push_back( elem_index );

        // now put this cell's nodes into local_node list
        for(int kk=0; kk<cell_nLocBas[ii]; ++kk)
          local_node.push_back( ebc->get_ien(ii, jj, kk) );
      }
    }
    VEC_T::sort_unique_resize( local_node );

    // local_node now stores all the node that belong to this portion of the
    // surface; local_global_cell[ii] stores all the cell that belong to the
    // portion of the surface.
    num_local_node[ii] = static_cast<int>( local_node.size() );
    num_local_cell[ii] = static_cast<int>( local_global_cell[ii].size() );

    // local_node[jj] gives the node index on the surface domain. We extract
    // their xyz-coordinates and their global volumetric mesh indicies now.
    local_pt_xyz[ii].resize(num_local_node[ii] * 3);
    local_global_node[ii].resize( num_local_node[ii] );
    local_node_pos[ii].resize( num_local_node[ii] );
    for(int jj=0; jj<num_local_node[ii]; ++jj)
    {
      local_pt_xyz[ii][3*jj]   = ebc->get_pt_xyz( ii, local_node[jj], 0 );
      local_pt_xyz[ii][3*jj+1] = ebc->get_pt_xyz( ii, local_node[jj], 1 );
      local_pt_xyz[ii][3*jj+2] = ebc->get_pt_xyz( ii, local_node[jj], 2 );
      local_global_node[ii][jj] = ebc->get_global_node( ii, local_node[jj] );
      local_node_pos[ii][jj] = part->get_nodeLocGhoIndex( mnindex->get_old2new( local_global_node[ii][jj] ) );
      assert(local_node_pos[ii][jj] >= 0);
    }

    // now create the new IEN & cell interior point
    local_tri_ien[ii].resize( num_local_cell[ii] * cell_nLocBas[ii] );
    local_intpt[ii].resize( num_local_cell[ii] * 3 );
    for(int jj=0; jj<num_local_cell[ii]; ++jj)
    {
      local_intpt[ii][jj*3+0] = ebc->get_intpt_xyz(ii, local_elem[jj], 0);
      local_intpt[ii][jj*3+1] = ebc->get_intpt_xyz(ii, local_elem[jj], 1);
      local_intpt[ii][jj*3+2] = ebc->get_intpt_xyz(ii, local_elem[jj], 2);

      for(int kk=0; kk<cell_nLocBas[ii]; ++kk)
      {
        int temp_node = ebc->get_ien(ii, local_elem[jj], kk);
        int temp_npos = VEC_T::get_pos( local_node, temp_node );
        SYS_T::print_exit_if( temp_npos < 0, 
            "Error: EBC_Partition_FEM, local_node is incomplete. \n" );
        local_tri_ien[ii][jj*cell_nLocBas[ii] + kk] = temp_npos;
      }
    }
  }// end-loop-over-ebc
}


EBC_Partition_FEM::~EBC_Partition_FEM()
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


void EBC_Partition_FEM::write_hdf5( const char * FileName ) const
{
  std::string fName(FileName);
  fName.append("_p");

  if( cpu_rank / 10 == 0 )
    fName.append("0000");
  else if( cpu_rank / 100 == 0 )
    fName.append("000");
  else if( cpu_rank / 1000 == 0 )
    fName.append("00");
  else if( cpu_rank / 10000 == 0 )
    fName.append("0");

  std::stringstream sstrm;
  sstrm<<cpu_rank;
  fName.append(sstrm.str());

  fName.append(".h5");

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, "/ebc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_intScalar( g_id, "num_ebc", num_ebc );

  h5w -> write_intVector( g_id, "num_local_node", num_local_node );

  h5w -> write_intVector( g_id, "num_local_cell", num_local_cell );

  h5w -> write_intVector( g_id, "cell_nLocBas", cell_nLocBas );

  const std::string groupbase("ebcid_");

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( SYS_T::to_string(ii) );

      hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(), 
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      h5w->write_doubleVector( group_id, "local_pt_xyz", local_pt_xyz[ii] );

      h5w->write_intVector( group_id, "local_tri_ien", local_tri_ien[ii] );

      h5w->write_intVector( group_id, "local_global_node", local_global_node[ii] );
      
      h5w->write_intVector( group_id, "local_node_pos", local_node_pos[ii] );

      h5w->write_intVector( group_id, "local_global_cell", local_global_cell[ii] );

      h5w->write_doubleVector( group_id, "local_intpt_xyz", local_intpt[ii] );

      H5Gclose( group_id );
    }
  }

  delete h5w;
  H5Gclose( g_id );
  H5Fclose( file_id );
}


void EBC_Partition_FEM::write_hdf5( const char * FileName,
   const char * GroupName ) const
{
  std::string fName(FileName);
  fName.append("_p");

  if( cpu_rank / 10 == 0 )
    fName.append("0000");
  else if( cpu_rank / 100 == 0 )
    fName.append("000");
  else if( cpu_rank / 1000 == 0 )
    fName.append("00");
  else if( cpu_rank / 10000 == 0 )
    fName.append("0");

  std::stringstream sstrm;
  sstrm<<cpu_rank;
  fName.append(sstrm.str());

  fName.append(".h5");

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, GroupName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_intScalar( g_id, "num_ebc", num_ebc );

  h5w -> write_intVector( g_id, "num_local_node", num_local_node );

  h5w -> write_intVector( g_id, "num_local_cell", num_local_cell );

  h5w -> write_intVector( g_id, "cell_nLocBas", cell_nLocBas );

  const std::string groupbase("ebcid_");

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( SYS_T::to_string(ii) );

      hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(), 
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      h5w->write_doubleVector( group_id, "local_pt_xyz", local_pt_xyz[ii] );

      h5w->write_intVector( group_id, "local_tri_ien", local_tri_ien[ii] );

      h5w->write_intVector( group_id, "local_global_node", local_global_node[ii] );
      
      h5w->write_intVector( group_id, "local_node_pos", local_node_pos[ii] );

      h5w->write_intVector( group_id, "local_global_cell", local_global_cell[ii] );

      h5w->write_doubleVector( group_id, "local_intpt_xyz", local_intpt[ii] );
      
      H5Gclose( group_id );
    }
  }

  delete h5w;
  H5Gclose( g_id );
  H5Fclose( file_id );
}


void EBC_Partition_FEM::print_info() const
{
  std::cout<<"=========================================== \n";
  std::cout<<"EBC_Partition_FEM : \n";
  std::cout<<"-- num_ebc = "<<num_ebc<<std::endl;
  std::cout<<"-- num_local_node : ";
  VEC_T::print(num_local_node);
  std::cout<<"-- num_local_cell : ";
  VEC_T::print(num_local_cell);
  std::cout<<"-- cell_nLocBas : ";
  VEC_T::print(cell_nLocBas);
  for(int ii=0; ii<num_ebc; ++ii)
  {
    std::cout<<"-- ebc_id = "<<ii<<std::endl;
    std::cout<<"   local_global_node : \n";
    VEC_T::print( local_global_node[ii] );
    std::cout<<"   local_pt_xyz : \n";
    VEC_T::print( local_pt_xyz[ii] );
    std::cout<<"   local_global_cell : \n";
    VEC_T::print( local_global_cell[ii] );
    std::cout<<"   local_tri_ien : \n";
    VEC_T::print( local_tri_ien[ii] );
  } 
  std::cout<<"=========================================== \n";
}


// EOF
