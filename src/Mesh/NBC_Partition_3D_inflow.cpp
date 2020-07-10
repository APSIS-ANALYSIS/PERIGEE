#include "NBC_Partition_3D_inflow.hpp"

NBC_Partition_3D_inflow::NBC_Partition_3D_inflow(
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const INodalBC * const &nbc ) 
: NBC_Partition_3D( part, mnindex, nbc )
{
  actarea  = nbc -> get_para_1();
  facearea = nbc -> get_para_6();

  outvec.clear();
  outvec.push_back( nbc->get_para_2(0) );
  outvec.push_back( nbc->get_para_2(1) );
  outvec.push_back( nbc->get_para_2(2) );

  num_out_bc_pts = nbc->get_para_3();

  centroid.resize(3);
  centroid[0] = nbc->get_para_4(0);
  centroid[1] = nbc->get_para_4(1);
  centroid[2] = nbc->get_para_4(2);

  outline_pts.resize( 3*num_out_bc_pts );
  for(int ii=0; ii<3*num_out_bc_pts; ++ii)
    outline_pts[ii] = nbc->get_para_5(ii);

  // Record the geometrical info of the inlet in this CPU
  cell_nLocBas = nbc -> get_nLocBas();

  std::vector<int> local_node, local_elem;
  local_node.clear(); local_elem.clear();

  local_global_cell.clear();

  // Loop over all cells on inlet surface
  for(int jj=0; jj < nbc -> get_num_cell(); ++jj)
  {
    const int elem_index = nbc -> get_global_cell(jj); // cell vol id

    if( part -> get_elemLocIndex( elem_index ) != -1  )
    {
      local_elem.push_back( jj );

      local_global_cell.push_back( elem_index );

      for( int kk=0; kk < cell_nLocBas; ++kk )
        local_node.push_back( nbc->get_ien(jj, kk) );
    }
  }

  VEC_T::sort_unique_resize( local_node );

  num_local_node = static_cast<int>( local_node.size() );
  num_local_cell = static_cast<int>( local_global_cell.size() );

  local_pt_xyz.resize( num_local_node * 3 );
  local_global_node.resize( num_local_node );
  local_node_pos.resize( num_local_node );

  for(int jj=0; jj<num_local_node; ++jj)
  {
    local_pt_xyz[3*jj+0] = nbc -> get_pt_xyz( local_node[jj], 0 );
    local_pt_xyz[3*jj+1] = nbc -> get_pt_xyz( local_node[jj], 1 );
    local_pt_xyz[3*jj+2] = nbc -> get_pt_xyz( local_node[jj], 2 );

    local_global_node[jj] = nbc -> get_global_node( local_node[jj] );

    local_node_pos[jj] = part->get_nodeLocGhoIndex( mnindex->get_old2new(local_global_node[jj]) );
  }

  // create new IEN
  local_tri_ien.resize( num_local_cell * cell_nLocBas );

  for(int jj=0; jj<num_local_cell; ++jj)
  {
    for(int kk=0; kk<cell_nLocBas; ++kk)
    {
      const int temp_node = nbc -> get_ien( local_elem[jj], kk );
      const int temp_npos = VEC_T::get_pos( local_node, temp_node );
      local_tri_ien[jj*cell_nLocBas + kk] = temp_npos;
    }
  }
}


NBC_Partition_3D_inflow:: ~NBC_Partition_3D_inflow()
{}


void NBC_Partition_3D_inflow::write_hdf5( const char * FileName ) const
{
  std::string filebname(FileName);
  std::string fName = SYS_T::gen_partfile_name( filebname, cpu_rank );
  
  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t group_id = H5Gcreate(file_id, "/inflow", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5writer = new HDF5_Writer(file_id);

  if(LDN.size() > 0)
    h5writer->write_intVector(group_id, "LDN", LDN);

  h5writer->write_intVector(group_id, "Num_LD", Num_LD);

  h5writer->write_doubleVector( group_id, "Outward_normal_vector", outvec );

  h5writer->write_doubleScalar( group_id, "Inflow_active_area", actarea );
  
  h5writer->write_doubleScalar( group_id, "Inflow_full_area", facearea );

  h5writer->write_intScalar( group_id, "num_out_bc_pts", num_out_bc_pts);

  h5writer->write_doubleVector( group_id, "centroid", centroid);

  h5writer->write_doubleVector( group_id, "outline_pts", outline_pts);

  h5writer->write_intScalar( group_id, "num_local_node", num_local_node );
  
  h5writer->write_intScalar( group_id, "num_local_cell", num_local_cell );
  
  h5writer->write_intScalar( group_id, "cell_nLocBas", cell_nLocBas );

  h5writer->write_doubleVector( group_id, "local_pt_xyz", local_pt_xyz );

  h5writer->write_intVector( group_id, "local_tri_ien", local_tri_ien );

  h5writer->write_intVector( group_id, "local_global_node", local_global_node );

  h5writer->write_intVector( group_id, "local_node_pos", local_node_pos );

  h5writer->write_intVector( group_id, "local_global_cell", local_global_cell );

  delete h5writer; H5Gclose(group_id); H5Fclose(file_id);
}

// EOF
