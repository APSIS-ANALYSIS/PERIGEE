#include "NBC_Partition_3D_inflow.hpp"

NBC_Partition_3D_inflow::NBC_Partition_3D_inflow(
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const INodalBC * const &nbc ) 
: NBC_Partition_3D( part, mnindex, nbc )
{
  actarea.resize(num_nbc); facearea.resize(num_nbc);
  outvec.resize(num_nbc);  centroid.resize(num_nbc);
  num_out_bc_pts.resize(num_nbc); outline_pts.resize(num_nbc);
  cell_nLocBas.resize(num_nbc);      local_tri_ien.resize(num_nbc);
  num_local_node.resize(num_nbc);    num_local_cell.resize(num_nbc);
  local_global_node.resize(num_nbc); local_global_cell.resize(num_nbc);
  local_node_pos.resize(num_nbc);    local_pt_xyz.resize(num_nbc);

  for(int ii=0; ii<num_nbc; ++ii)
  {
    // Area of the cap surface
    actarea[ii]  = nbc -> get_para_1(ii);
    facearea[ii] = nbc -> get_para_6(ii);

    // Outward normal vector
    outvec[ii] = nbc->get_para_2(ii);

    // Centroid and outline points
    num_out_bc_pts[ii] = nbc->get_para_3(ii);

    centroid[ii] = nbc->get_para_4(ii);

    outline_pts[ii].resize( 3*num_out_bc_pts[ii] );
    for(int jj=0; jj<3*num_out_bc_pts[ii]; ++jj)
      outline_pts[ii][jj] = nbc->get_para_5(ii, jj);

    // Record the geometrical info of the inlet in this CPU
    cell_nLocBas[ii] = nbc -> get_nLocBas(ii);

    std::vector<int> local_node, local_elem;
    local_node.clear(); local_elem.clear();

    local_global_cell[ii].clear();

    // Loop over all cells on inlet cap surface
    for(int jj=0; jj < nbc->get_num_cell(ii); ++jj)
    {
      const int elem_index = nbc -> get_global_cell(ii, jj); // cell vol id

      if( part -> get_elemLocIndex( elem_index ) != -1  )
      {
        local_elem.push_back( jj );

        local_global_cell[ii].push_back( elem_index );

        for( int kk=0; kk<cell_nLocBas[ii]; ++kk )
          local_node.push_back( nbc->get_ien(ii, jj, kk) );
      }
    }

    VEC_T::sort_unique_resize( local_node );

    num_local_node[ii] = static_cast<int>( local_node.size() );
    num_local_cell[ii] = static_cast<int>( local_global_cell[ii].size() );

    local_pt_xyz[ii].resize(      num_local_node[ii] * 3 );
    local_global_node[ii].resize( num_local_node[ii] );
    local_node_pos[ii].resize(    num_local_node[ii] );

    for(int jj=0; jj<num_local_node[ii]; ++jj)
    {
      local_pt_xyz[ii][3*jj+0] = nbc -> get_pt_xyz( ii, local_node[jj], 0 );
      local_pt_xyz[ii][3*jj+1] = nbc -> get_pt_xyz( ii, local_node[jj], 1 );
      local_pt_xyz[ii][3*jj+2] = nbc -> get_pt_xyz( ii, local_node[jj], 2 );

      local_global_node[ii][jj] = nbc -> get_global_node( ii, local_node[jj] );

      local_node_pos[ii][jj] = part->get_nodeLocGhoIndex( mnindex->get_old2new(local_global_node[ii][jj]) );
    }

    // create new IEN
    local_tri_ien[ii].resize( num_local_cell[ii] * cell_nLocBas[ii] );

    for(int jj=0; jj<num_local_cell[ii]; ++jj)
    {
      for(int kk=0; kk<cell_nLocBas[ii]; ++kk)
      {
        const int temp_node = nbc -> get_ien( ii, local_elem[jj], kk );
        const int temp_npos = VEC_T::get_pos( local_node, temp_node );
        local_tri_ien[ii][jj*cell_nLocBas[ii] + kk] = temp_npos;
      }
    }
  } // end ii-loop over num_nbc
}

NBC_Partition_3D_inflow::~NBC_Partition_3D_inflow()
{}

void NBC_Partition_3D_inflow::write_hdf5( const std::string &FileName ) const
{
  std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );
  
  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, "/inflow", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  h5w -> write_intScalar( g_id, "num_nbc", num_nbc );

  h5w->write_doubleVector( g_id, "Inflow_active_area", actarea );
  
  h5w->write_doubleVector( g_id, "Inflow_full_area", facearea );

  h5w->write_intVector( g_id, "num_out_bc_pts", num_out_bc_pts);

  h5w->write_intVector( g_id, "num_local_node", num_local_node );
  
  h5w->write_intVector( g_id, "num_local_cell", num_local_cell );
  
  h5w->write_intVector( g_id, "cell_nLocBas", cell_nLocBas );

  for(int ii=0; ii<num_nbc; ++ii)
  {
    std::string subgroup_name( "nbcid_" );
    subgroup_name.append( SYS_T::to_string(ii) );

    hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(),
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if( LDN[ii].size() > 0 )
      h5w->write_intVector( group_id, "LDN", LDN[ii] );

    h5w->write_intVector( group_id, "Num_LD", Num_LD[ii] );

    h5w->write_Vector_3( group_id, "Outward_normal_vector", outvec[ii] );

    h5w->write_Vector_3( group_id, "centroid", centroid[ii] );

    h5w->write_doubleVector( group_id, "outline_pts", outline_pts[ii] );

    h5w->write_doubleVector( group_id, "local_pt_xyz", local_pt_xyz[ii] );

    h5w->write_intVector( group_id, "local_tri_ien", local_tri_ien[ii] );

    h5w->write_intVector( group_id, "local_global_node", local_global_node[ii] );

    h5w->write_intVector( group_id, "local_node_pos", local_node_pos[ii] );

    h5w->write_intVector( group_id, "local_global_cell", local_global_cell[ii] );

    H5Gclose( group_id );
  }

  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}

// EOF
