#include "NBC_Partition_inflow.hpp"
#include "HDF5_Group.hpp"

NBC_Partition_inflow::NBC_Partition_inflow(
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const INodalBC * const &nbc ) 
: cpu_rank(part->get_cpu_rank()), num_nbc( nbc -> get_num_nbc() )
{
  LDN.resize(num_nbc);
  for(int ii=0; ii<num_nbc; ++ii) LDN[ii].clear();

  Num_LD.resize(num_nbc);

  actarea.resize(num_nbc); facearea.resize(num_nbc);
  outvec.resize(num_nbc);  centroid.resize(num_nbc);
  num_out_bc_pts.resize(num_nbc); outline_pts.resize(num_nbc);
  cell_nLocBas.resize(num_nbc);      local_cell_ien.resize(num_nbc);
  num_local_node.resize(num_nbc);    num_local_cell.resize(num_nbc);
  local_global_node.resize(num_nbc); local_global_cell.resize(num_nbc);
  local_node_pos.resize(num_nbc);    local_pt_xyz.resize(num_nbc);

  for(int ii=0; ii<num_nbc; ++ii)
  {
    // Collect the Dirichlet nodes on the ii-th inlet surface
    Num_LD[ii] = 0;

    for(unsigned int jj=0; jj<nbc->get_num_dir_nodes_on_inlet(ii); ++jj)
    {
      int node_index = static_cast<int>( nbc -> get_dir_nodes_on_inlet(ii, jj) );
      node_index = mnindex -> get_old2new(node_index);
      if(part->isNodeInPart(node_index))
      {
        LDN[ii].push_back(node_index);
        Num_LD[ii] += 1;
      }
    }

    // Area of the cap surface
    actarea[ii]  = nbc -> get_inf_active_area(ii);
    facearea[ii] = nbc -> get_face_area(ii);

    // Outward normal vector
    outvec[ii] = nbc->get_outnormal(ii);

    // Centroid and outline points
    num_out_bc_pts[ii] = nbc->get_num_out_bc_pts(ii);

    centroid[ii] = nbc->get_centroid(ii);

    outline_pts[ii].resize( 3*num_out_bc_pts[ii] );
    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<3*num_out_bc_pts[ii]; ++jj)
      outline_pts[ii][jj] = nbc->get_outline_pts(ii, jj);

    // Record the geometrical info of the inlet in this CPU
    cell_nLocBas[ii] = nbc -> get_nLocBas(ii);

    std::vector<int> local_node {}, local_elem {};

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

    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_local_node[ii]; ++jj)
    {
      local_pt_xyz[ii][3*jj+0] = nbc -> get_pt_xyz( ii, local_node[jj], 0 );
      local_pt_xyz[ii][3*jj+1] = nbc -> get_pt_xyz( ii, local_node[jj], 1 );
      local_pt_xyz[ii][3*jj+2] = nbc -> get_pt_xyz( ii, local_node[jj], 2 );

      local_global_node[ii][jj] = nbc -> get_global_node( ii, local_node[jj] );

      local_node_pos[ii][jj] = part->get_nodeLocGhoIndex( mnindex->get_old2new(local_global_node[ii][jj]) );
    }

    // create new IEN
    local_cell_ien[ii].resize( num_local_cell[ii] * cell_nLocBas[ii] );

    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_local_cell[ii]; ++jj)
    {
      for(int kk=0; kk<cell_nLocBas[ii]; ++kk)
      {
        const int temp_node = nbc -> get_ien( ii, local_elem[jj], kk );
        const int temp_npos = VEC_T::get_pos( local_node, temp_node );
        local_cell_ien[ii][jj*cell_nLocBas[ii] + kk] = temp_npos;
      }
    }
  } // end ii-loop over num_nbc
}

void NBC_Partition_inflow::write_hdf5( const std::string &FileName ) const
{
  std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  auto h5w = SYS_T::make_unique<HDF5_Writer>(fName, H5F_ACC_RDWR);
  const hid_t file_id = h5w->get_file_id();

  auto inflow_group = HDF5_Group::create( file_id, "/inflow" );

  h5w -> write_intScalar( inflow_group.id(), "num_nbc", num_nbc );

  for(int ii=0; ii<num_nbc; ++ii)
  {
    const std::string subgroup_name = std::string("nbcid_") + std::to_string(ii);

    auto sub_group = HDF5_Group::create( inflow_group.id(), subgroup_name );

    h5w->write_intScalar( sub_group.id(), "Num_LD", Num_LD[ii] );
    
    h5w->write_intVector( sub_group.id(), "LDN", LDN[ii] );

    h5w->write_doubleScalar( sub_group.id(), "Inflow_active_area", actarea[ii] );

    h5w->write_doubleScalar( sub_group.id(), "Inflow_full_area", facearea[ii] );

    h5w->write_intScalar( sub_group.id(), "num_out_bc_pts", num_out_bc_pts[ii] );

    h5w->write_intScalar( sub_group.id(), "num_local_node", num_local_node[ii] );

    h5w->write_intScalar( sub_group.id(), "num_local_cell", num_local_cell[ii] );

    h5w->write_intScalar( sub_group.id(), "cell_nLocBas", cell_nLocBas[ii] );

    h5w->write_Vector_3( sub_group.id(), "Outward_normal_vector", outvec[ii].to_std_array() );

    h5w->write_Vector_3( sub_group.id(), "centroid", centroid[ii].to_std_array() );

    h5w->write_doubleVector( sub_group.id(), "outline_pts", outline_pts[ii] );

    h5w->write_doubleVector( sub_group.id(), "local_pt_xyz", local_pt_xyz[ii] );

    h5w->write_intVector( sub_group.id(), "local_cell_ien", local_cell_ien[ii] );

    h5w->write_intVector( sub_group.id(), "local_global_node", local_global_node[ii] );

    h5w->write_intVector( sub_group.id(), "local_node_pos", local_node_pos[ii] );

    h5w->write_intVector( sub_group.id(), "local_global_cell", local_global_cell[ii] );
  }
}

// EOF
