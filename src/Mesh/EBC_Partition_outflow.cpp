#include "EBC_Partition_outflow.hpp"
#include "HDF5_Group.hpp"

EBC_Partition_outflow::EBC_Partition_outflow( 
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const ElemBC * const &ebc,
    const std::vector<INodalBC *> &nbc_list )
: EBC_Partition(part, mnindex, ebc)
{
  face_int_NA.resize(num_ebc);
  LID_all_face_nodes.resize(num_ebc);
  outvec.resize(num_ebc);

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      face_int_NA[ii] = ebc -> get_intNA(ii);
    
      // Obtain the old indices of the face nodes 
      const std::vector<int> old_node_idx = ebc -> get_global_node(ii);
      
      LID_all_face_nodes[ii].resize( old_node_idx.size() * 3 );
      for(unsigned int jj=0; jj<old_node_idx.size(); ++jj)
      {
        int idx = nbc_list[1] -> get_ID( old_node_idx[jj] );
        int idy = nbc_list[2] -> get_ID( old_node_idx[jj] );
        int idz = nbc_list[3] -> get_ID( old_node_idx[jj] );
        
        // If the id is not a dirichlet node, obtain the new numbering
        if(idx != -1) idx = mnindex -> get_old2new(idx);
        if(idy != -1) idy = mnindex -> get_old2new(idy);
        if(idz != -1) idz = mnindex -> get_old2new(idz);

        LID_all_face_nodes[ii][3*jj+0] = idx;
        LID_all_face_nodes[ii][3*jj+1] = idy;
        LID_all_face_nodes[ii][3*jj+2] = idz;
      }

      // Outward normal vector
      outvec[ii] = ebc -> get_normal_vec(ii);
    }
    else
    {
      face_int_NA[ii].clear();
      LID_all_face_nodes[ii].clear();
      outvec[ii] = Vector_3( 0.0, 0.0, 0.0 );
    }
  }
}

void EBC_Partition_outflow::write_hdf5( const std::string &FileName, 
    const std::string &GroupName ) const 
{
  // --------------------------------------------------------------------------
  // Call the base class writer to write the base class data
  EBC_Partition::write_hdf5( FileName, GroupName );
  // --------------------------------------------------------------------------

  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  auto h5w = SYS_T::make_unique<HDF5_Writer>( fName, H5F_ACC_RDWR );
  const hid_t file_id = h5w->get_file_id();
  auto root_group = HDF5_Group::open( file_id, GroupName );

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      std::string subgroup_name( "ebcid_" );
      subgroup_name.append( std::to_string(ii) );
      auto sub_group = HDF5_Group::open( root_group.id(), subgroup_name );

      h5w->write_doubleVector( sub_group.id(), "intNA", face_int_NA[ii] );
      h5w->write_intVector( sub_group.id(), "LID_all_face_nodes", LID_all_face_nodes[ii] );
      h5w->write_Vector_3( sub_group.id(), "out_normal", outvec[ii].to_std_array() );

    }
  }

}

// EOF
