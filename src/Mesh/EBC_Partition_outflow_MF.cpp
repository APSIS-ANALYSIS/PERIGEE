#include "EBC_Partition_outflow_MF.hpp"

EBC_Partition_outflow_MF::EBC_Partition_outflow_MF(
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const ElemBC * const &ebc,
    const std::vector<INodalBC *> &nbc_list,
    const std::vector< std::vector<int> > &grid2id )
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
        int idx = nbc_list[0] -> get_ID( old_node_idx[jj] );
        int idy = nbc_list[1] -> get_ID( old_node_idx[jj] );
        int idz = nbc_list[2] -> get_ID( old_node_idx[jj] );

        // If the id is not a dirichlet node, obtain the new numbering
        if(idx != -1)
        {
          idx = mnindex -> get_old2new(idx);
          LID_all_face_nodes[ii][3*jj+0] = grid2id[0][idx];
        }
        else
          LID_all_face_nodes[ii][3*jj+0] = -1;

        if(idy != -1)
        {
          idy = mnindex -> get_old2new(idy);
          LID_all_face_nodes[ii][3*jj+1] = grid2id[1][idy];
        }
        else
          LID_all_face_nodes[ii][3*jj+1] = -1;

        if(idz != -1)
        {
          idz = mnindex -> get_old2new(idz);
          LID_all_face_nodes[ii][3*jj+2] = grid2id[2][idz];
        }
        else
          LID_all_face_nodes[ii][3*jj+2] = -1;
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

EBC_Partition_outflow_MF::~EBC_Partition_outflow_MF()
{
  for(int ii=0; ii<num_ebc; ++ii)
  {
    VEC_T::clean( face_int_NA[ii] );
    VEC_T::clean( LID_all_face_nodes[ii] );
  }
  VEC_T::clean( face_int_NA );
  VEC_T::clean( LID_all_face_nodes );
  VEC_T::clean( outvec );
}

void EBC_Partition_outflow_MF::write_hdf5( const std::string &FileName,
    const std::string &GroupName ) const
{
  // --------------------------------------------------------------------------
  // Call the base class writer to write the base class data
  EBC_Partition::write_hdf5( FileName, GroupName );
  // --------------------------------------------------------------------------

  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t g_id = H5Gopen( file_id, GroupName.c_str(), H5P_DEFAULT );

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  for(int ii=0; ii<num_ebc; ++ii)
  {
    if( num_local_cell[ii] > 0 )
    {
      std::string subgroup_name( "ebcid_" );
      subgroup_name.append( std::to_string(ii) );
      hid_t subgroup_id = H5Gopen(g_id, subgroup_name.c_str(), H5P_DEFAULT );

      h5w->write_doubleVector( subgroup_id, "intNA", face_int_NA[ii] );
      h5w->write_intVector( subgroup_id, "LID_all_face_nodes", LID_all_face_nodes[ii] );
      h5w->write_Vector_3( subgroup_id, "out_normal", outvec[ii] );

      H5Gclose( subgroup_id );
    }
  }

  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}

// EOF
