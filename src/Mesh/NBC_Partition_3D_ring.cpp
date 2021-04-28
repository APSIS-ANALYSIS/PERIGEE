#include "NBC_Partition_3D_ring.hpp"

NBC_Partition_3D_ring::NBC_Partition_3D_ring(
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const INodalBC * const &nbc ) 
: NBC_Partition_3D( part, mnindex, nbc ), 
  num_caps( nbc -> get_para_3() )
{
  local_cap_id.clear();
  local_pt_xyz.clear();

  if( LDN.size() > 0 )
  {
    // Access all (unpartitioned) ring node's cap_ids
    std::vector<int> cap_id; 
    nbc -> get_cap_id( cap_id );

    for(unsigned int ii=0; ii<nbc->get_num_dir_nodes(); ++ii)
    {
      unsigned int node_index = nbc -> get_dir_nodes(ii);
      node_index = mnindex -> get_old2new(node_index);

      if( part->isNodeInPart(node_index) )
      {
        local_cap_id.push_back( cap_id[ii] );
        local_pt_xyz.push_back( nbc -> get_pt_xyz(ii, 0) );
        local_pt_xyz.push_back( nbc -> get_pt_xyz(ii, 1) );
        local_pt_xyz.push_back( nbc -> get_pt_xyz(ii, 2) );
      }
    } 
  }

  nbc -> get_dominant_comp( dominant_comp );
  nbc -> get_outnormal( outnormal );
  nbc -> get_centroid( centroid );
}

NBC_Partition_3D_ring::~NBC_Partition_3D_ring()
{}

void NBC_Partition_3D_ring::write_hdf5( const char * FileName ) const
{
  std::string filebname(FileName);
  std::string fName = SYS_T::gen_partfile_name( filebname, cpu_rank );
  
  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t group_id = H5Gcreate(file_id, "/ring", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5writer = new HDF5_Writer(file_id);

  if(LDN.size() > 0)
  {
    h5writer->write_intVector( group_id, "LDN", LDN );
    h5writer->write_intVector( group_id, "local_cap_id", local_cap_id );
    h5writer->write_doubleVector( group_id, "local_pt_xyz", local_pt_xyz );
  }

  h5writer->write_intVector( group_id, "Num_LD", Num_LD );

  h5writer->write_intScalar( group_id, "num_caps", num_caps );

  h5writer->write_intVector( group_id, "cap_dominant_comp", dominant_comp );

  h5writer->write_doubleVector( group_id, "cap_out_normal", outnormal );

  h5writer->write_doubleVector( group_id, "cap_centroid", centroid );

  delete h5writer; H5Gclose(group_id); H5Fclose(file_id);
}

// EOF
