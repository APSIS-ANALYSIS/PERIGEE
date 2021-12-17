#include "NBC_Partition_MF.hpp"

NBC_Partition_MF::NBC_Partition_MF( const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const std::vector<INodalBC *> &nbc_list,
    const std::vector< std::vector<int> > &grid2id ) 
: NBC_Partition(part, mnindex, nbc_list)
{ 
  const int dof = (int) nbc_list.size();
  const int totnode = part->get_nlocghonode(); 

  SYS_T::print_fatal_if( VEC_T::get_size(grid2id) != dof, "Error: NBC_Partition_MF, the grid2id array size should math that of nbc list.\n" );
 
  // Generate an offset for accessing different dof's node 
  std::vector<int> LD_offset {0}, LPS_offset {0}, LPM_offset {0};

  for(int ii=1; ii<dof; ++ii)
  {
    LD_offset.push_back(  LD_offset[ii-1]  + Num_LD[ii-1]  );
    LPS_offset.push_back( LPS_offset[ii-1] + Num_LPS[ii-1] );
    LPM_offset.push_back( LPM_offset[ii-1] + Num_LPM[ii-1] );
  }

  LDN_MF.resize(LDN.size());
  LPSN_MF.resize(LPSN.size());
  LPMN_MF.resize(LPMN.size());
  LocalMaster_MF.resize(LocalMaster.size());
  LocalMasterSlave_MF.resize(LocalMasterSlave.size());
  LID_MF.resize(LID.size());

  for(int ii=0; ii<dof; ++ii)
  { 
    for(int jj=0; jj<Num_LD[ii]; ++jj)
    {
      const int loc = LD_offset[ii] + jj;
      LDN_MF[loc] = grid2id[ii][ LDN[loc] ];
    }

    for(int jj=0; jj<Num_LPS[ii]; ++jj)
    {
      const int loc = LPS_offset[ii] + jj;
      LPSN_MF[loc] = grid2id[ii][ LPSN[ loc ] ];
      LPMN_MF[loc] = grid2id[ii][ LPMN[ loc ] ];
    }

    for(int jj=0; jj<Num_LPM[ii]; ++jj)
    {
      const int loc = LPM_offset[ii] + jj;
      LocalMaster_MF[loc] = grid2id[ii][ LocalMaster[ loc ] ];
      LocalMasterSlave_MF[loc] = grid2id[ii][ LocalMasterSlave[ loc ] ];
    }

    for(int jj=0; jj<totnode; ++jj)
    {
      const int loc = ii * totnode + jj;
      if(LID[loc] != -1) LID_MF[loc] = grid2id[ii][ LID[loc] ];
    }
  } // end ii-loop over dof
  
  VEC_T::shrink2fit( LDN_MF ); 
  VEC_T::shrink2fit( LPSN_MF ); 
  VEC_T::shrink2fit( LPMN_MF ); 
  VEC_T::shrink2fit( LocalMaster_MF ); 
  VEC_T::shrink2fit( LocalMasterSlave_MF ); 
  VEC_T::shrink2fit( LID_MF );
}

NBC_Partition_MF::~NBC_Partition_MF()
{
  VEC_T::clean( LDN_MF ); 
  VEC_T::clean( LPSN_MF ); 
  VEC_T::clean( LPMN_MF ); 
  VEC_T::clean( LocalMaster_MF ); 
  VEC_T::clean( LocalMasterSlave_MF ); 
  VEC_T::clean( LID_MF );
}

void NBC_Partition_MF::write_hdf5( const std::string &FileName,
    const std::string &GroupName ) const
{
  // --------------------------------------------------------------------------
  // Call the base class writer to write the base class data
  NBC_Partition::write_hdf5( FileName, GroupName );
  // --------------------------------------------------------------------------
  
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  HDF5_Writer * h5writer = new HDF5_Writer(file_id);

  hid_t g_id = H5Gopen(file_id, GroupName.c_str(), H5P_DEFAULT);

  h5writer->write_intVector( g_id, "LID_MF", LID_MF );

  if( LDN_MF.size() > 0 ) h5writer->write_intVector( g_id, "LDN_MF", LDN_MF );

  if( LPSN_MF.size() > 0 )
  {
    h5writer->write_intVector( g_id, "LPSN_MF", LPSN_MF );
    h5writer->write_intVector( g_id, "LPMN_MF", LPMN_MF );
  }

  if( LocalMaster_MF.size() > 0 )
  {
    h5writer->write_intVector( g_id, "LocalMaster_MF",      LocalMaster_MF );
    h5writer->write_intVector( g_id, "LocalMasterSlave_MF", LocalMasterSlave_MF );
  }

  delete h5writer; H5Gclose(g_id); H5Fclose(file_id);
}

// EOF
