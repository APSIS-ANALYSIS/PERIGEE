#include "NBC_Partition_3D.hpp"

NBC_Partition_3D::NBC_Partition_3D( const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const std::vector<INodalBC *> &nbc_list )
: cpu_rank(part->get_cpu_rank()), num_nbc( nbc_list[0]->get_num_nbc() )
{
  const int dof = (int) nbc_list.size();
  
  for(int ii=1; ii<dof; ++ii)
    SYS_T::print_exit_if( num_nbc != nbc_list[ii]->get_num_nbc(), "Error: NBC_Partition_3D assumes all input INodalBC have the same num_nbc.\n");

  LID.clear();
  
  LDN.resize(num_nbc);         Num_LD.resize(num_nbc);
  LPSN.resize(num_nbc);        LPMN.resize(num_nbc);
  LocalMaster.resize(num_nbc); LocalMasterSlave.resize(num_nbc);
  Num_LPS.resize(num_nbc);     Num_LPM.resize(num_nbc);

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    LDN[nbc_id].clear(); LPSN[nbc_id].clear(); LPMN[nbc_id].clear();
    LocalMaster[nbc_id].clear(); LocalMasterSlave[nbc_id].clear();

    Num_LD[nbc_id].resize(dof);
    Num_LPS[nbc_id].resize(dof); Num_LPM[nbc_id].resize(dof);

    for(int ii=0; ii<dof; ++ii)
    {
      unsigned int node_num = 0;
      unsigned int ps_num = 0;
      unsigned int pm_num = 0;

      for(unsigned int jj=0; jj<nbc_list[ii]->get_num_dir_nodes(nbc_id); ++jj)
      {
        unsigned int node_index = nbc_list[ii]->get_dir_nodes(nbc_id, jj);
        node_index = mnindex->get_old2new(node_index);

        if(part->isNodeInPart(node_index))
        {
          LDN[nbc_id].push_back(node_index);
          node_num += 1;
        }
      } // end jj-loop

      Num_LD[nbc_id][ii] = node_num;

      for(unsigned int jj=0; jj<nbc_list[ii]->get_num_per_nodes(nbc_id); ++jj)
      {
        unsigned int node_ps = nbc_list[ii]->get_per_slave_nodes( nbc_id, jj);
        unsigned int node_pm = nbc_list[ii]->get_per_master_nodes(nbc_id, jj);

        node_ps = mnindex->get_old2new(node_ps);
        node_pm = mnindex->get_old2new(node_pm);
        
        if(part->isNodeInPart(node_ps))
        {
          LPSN[ii].push_back(node_ps);
          LPMN[ii].push_back(node_pm);
          ps_num += 1;
        }

        if(part->isNodeInPart(node_pm))
        {
          LocalMaster[ii].push_back(node_pm);
          LocalMasterSlave[ii].push_back(node_ps);
          pm_num += 1;
        }
      } // end jj-loop

      Num_LPS[nbc_id][ii] = ps_num;
      Num_LPM[nbc_id][ii] = pm_num;

    } // end ii-loop over dof

    VEC_T::shrink2fit( LDN[nbc_id] );
    VEC_T::shrink2fit( LPSN[nbc_id] ); VEC_T::shrink2fit( LPMN[nbc_id] );
    VEC_T::shrink2fit( LocalMaster[nbc_id] ); VEC_T::shrink2fit( LocalMasterSlave[nbc_id] );

  } // end nbc_id-loop

  const int totnode = part->get_nlocghonode();

  LID.resize(totnode * dof);

  for(int ii=0; ii<dof; ++ii)
  {
    for(int jj=0; jj<totnode; ++jj)
    {
      const int new_index = part->get_local_to_global(jj);
      const int old_index = mnindex->get_new2old(new_index);
      LID[ii*totnode + jj] = nbc_list[ii]->get_ID(old_index);
    }
  }

  for(int ii=0; ii<dof*totnode; ++ii)
  {
    if(LID[ii] != -1)
    {
      LID[ii] = mnindex->get_old2new(LID[ii]);
    }
  }

  VEC_T::shrink2fit( LID );
}


NBC_Partition_3D::NBC_Partition_3D( const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const INodalBC * const &nbc )
: cpu_rank(part->get_cpu_rank()), num_nbc( nbc->get_num_nbc() )
{
  const int dof = 1;

  LID.clear();
  
  LDN.resize(num_nbc);         Num_LD.resize(num_nbc);
  LPSN.resize(num_nbc);        LPMN.resize(num_nbc);
  LocalMaster.resize(num_nbc); LocalMasterSlave.resize(num_nbc);
  Num_LPS.resize(num_nbc);     Num_LPM.resize(num_nbc);

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    LDN[nbc_id].clear(); LPSN[nbc_id].clear(); LPMN[nbc_id].clear();
    LocalMaster[nbc_id].clear(); LocalMasterSlave[nbc_id].clear();

    Num_LD[nbc_id].resize(dof);
    Num_LPS[nbc_id].resize(dof); Num_LPM[nbc_id].resize(dof);

    unsigned int node_num = 0;
    unsigned int ps_num = 0;
    unsigned int pm_num = 0;

    for(unsigned int jj=0; jj<nbc->get_num_dir_nodes(nbc_id); ++jj)
    {
      unsigned int node_index = nbc -> get_dir_nodes(nbc_id, jj);
      node_index = mnindex -> get_old2new(node_index);
      if(part->isNodeInPart(node_index))
      {
        LDN[nbc_id].push_back(node_index);
        node_num += 1;
      }
    } // end jj-loop

    Num_LD[nbc_id][0] = node_num;

    for(unsigned int jj=0; jj<nbc->get_num_per_nodes(nbc_id); ++jj)
    {
      unsigned int node_ps = nbc -> get_per_slave_nodes( nbc_id, jj);
      unsigned int node_pm = nbc -> get_per_master_nodes(nbc_id, jj);

      node_ps = mnindex -> get_old2new(node_ps);
      node_pm = mnindex -> get_old2new(node_pm);

      if(part->isNodeInPart(node_ps))
      {
        LPSN[nbc_id].push_back(node_ps);
        LPMN[nbc_id].push_back(node_pm);
        ps_num += 1;
      }

      if(part->isNodeInPart(node_pm))
      {
        LocalMaster[nbc_id].push_back(node_pm);
        LocalMasterSlave[nbc_id].push_back(node_ps);
        pm_num += 1;
      }
    } // end jj-loop

    Num_LPS[nbc_id][0] = ps_num;
    Num_LPM[nbc_id][0] = pm_num;

    VEC_T::shrink2fit( LDN[nbc_id] );
    VEC_T::shrink2fit( LPSN[nbc_id] ); VEC_T::shrink2fit( LPMN[nbc_id] );
    VEC_T::shrink2fit( LocalMaster[nbc_id] ); VEC_T::shrink2fit( LocalMasterSlave[nbc_id] );

  } // end nbc_id-loop

  const int totnode = part->get_nlocghonode();

  LID.resize(totnode * dof);

  for(int jj=0; jj<totnode; ++jj)
  {
    const int new_index = part->get_local_to_global(jj);
    const int old_index = mnindex->get_new2old(new_index);
    LID[jj] = nbc -> get_ID(old_index);
  }

  for(int ii=0; ii<dof*totnode; ++ii)
  {
    if(LID[ii] != -1)
    {
      LID[ii] = mnindex->get_old2new(LID[ii]);
    }
  }

  VEC_T::shrink2fit( LID );
}


NBC_Partition_3D::~NBC_Partition_3D()
{
  VEC_T::clean(LID);
  VEC_T::clean(LDN);
  VEC_T::clean(LPSN);
  VEC_T::clean(LPMN);
  VEC_T::clean(LocalMaster);
  VEC_T::clean(LocalMasterSlave);
  VEC_T::clean(Num_LD);
  VEC_T::clean(Num_LPS);
  VEC_T::clean(Num_LPM);
}


void NBC_Partition_3D::write_hdf5(const char * FileName) const
{
  const std::string input_fName(FileName);
  const std::string fName = SYS_T::gen_partfile_name( input_fName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, "/nbc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5writer = new HDF5_Writer(file_id);

  h5writer->write_intScalar( g_id, "num_nbc", num_nbc );

  h5writer->write_intVector( g_id, "LID", LID );

  const std::string groupbase("nbcid_");

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    std::string subgroup_name(groupbase);
      subgroup_name.append( SYS_T::to_string(nbc_id) );

    hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(),
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if( LDN[nbc_id].size() > 0 )
      h5writer->write_intVector( group_id, "LDN", LDN[nbc_id] );

    if( LPSN[nbc_id].size() > 0)
    {
      h5writer->write_intVector( group_id, "LPSN", LPSN[nbc_id] );
      h5writer->write_intVector( group_id, "LPMN", LPMN[nbc_id] );
    }

    if( LocalMaster[nbc_id].size() > 0 )
    {
      h5writer->write_intVector( group_id, "LocalMaster",      LocalMaster[nbc_id] );
      h5writer->write_intVector( group_id, "LocalMasterSlave", LocalMasterSlave[nbc_id] );
    }

    h5writer->write_intVector(group_id, "Num_LD",  Num_LD[nbc_id]);
    h5writer->write_intVector(group_id, "Num_LPS", Num_LPS[nbc_id]);
    h5writer->write_intVector(group_id, "Num_LPM", Num_LPM[nbc_id]);

    H5Gclose( group_id );
  }

  delete h5writer; H5Gclose(g_id); H5Fclose(file_id);
}


void NBC_Partition_3D::write_hdf5( const char * FileName,
    const char * GroupName ) const
{
  const std::string input_fName(FileName);
  const std::string fName = SYS_T::gen_partfile_name( input_fName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, GroupName, H5P_DEFAULT, H5P_DEFAULT, 
      H5P_DEFAULT);

  HDF5_Writer * h5writer = new HDF5_Writer(file_id);

  h5writer->write_intScalar( g_id, "num_nbc", num_nbc );

  h5writer->write_intVector( g_id, "LID", LID );

  const std::string groupbase("nbcid_");

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    std::string subgroup_name(groupbase);
      subgroup_name.append( SYS_T::to_string(nbc_id) );

    hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(),
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if( LDN[nbc_id].size() > 0 )
      h5writer->write_intVector( group_id, "LDN", LDN[nbc_id] );

    if( LPSN[nbc_id].size() > 0)
    {
      h5writer->write_intVector( group_id, "LPSN", LPSN[nbc_id] );
      h5writer->write_intVector( group_id, "LPMN", LPMN[nbc_id] );
    }

    if( LocalMaster[nbc_id].size() > 0 )
    {
      h5writer->write_intVector( group_id, "LocalMaster",      LocalMaster[nbc_id] );
      h5writer->write_intVector( group_id, "LocalMasterSlave", LocalMasterSlave[nbc_id] );
    }

    h5writer->write_intVector(group_id, "Num_LD",  Num_LD[nbc_id]);
    h5writer->write_intVector(group_id, "Num_LPS", Num_LPS[nbc_id]);
    h5writer->write_intVector(group_id, "Num_LPM", Num_LPM[nbc_id]);

    H5Gclose( group_id );
  }

  delete h5writer; H5Gclose(g_id); H5Fclose(file_id);
}


void NBC_Partition_3D::print_info() const
{
  std::cout<<"=========================================== \n";
  std::cout<<"NBC_Partition_3D : \n";
  std::cout<<"-- num_nbc = "<<num_nbc<<std::endl;
  std::cout<<"--- LID : \n";
  VEC_T::print(LID);

  for(int ii=0; ii<num_nbc; ++ii)
  {
    std::cout<<"-- nbc_id = "<<ii<<std::endl;
    std::cout<<"\n   LDN : \n";
    VEC_T::print(LDN[ii]);
    std::cout<<"\n   LPSN : \n";
    VEC_T::print(LPSN[ii]);
    std::cout<<"\n   LPMN : \n";
    VEC_T::print(LPMN[ii]);
    std::cout<<"\n   LocalMaster : \n";
    VEC_T::print(LocalMaster[ii]);
    std::cout<<"\n   LocalMasterSlave : \n";
    VEC_T::print(LocalMasterSlave[ii]);
    std::cout<<"\n   Num_LD : \n";
    VEC_T::print(Num_LD[ii]);
    std::cout<<"\n   Num_LPS : \n";
    VEC_T::print(Num_LPS[ii]);
    std::cout<<"\n   Num_LPM : \n";
    VEC_T::print(Num_LPM[ii]);
  }
  std::cout<<"=========================================== \n";
}

// EOF
