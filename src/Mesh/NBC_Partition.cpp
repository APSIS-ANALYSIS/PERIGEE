#include "NBC_Partition.hpp"

NBC_Partition::NBC_Partition( const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const std::vector<INodalBC *> &nbc_list ) : cpu_rank(part->get_cpu_rank())
{
  const int dof = (int) nbc_list.size();

  LID.clear(); LDN.clear(); LPSN.clear(); LPMN.clear();
  LocalMaster.clear(); LocalMasterSlave.clear();

  Num_LD.resize(dof); Num_LPS.resize(dof); Num_LPM.resize(dof);

  for(int ii=0; ii<dof; ++ii)
  {
    Num_LD[ii] = 0;
    
    for(unsigned int jj=0; jj<nbc_list[ii]->get_num_dir_nodes(); ++jj)
    {
      unsigned int node_index = nbc_list[ii]->get_dir_nodes(jj);
      node_index = mnindex->get_old2new(node_index);

      if(part->isNodeInPart(node_index))
      {
        LDN.push_back(node_index);
        Num_LD[ii] += 1;
      }
    } // end jj-loop

    Num_LPS[ii] = 0; Num_LPM[ii] = 0;
    PERIGEE_OMP_PARALLEL_FOR
    for(unsigned int jj=0; jj<nbc_list[ii]->get_num_per_nodes(); ++jj)
    {
      unsigned int node_ps = nbc_list[ii]->get_per_slave_nodes(jj);
      unsigned int node_pm = nbc_list[ii]->get_per_master_nodes(jj);

      node_ps = mnindex->get_old2new(node_ps);
      node_pm = mnindex->get_old2new(node_pm);

      if(part->isNodeInPart(node_ps))
      {
        LPSN.push_back(node_ps);
        LPMN.push_back(node_pm);
        Num_LPS[ii] += 1;
      }

      if(part->isNodeInPart(node_pm))
      {
        LocalMaster.push_back(node_pm);
        LocalMasterSlave.push_back(node_ps);
        Num_LPM[ii] += 1;
      }
    } // end jj-loop
  } // end ii-loop over dof

  VEC_T::shrink2fit( LDN ); VEC_T::shrink2fit( LPSN ); VEC_T::shrink2fit( LPMN );
  VEC_T::shrink2fit( LocalMaster ); VEC_T::shrink2fit( LocalMasterSlave );

  const int totnode = part->get_nlocghonode();

  LID.resize(totnode * dof);
  PERIGEE_OMP_PARALLEL_FOR
  for(int ii=0; ii<dof; ++ii)
  {
    for(int jj=0; jj<totnode; ++jj)
    {
      const int new_index = part->get_local_to_global(jj);
      const int old_index = mnindex->get_new2old(new_index);
      LID[ii*totnode + jj] = nbc_list[ii]->get_ID(old_index);
    }
  }
  PERIGEE_OMP_PARALLEL_FOR
  for(int ii=0; ii<dof*totnode; ++ii)
  {
    if(LID[ii] != -1) LID[ii] = mnindex->get_old2new(LID[ii]);
  }

  VEC_T::shrink2fit( LID );
}

void NBC_Partition::write_hdf5( const std::string &FileName, 
    const std::string &GroupName ) const
{
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, GroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5writer = new HDF5_Writer(file_id);

  h5writer->write_intVector( g_id, "LID", LID );

  if( LDN.size() > 0 ) h5writer->write_intVector( g_id, "LDN", LDN );

  if( LPSN.size() > 0 )
  {
    h5writer->write_intVector( g_id, "LPSN", LPSN );
    h5writer->write_intVector( g_id, "LPMN", LPMN );
  }

  if( LocalMaster.size() > 0 )
  {
    h5writer->write_intVector( g_id, "LocalMaster",      LocalMaster );
    h5writer->write_intVector( g_id, "LocalMasterSlave", LocalMasterSlave );
  }

  h5writer->write_intVector(g_id, "Num_LD",  Num_LD);
  h5writer->write_intVector(g_id, "Num_LPS", Num_LPS);
  h5writer->write_intVector(g_id, "Num_LPM", Num_LPM);

  delete h5writer; H5Gclose(g_id); H5Fclose(file_id);
}

void NBC_Partition::print_info() const
{
  std::cout<<"=========================================== \n";
  std::cout<<"NBC_Partition : \n";
  std::cout<<"--- LID : \n";
  VEC_T::print(LID);
  std::cout<<"\n   LDN : \n";
  VEC_T::print(LDN);
  std::cout<<"\n   LPSN : \n";
  VEC_T::print(LPSN);
  std::cout<<"\n   LPMN : \n";
  VEC_T::print(LPMN);
  std::cout<<"\n   LocalMaster : \n";
  VEC_T::print(LocalMaster);
  std::cout<<"\n   LocalMasterSlave : \n";
  VEC_T::print(LocalMasterSlave);
  std::cout<<"\n   Num_LD : \n";
  VEC_T::print(Num_LD);
  std::cout<<"\n   Num_LPS : \n";
  VEC_T::print(Num_LPS);
  std::cout<<"\n   Num_LPM : \n";
  VEC_T::print(Num_LPM);
  std::cout<<"=========================================== \n";
}

// EOF
