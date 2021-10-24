#include "ALocal_NodalBC.hpp"

ALocal_NodalBC::ALocal_NodalBC( const std::string &fileBaseName, 
    const int &cpu_rank, const std::string &gname )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  num_nbc = h5r -> read_intScalar(gname.c_str(), "num_nbc");

  nlocghonode = h5r->read_intScalar( "Local_Node", "nlocghonode" );

  LID = h5r->read_intVector( gname.c_str(), "LID" );

  SYS_T::print_fatal_if( LID.size() % nlocghonode != 0, "Error:ALocal_NodalBC, LID length is not compatible with local and ghost node number. \n");

  dof = LID.size() / nlocghonode;

  std::string groupbase(gname);
  groupbase.append("/nbcid_");

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( SYS_T::to_string(nbc_id) );

    // Read local dirichlet nodes
    Num_LD[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "Num_LD" );

    int size_ldn = 0;
    for( unsigned int ii=0; ii<Num_LD[nbc_id].size(); ++ii )
      size_ldn += Num_LD[nbc_id][ii];

    LDN[nbc_id].clear(); LPSN[nbc_id].clear(); LPMN[nbc_id].clear();
    LocalMaster[nbc_id].clear(); LocalMasterSlave[nbc_id].clear();

    if( size_ldn > 0 )
      LDN[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "LDN" );

    SYS_T::print_fatal_if( int(LDN[nbc_id].size()) != size_ldn, "Error:ALocal_NodalBC, LDN length does not match Num_LD. \n" );
    
    // Read local periodic nodes
    Num_LPS[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "Num_LPS" );

    int size_lp = 0;
    for(unsigned int ii=0; ii<Num_LPS[nbc_id].size(); ++ii)
      size_lp += Num_LPS[nbc_id][ii];

    if( size_lp > 0 )
    {
      LPSN[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "LPSN" );
      LPMN[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "LPMN" );
    }

    SYS_T::print_fatal_if( int(LPSN[nbc_id].size()) != size_lp, "Error: ALocal_NodalBC, LPSN length does not match Num_LPS. \n" );
    SYS_T::print_fatal_if( int(LPMN[nbc_id].size()) != size_lp, "Error: ALocal_NodalBC, LPMN length does not match Num_LPS. \n" );

    // Read local periodic master nodes
    Num_LPM[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "Num_LPM" );

    size_lp = 0;
    for(unsigned int ii=0; ii<Num_LPM[nbc_id].size(); ++ii)
      size_lp += Num_LPM[nbc_id][ii];

    if( size_lp > 0 )
    {
      LocalMaster[nbc_id]      = h5r->read_intVector( subgroup_name.c_str(), "LocalMaster" );
      LocalMasterSlave[nbc_id] = h5r->read_intVector( subgroup_name.c_str(), "LocalMasterSlave" );
    }
    
    SYS_T::print_fatal_if( int(LocalMaster[nbc_id].size()) != size_lp, "Error: ALocal_NodalBC, LocalMaster length does not match Num_LPM. \n" );
    SYS_T::print_fatal_if( int(LocalMasterSlave[nbc_id].size()) != size_lp, "Error: ALocal_NodalBC, LocalMasterSlave length does not match Num_LPM. \n" );
  } // end nbc_id-loop

  delete h5r; H5Fclose( file_id );
  
  // Generate the offsets 
  LD_offset.resize(num_nbc); LPS_offset.resize(num_nbc); LPM_offset.resize(num_nbc);
  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    LD_offset[nbc_id].clear(); LPS_offset[nbc_id].clear(); LPM_offset[nbc_id].clear();

    LD_offset[nbc_id].push_back(0);
    LPS_offset[nbc_id].push_back(0);
    LPM_offset[nbc_id].push_back(0);

    for(int ii=1; ii<dof; ++ii)
    {
      LD_offset[nbc_id].push_back(  LD_offset[nbc_id][ii-1]  + Num_LD[nbc_id][ii-1]  );
      LPS_offset[nbc_id].push_back( LPS_offset[nbc_id][ii-1] + Num_LPS[nbc_id][ii-1] );
      LPM_offset[nbc_id].push_back( LPM_offset[nbc_id][ii-1] + Num_LPM[nbc_id][ii-1] );
    }

    VEC_T::shrink2fit( LD_offset[nbc_id]  );
    VEC_T::shrink2fit( LPS_offset[nbc_id] );
    VEC_T::shrink2fit( LPM_offset[nbc_id] );
  }
}

ALocal_NodalBC::~ALocal_NodalBC()
{
  clean_LocalMaster();
  VEC_T::clean(LID);
  VEC_T::clean(LDN); VEC_T::clean(LPSN); VEC_T::clean(LPMN);
  VEC_T::clean(Num_LD); VEC_T::clean(Num_LPS);
  VEC_T::clean(LD_offset); VEC_T::clean(LPS_offset);
}

void ALocal_NodalBC::print_info() const
{
  std::cout<<"ALocal_NodalBC: \n";

  std::cout<<"LID: \n";

  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof "<<ii<<'\n';
    for(int jj=0; jj<nlocghonode; ++jj)
      std::cout<<get_LID(ii,jj)<<'\t';
    std::cout<<"\n \n";
  }

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    std::cout<<"-- nbc_id = " << nbc_id << std::endl;
    std::cout<<"   LDN: \n";
    for(int ii=0; ii<dof; ++ii)
    {
      std::cout<<"   dof "<<ii<<'\n';
      for(int jj=0; jj<Num_LD[nbc_id][ii]; ++jj)
        std::cout<<"   "<<get_LDN(nbc_id, ii, jj)<<'\t';
      std::cout<<"\n \n";
    }

    for(int ii=0; ii<dof; ++ii)
    {
      std::cout<<"   dof: "<<ii<<'\t';
      std::cout<<"   Num_LD: " <<get_Num_LD( nbc_id, ii)<<'\t';
      std::cout<<"   Num_LPS: "<<get_Num_LPS(nbc_id, ii)<<'\t';
      std::cout<<"   Num_LPM: "<<get_Num_LPM(nbc_id, ii)<<'\n';
    }

    std::cout<<std::endl;

    std::cout<<"   LPSN: \n";
    for(int ii=0; ii<dof; ++ii)
    {
      std::cout<<"   dof "<<ii<<'\n';
      for(int jj=0; jj<Num_LPS[nbc_id][ii]; ++jj)
        std::cout<<"   "<<get_LPSN(nbc_id, ii, jj)<<'\t';
      std::cout<<"\n \n";
    }

    std::cout<<"   LPMN: \n";
    for(int ii=0; ii<dof; ++ii)
    {
      std::cout<<"   dof "<<ii<<'\n';
      for(int jj=0; jj<Num_LPS[nbc_id][ii]; ++jj)
        std::cout<<"   "<<get_LPMN(nbc_id, ii, jj)<<'\t';
      std::cout<<"\n \n";
    }

    std::cout<<"   LocalMaster: \n";
    for(int ii=0; ii<dof; ++ii)
    {
      std::cout<<"   dof "<<ii<<'\n';
      for(int jj=0; jj<Num_LPM[nbc_id][ii]; ++jj)
        std::cout<<get_LocalMaster(nbc_id, ii, jj)<<'\t';
      std::cout<<"\n \n";
    }

    std::cout<<"   LocalMasterSlave: \n";
    for(int ii=0; ii<dof; ++ii)
    {
      std::cout<<"   dof "<<ii<<'\n';
      for(int jj=0; jj<Num_LPM[nbc_id][ii]; ++jj)
        std::cout<<"   "<<get_LocalMasterSlave(nbc_id, ii, jj)<<'\t';
      std::cout<<"\n \n";
    }
  }
}

// EOF
