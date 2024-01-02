#include "ALocal_NBC.hpp"

ALocal_NBC::ALocal_NBC( const std::string &fileBaseName, 
    const int &cpu_rank, const std::string &gname )
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  nlocghonode = h5r->read_intScalar( "Local_Node", "nlocghonode" );
  LID = h5r->read_intVector( gname.c_str(), "LID" );

  SYS_T::print_fatal_if( LID.size() % nlocghonode != 0, "Error:ALocal_NBC, LID length is not compatible with local and ghost node number. \n");

  dof = LID.size() / nlocghonode;

  // Read local dirichlet nodes
  Num_LD = h5r->read_intVector( gname.c_str(), "Num_LD" );

  LDN.clear(); LPSN.clear(); LPMN.clear();
  LocalMaster.clear(); LocalMasterSlave.clear();

  if( VEC_T::sum( Num_LD ) > 0 ) LDN = h5r->read_intVector( gname.c_str(), "LDN" );

  SYS_T::print_fatal_if( int(LDN.size()) != VEC_T::sum( Num_LD ), "Error:ALocal_NBC, LDN length does not match Num_LD. \n" );

  // Read local periodic nodes
  Num_LPS = h5r->read_intVector( gname.c_str(), "Num_LPS" );

  if( VEC_T::sum( Num_LPS ) > 0 )
  {
    LPSN = h5r->read_intVector( gname.c_str(), "LPSN" );
    LPMN = h5r->read_intVector( gname.c_str(), "LPMN" );
  }

  SYS_T::print_fatal_if( int(LPSN.size()) != VEC_T::sum( Num_LPS ), "Error: ALocal_NBC, LPSN length does not match Num_LPS. \n" );
  SYS_T::print_fatal_if( int(LPMN.size()) != VEC_T::sum( Num_LPS ), "Error: ALocal_NBC, LPMN length does not match Num_LPS. \n" );

  // Read local periodic master nodes
  Num_LPM = h5r->read_intVector( gname.c_str(), "Num_LPM" );

  if( VEC_T::sum( Num_LPM ) > 0 )
  {
    LocalMaster      = h5r->read_intVector( gname.c_str(), "LocalMaster" );
    LocalMasterSlave = h5r->read_intVector( gname.c_str(), "LocalMasterSlave" );
  }

  SYS_T::print_fatal_if( int(LocalMaster.size()) != VEC_T::sum( Num_LPM ), "Error: ALocal_NBC, LocalMaster length does not match Num_LPM. \n" );
  SYS_T::print_fatal_if( int(LocalMasterSlave.size()) != VEC_T::sum( Num_LPM ), "Error: ALocal_NBC, LocalMasterSlave length does not match Num_LPM. \n" );

  delete h5r; H5Fclose( file_id );

  // Generate the offsets 
  LD_offset.clear(); LPS_offset.clear(); LPM_offset.clear();
  LD_offset.push_back(0); LPS_offset.push_back(0); LPM_offset.push_back(0);

  for(int ii=1; ii<dof; ++ii)
  {
    LD_offset.push_back(  LD_offset[ii-1]  + Num_LD[ii-1]  );
    LPS_offset.push_back( LPS_offset[ii-1] + Num_LPS[ii-1] );
    LPM_offset.push_back( LPM_offset[ii-1] + Num_LPM[ii-1] );
  }

  VEC_T::shrink2fit( LD_offset  );
  VEC_T::shrink2fit( LPS_offset );
  VEC_T::shrink2fit( LPM_offset );
}

void ALocal_NBC::print_info() const
{
  std::cout<<"ALocal_NBC: \n";

  std::cout<<"LID: \n";

  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof "<<ii<<'\n';
    for(int jj=0; jj<nlocghonode; ++jj)
      std::cout<<get_LID(ii,jj)<<'\t';
    std::cout<<"\n \n";
  }

  std::cout<<"LDN: \n";
  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof "<<ii<<'\n';
    for(int jj=0; jj<Num_LD[ii]; ++jj)
      std::cout<<get_LDN(ii, jj)<<'\t';
    std::cout<<"\n \n";
  }

  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof: "<<ii<<'\t';
    std::cout<<"Num_LD: "<<get_Num_LD(ii)<<'\t';
    std::cout<<"Num_LPS: "<<get_Num_LPS(ii)<<'\t';
    std::cout<<"Num_LPM: "<<get_Num_LPM(ii)<<'\n';
  }

  std::cout<<std::endl;

  std::cout<<"LPSN: \n";
  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof "<<ii<<'\n';
    for(int jj=0; jj<Num_LPS[ii]; ++jj)
      std::cout<<get_LPSN(ii, jj)<<'\t';
    std::cout<<"\n \n";
  }

  std::cout<<"LPMN: \n";
  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof "<<ii<<'\n';
    for(int jj=0; jj<Num_LPS[ii]; ++jj)
      std::cout<<get_LPMN(ii, jj)<<'\t';
    std::cout<<"\n \n";
  }

  std::cout<<"LocalMaster: \n";
  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof "<<ii<<'\n';
    for(int jj=0; jj<Num_LPM[ii]; ++jj)
      std::cout<<get_LocalMaster(ii, jj)<<'\t';
    std::cout<<"\n \n";
  }

  std::cout<<"LocalMasterSlave: \n";
  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof "<<ii<<'\n';
    for(int jj=0; jj<Num_LPM[ii]; ++jj)
      std::cout<<get_LocalMasterSlave(ii, jj)<<'\t';
    std::cout<<"\n \n";
  }
}

// EOF
