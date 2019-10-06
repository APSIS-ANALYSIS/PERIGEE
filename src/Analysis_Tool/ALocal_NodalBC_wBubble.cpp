#include "ALocal_NodalBC_wBubble.hpp"

ALocal_NodalBC_wBubble::ALocal_NodalBC_wBubble( 
    const std::string &fileBaseName,
    const int &cpu_rank, const APart_Node * const &pNode )
: ALocal_NodalBC(fileBaseName, cpu_rank),
  ntotnode( pNode -> get_ntotalnode() )
{
  SYS_T::print_fatal("Warning: this class needs more testing.\n");
  // ----------------------------------------------------------------
  // Read in the geo2phy mapping
  hid_t file_id_g2p = H5Fopen("GlobalPartInfo.h5",H5F_ACC_RDONLY,H5P_DEFAULT);

  HDF5_Reader * h5r_g2p = new HDF5_Reader( file_id_g2p );

  std::vector<int> g2p;

  h5r_g2p -> read_intVector("/", "geo2phy", g2p);

  delete h5r_g2p; H5Fclose(file_id_g2p);
  // ----------------------------------------------------------------

  // Copy LID to temporary vector
  std::vector<int> temp = LID;

  // Clean LID
  LID.clear();

  const int nbub = pNode -> get_nbubblenode();

  LID.resize( ntotnode * dof );

  // Now we first use g2p to update the temp to the phy indices
  for(int ii=0; ii<(int)temp.size(); ++ii){
    if( temp[ii] >= 0 ) temp[ii] = g2p[ temp[ii] ];
  }

  for(int dd=0; dd<dof; ++dd)
  {
    // the first nlocghonode nodes are copied from temp
    for(int ii=0; ii<nlocghonode; ++ii)
      LID[dd*ntotnode+ii] = temp[dd*nlocghonode+ii];

    // the next nbub nodes are the bubble nodes
    for(int ii=0; ii<nbub; ++ii)
      LID[dd*ntotnode+nlocghonode+ii] = pNode->get_bubble_node(ii);
  }

  // Update LDN
  for(int ii=0; ii<(int)LDN.size(); ++ii)
    LDN[ii] = g2p[ LDN[ii] ];

  // Update LPSN
  for(int ii=0; ii<(int)LPSN.size(); ++ii)
    LPSN[ii] = g2p[ LPSN[ii] ];

  // Update LPMN
  for(int ii=0; ii<(int)LPMN.size(); ++ii)
    LPMN[ii] = g2p[ LPMN[ii] ];

  // Update LocalMaster
  for(int ii=0; ii<(int)LocalMaster.size(); ++ii)
    LocalMaster[ii] = g2p[ LocalMaster[ii] ];

  // Update LocalMasterSlave
  for(int ii=0; ii<(int)LocalMasterSlave.size(); ++ii)
    LocalMasterSlave[ii] = g2p[ LocalMasterSlave[ii] ];

}


ALocal_NodalBC_wBubble::~ALocal_NodalBC_wBubble()
{}


void ALocal_NodalBC_wBubble::print_info() const
{
  std::cout<<"ALocal_NodalBC_wBubble: \n";

  std::cout<<"LID: \n";

  for(int ii=0; ii<dof; ++ii)
  {
    std::cout<<"dof "<<ii<<'\n';
    for(int jj=0; jj<ntotnode; ++jj)
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
