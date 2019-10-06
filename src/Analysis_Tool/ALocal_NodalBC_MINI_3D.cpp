#include "ALocal_NodalBC_MINI_3D.hpp"

ALocal_NodalBC_MINI_3D::ALocal_NodalBC_MINI_3D( 
    const std::string &fileBaseName,
    const int &cpu_rank, const APart_Node * const &pNode_disp )
: ntotnode( pNode_disp -> get_ntotalnode() )
{
  // ----------------------------------------------------------------
  // Read in the geo2phy mapping
  hid_t file_id_g2p = H5Fopen("GlobalPartInfo.h5",H5F_ACC_RDONLY,H5P_DEFAULT);

  HDF5_Reader * h5r_g2p = new HDF5_Reader( file_id_g2p );

  std::vector<int> g2p;

  h5r_g2p -> read_intVector("/", "geo2phy", g2p);

  delete h5r_g2p; H5Fclose(file_id_g2p);
  // ----------------------------------------------------------------

  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const int nlocghonode = h5r->read_intScalar( "Local_Node", "nlocghonode" );
  const int nbubnode = pNode_disp -> get_nbubblenode();

  std::vector<int> LID;

  h5r->read_intVector( "/nbc", "LID", LID );

  // Format check 
  SYS_T::print_fatal_if( LID.size() % nlocghonode != 0 ,"Error:ALocal_NodalBC, LID length is not compatible with local and ghost node number. \n");

  SYS_T::print_fatal_if( LID.size() / nlocghonode != 4 ,"Error:ALocal_NodalBC, LID's dof should be 4 for ALocal_NodalBC_MINI_3D. \n");

  SYS_T::print_fatal_if( ntotnode != nlocghonode + nbubnode, "Error:ALocal_NodalBC, number of nodes does not match. Check the enriched APart_Node_ class.\n" );

  LID_p.resize( nlocghonode );
  LID_u.resize( 3 * ntotnode );

  for(int ii=0; ii<nlocghonode; ++ii) LID_p[ii] = LID[ii];

  for(int dd=0; dd<3; ++dd)
  {
    // Copy the continuous nodes from LID, and if not -1, map by g2p
    for(int ii=0; ii<nlocghonode; ++ii){
      if( LID[(dd+1)*nlocghonode + ii] >= 0 )
        LID_u[dd*ntotnode+ii] = g2p[ LID[(dd+1)*nlocghonode + ii] ];
      else
        LID_u[dd*ntotnode+ii] = -1;
    }
    // Append the bubble node
    for(int ii=0; ii<nbubnode; ++ii)
      LID_u[dd*ntotnode+nlocghonode+ii] = pNode_disp->get_bubble_node(ii);
  }

  std::vector<int> Num_LD;
  h5r->read_intVector( "/nbc", "Num_LD", Num_LD );

  Num_LD_u.resize(3);

  // Define the number of LD nodes. 
  Num_LD_p = Num_LD[0];
  Num_LD_u[0] = Num_LD[1]; 
  Num_LD_u[1] = Num_LD[2]; 
  Num_LD_u[2] = Num_LD[3];

  int size_ldn = Num_LD[0] + Num_LD[1] + Num_LD[2] + Num_LD[3];

  // Define LDN array
  LDN.clear();

  if(size_ldn > 0) h5r->read_intVector("/nbc", "LDN", LDN);

  SYS_T::print_fatal_if( int(LDN.size())!= size_ldn, "Error:ALocal_Nodalbc, LDN length does not match Num_LD.\n" );

  delete h5r; H5Fclose(file_id);

  // Generate the offsets
  LD_offset_u.resize(3);
  LD_offset_u[0] = Num_LD_p;
  LD_offset_u[1] = LD_offset_u[0] + Num_LD_u[0];
  LD_offset_u[2] = LD_offset_u[1] + Num_LD_u[1];

  // Update the LDN indices by g2p
  for(int ii=Num_LD_p; ii<(int)LDN.size(); ++ii)
    LDN[ii] = g2p[LDN[ii]];
}


ALocal_NodalBC_MINI_3D::~ALocal_NodalBC_MINI_3D()
{
  VEC_T::clean(LID_p);
  VEC_T::clean(LID_u);
  VEC_T::clean(Num_LD_u);
  VEC_T::clean(LD_offset_u);
  VEC_T::clean(LDN);
}


void ALocal_NodalBC_MINI_3D::print_info() const
{
  std::cout<<"ALocal_NodalBC: \n";
  std::cout<<"LID_p: \n";
  VEC_T::print(LID_p);
  std::cout<<"LID_u: \n";
  VEC_T::print(LID_u);

  std::cout<<"LDN: \n";
  VEC_T::print(LDN);

  std::cout<<"Num_LD_p: "<<Num_LD_p<<std::endl;
  std::cout<<"Num_LD_u: ";
  VEC_T::print(Num_LD_u);
  std::cout<<"LD_offset_u: ";
  VEC_T::print(LD_offset_u);
}


// EOF
