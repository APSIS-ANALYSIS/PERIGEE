#include "ALocal_NodalBC_P2P1_3D.hpp"

ALocal_NodalBC_P2P1_3D::ALocal_NodalBC_P2P1_3D( 
    const std::string &fileBaseName, const int &cpu_rank, 
    const APart_Node * const &pNode_pres )
{
  // ----------------------------------------------------------------
  // Read in the nlp in GlobalPartInfo.h5
  hid_t file_id_g = H5Fopen("GlobalPartInfo.h5",H5F_ACC_RDONLY,H5P_DEFAULT);

  HDF5_Reader * h5r_g = new HDF5_Reader( file_id_g );

  std::vector<int> np_loc, np_gho, mapper;

  h5r_g -> read_intVector("/", "np_loc", np_loc);
  h5r_g -> read_intVector("/", "np_gho", np_gho);
  h5r_g -> read_intVector("/", "g2p", mapper);

  delete h5r_g; H5Fclose(file_id_g);
  // ----------------------------------------------------------------
  // Read in the LID array 
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  nlgnode = h5r->read_intScalar( "Local_Node", "nlocghonode" );

  std::vector<int> LID;

  h5r->read_intVector( "/nbc", "LID", LID );

  // Format check
  SYS_T::print_fatal_if( LID.size() % nlgnode != 0 ,"Error:ALocal_NodalBC, LID length is not compatible with local and ghost node number. \n");

  SYS_T::print_fatal_if( LID.size() / nlgnode != 4 ,"Error:ALocal_NodalBC, LID's dof should be 4 for ALocal_NodalBC_MINI_3D. \n");

  SYS_T::print_fatal_if( np_loc[cpu_rank] != pNode_pres->get_nlocalnode(), "Error: the number of pressure local node does not match with the one in APart_Node_P2P1. \n");

  SYS_T::print_fatal_if( np_gho[cpu_rank] != pNode_pres->get_nghostnode(), "Error: the number of pressure ghost node does not match with the one in APart_Node_P2P1. \n");

  // Create LID_p
  LID_p.resize( np_loc[cpu_rank] + np_gho[cpu_rank] );

  for(int ii=0; ii<np_loc[cpu_rank] + np_gho[cpu_rank]; ++ii)
  {
    int p2idx = pNode_pres -> get_local_to_global( ii );
    LID_p[ii] = mapper[ p2idx ];
  }

  // Create LID_u
  LID_u.resize(3 * nlgnode);
  for(int dd=0; dd<3; ++dd)
  {
    for(int ii=0; ii<nlgnode; ++ii)
      LID_u[dd*nlgnode+ii] = LID[ (dd+1) * nlgnode + ii ];
  }
  
  // Num_LD
  std::vector<int> Num_LD;
  h5r->read_intVector( "/nbc", "Num_LD", Num_LD );

  Num_LD_p = Num_LD[0];
  Num_LD_u[0] = Num_LD[1];
  Num_LD_u[1] = Num_LD[2];
  Num_LD_u[2] = Num_LD[3];

  SYS_T::print_fatal_if(Num_LD_p != 0, "Error: We do not allow apply Nodal BC for pressure here.\n");

  int size_ldn = Num_LD[0] + Num_LD[1] + Num_LD[2] + Num_LD[3];

  // Define LDN array
  LDN.clear();

  if(size_ldn > 0) h5r->read_intVector("/nbc", "LDN", LDN);

  SYS_T::print_fatal_if( int(LDN.size())!= size_ldn, "Error:ALocal_Nodalbc, LDN length does not match Num_LD.\n" );

  delete h5r; H5Fclose(file_id);

  // Generate the offsets
  LD_offset_u[0] = Num_LD_p;
  LD_offset_u[1] = LD_offset_u[0] + Num_LD_u[0];
  LD_offset_u[2] = LD_offset_u[1] + Num_LD_u[1];
}


ALocal_NodalBC_P2P1_3D::~ALocal_NodalBC_P2P1_3D()
{
  VEC_T::clean(LID_p);
  VEC_T::clean(LID_u);
  VEC_T::clean(LDN);
}


void ALocal_NodalBC_P2P1_3D::print_info() const
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
  std::cout<<Num_LD_u[0]<<'\t';
  std::cout<<Num_LD_u[1]<<'\t';
  std::cout<<Num_LD_u[2]<<'\n';
  std::cout<<"LD_offset_u: ";
  std::cout<<LD_offset_u[0]<<'\t';
  std::cout<<LD_offset_u[1]<<'\t';
  std::cout<<LD_offset_u[2]<<'\n';
}

// EOF
