#include "APart_Node_FSI.hpp"

APart_Node_FSI::APart_Node_FSI(const std::string &fileBaseName, const int &rank )
: APart_Node(fileBaseName, rank)
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

  const std::string gname("Local_Node");

  nlocalnode_fluid = h5r -> read_intScalar( gname.c_str(), "nlocalnode_fluid" );
  nlocalnode_solid = h5r -> read_intScalar( gname.c_str(), "nlocalnode_solid" );

  if(nlocalnode_fluid > 0)
    node_loc_fluid = h5r -> read_intVector( gname.c_str(), "node_loc_fluid" );

  if(nlocalnode_solid > 0)
    node_loc_solid = h5r -> read_intVector( gname.c_str(), "node_loc_solid" );

  H5Fclose( file_id );
}

void APart_Node_FSI::print_info() const
{
  APart_Node::print_info();
  std::cout<<"\n node_loc_solid: \n";
  VEC_T::print(node_loc_solid);
  std::cout<<"\n node_loc_fluid: \n";
  VEC_T::print(node_loc_fluid);
}

// EOF
