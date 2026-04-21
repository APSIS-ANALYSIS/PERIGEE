#include "APart_Node_Rotated.hpp"

APart_Node_Rotated::APart_Node_Rotated(const std::string &fileBaseName, int rank )
: APart_Node(fileBaseName, rank)
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  auto h5r = SYS_T::make_unique<HDF5_Reader>(fName);

  const std::string gname("Local_Node");

  nlocalnode_fixed = h5r -> read_intScalar( gname.c_str(), "nlocalnode_fixed" );
  nlocalnode_rotated = h5r -> read_intScalar( gname.c_str(), "nlocalnode_rotated" );

  if(nlocalnode_fixed > 0)
    node_loc_fixed = h5r -> read_intVector( gname.c_str(), "node_loc_fixed" );

  if(nlocalnode_rotated > 0)
    node_loc_rotated = h5r -> read_intVector( gname.c_str(), "node_loc_rotated" );
}

APart_Node_Rotated::APart_Node_Rotated(const HDF5_Reader * const &h5r)
: APart_Node(h5r)
{
  const std::string gname("Local_Node");

  nlocalnode_fixed = h5r -> read_intScalar( gname.c_str(), "nlocalnode_fixed" );
  nlocalnode_rotated = h5r -> read_intScalar( gname.c_str(), "nlocalnode_rotated" );

  if(nlocalnode_fixed > 0)
    node_loc_fixed = h5r -> read_intVector( gname.c_str(), "node_loc_fixed" );

  if(nlocalnode_rotated > 0)
    node_loc_rotated = h5r -> read_intVector( gname.c_str(), "node_loc_rotated" );
}

void APart_Node_Rotated::print_info() const
{
  APart_Node::print_info();
  std::cout<<"\n node_loc_rotated: \n";
  VEC_T::print(node_loc_rotated);
  std::cout<<"\n node_loc_fixed: \n";
  VEC_T::print(node_loc_fixed);
}

// EOF
