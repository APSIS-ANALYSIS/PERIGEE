#include "APart_Node.hpp"

APart_Node::APart_Node( const std::string &fbasename, const int &rank )
: cpu_rank( rank )
{
  const std::string fName = SYS_T::gen_partfile_name( fbasename, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );
  
  nlocalnode = h5r->read_intScalar("Local_Node", "nlocalnode");
  nghostnode = h5r->read_intScalar("Local_Node", "nghostnode");
  nbadnode   = h5r->read_intScalar("Local_Node", "nbadnode");
  nlocghonode = h5r->read_intScalar("Local_Node", "nlocghonode");
  ntotalnode = h5r->read_intScalar("Local_Node", "ntotalnode");

  local_to_global = h5r->read_intVector("Local_Node", "local_to_global");
  
  if( nghostnode > 0)
    node_ghost = h5r->read_intVector("Local_Node", "node_ghost");
  else
    node_ghost.clear();
  
  node_loc = h5r->read_intVector("Local_Node", "node_loc");

  dof = h5r->read_intScalar("Global_Mesh_Info", "dofNum");

  delete h5r; H5Fclose( file_id );
}

APart_Node::~APart_Node()
{
  VEC_T::clean(local_to_global);
  VEC_T::clean(node_ghost);
  VEC_T::clean(node_loc);
}

void APart_Node::print_info() const
{
  std::cout<<"cpu "<<cpu_rank<<" node info: \n";
  std::cout<<"dof of this mesh "<<dof<<std::endl;
  std::cout<<"nlocalnode: "<<nlocalnode<<std::endl;
  std::cout<<"nghostnode: "<<nghostnode<<std::endl;
  std::cout<<"nbadnode: "<<nbadnode<<std::endl;
  std::cout<<"nlocghonode: "<<nlocghonode<<std::endl;
  std::cout<<"ntotalnode: "<<ntotalnode<<std::endl;
  std::cout<<" local_to_global: \n";
  VEC_T::print(local_to_global);
  std::cout<<"\n node_ghost: \n";
  VEC_T::print(node_ghost);
  std::cout<<"\n node_loc: \n";
  VEC_T::print(node_loc);
}

// EOF
