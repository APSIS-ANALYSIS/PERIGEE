#include "INodalBC.hpp"

INodalBC::INodalBC()
{}

INodalBC::~INodalBC()
{
  VEC_T::clean(ID);
}

void INodalBC::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"======== BC info ======="<<std::endl;
  if( get_num_dir_nodes() > 0 )
  {
    std::cout<<"Dirichlet nodes: "<<std::endl;
    for(unsigned int ii=0; ii<get_num_dir_nodes(); ++ii)
      std::cout<<get_dir_nodes(ii)<<'\t';
    std::cout<<std::endl;
  }

  if( get_num_per_nodes() > 0 )
  {
    std::cout<<"Periodic master - slave nodes: "<<std::endl;
    for(unsigned int ii=0; ii<get_num_per_nodes(); ++ii)
      std::cout<<get_per_master_nodes(ii)<<'\t'<<get_per_slave_nodes(ii)<<std::endl;
  }

  std::cout<<std::endl<<"ID array: "<<std::endl;
  for(unsigned int ii=0; ii<ID.size(); ++ii)
    std::cout<<ID[ii]<<'\t';
  std::cout<<'\n';
  std::cout<<std::endl<<"========================"<<std::endl;
}

void INodalBC::Create_ID(const unsigned int &num_node)
{
  ID.resize(num_node); VEC_T::shrink2fit(ID);

  for(unsigned int ii = 0; ii<ID.size(); ++ii) ID[ii] = ii;

  for(unsigned int ii = 0; ii<get_num_dir_nodes(); ++ii)
    ID[ get_dir_nodes(ii) ] = -1;

  for(unsigned int ii = 0; ii<get_num_per_nodes(); ++ii)
    ID[ get_per_slave_nodes(ii) ] = get_per_master_nodes(ii);
}

// EOF
