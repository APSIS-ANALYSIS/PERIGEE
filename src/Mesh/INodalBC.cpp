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
  if(num_dir_nodes.size() > 0)
  {
    std::cout<<"Dirichlet nodes: "<<std::endl;
    for(unsigned int ii=0; ii<num_dir_nodes.size(); ++ii)
    {
      std::cout<<"nbc_id " << ii << ": ";
      for(unsigned int jj=0; jj<num_dir_nodes[ii]; ++jj)
        std::cout<<dir_nodes[ii][jj]<<'\t';
      std::cout<<std::endl;
    }
  }

  if(num_per_nodes.size() > 0)
  {
    std::cout<<"Periodic master - slave nodes: "<<std::endl;
    for(unsigned int ii=0; ii<num_per_nodes.size(); ++ii)
    {
      std::cout<<"nbc_id " << ii << ": ";
      for(unsigned int jj=0; jj<num_per_nodes[ii]; ++jj)
        std::cout<<per_master_nodes[ii][jj]<<'\t'<<per_slave_nodes[ii][jj]<<std::endl;
    }
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
