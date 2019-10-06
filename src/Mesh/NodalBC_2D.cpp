#include "NodalBC_2D.hpp"

NodalBC_2D::NodalBC_2D( const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;
  num_per_nodes = 0;

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_2D: No nodal BC is generated. \n";
}


NodalBC_2D::NodalBC_2D( const int &ed_idx,
    const std::vector<std::vector<int> > &in_pt, const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
  num_dir_nodes = in_pt[ed_idx].size();

  dir_nodes.resize(num_dir_nodes);

  for(unsigned int ii=0; ii<num_dir_nodes; ++ii)
    dir_nodes[ii] = static_cast<unsigned int>( in_pt[ed_idx][ii] );
  
  Create_ID(nFunc);

  std::cout<<"===> NodalBC_2D: Nodes from edge "<<ed_idx<<" is generated. \n";
}


NodalBC_2D::NodalBC_2D( const std::vector<int> &edge_idx,
    const std::vector<std::vector<int> > &in_pt, const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  const unsigned int number_of_edge = edge_idx.size();

  for(unsigned int ii=0; ii<number_of_edge; ++ii)
  {
    const unsigned int ed_idx = edge_idx[ii];
    
    if( ed_idx >= in_pt.size() )
      SYS_T::print_fatal("Error: NodalBC_2D:: the edge_idx is wrong. \n");
    
    const unsigned int edge_ii_len = in_pt[ed_idx].size();
    for(unsigned int jj=0; jj<edge_ii_len; ++jj)
      dir_nodes.push_back( static_cast<unsigned int>(in_pt[ed_idx][jj]) );
  }

  VEC_T::sort_unique_resize( dir_nodes ); // there may be redundant nodes

  num_dir_nodes = dir_nodes.size();

  Create_ID(nFunc);

  std::cout<<"===> NodalBC_2D: Nodes from edges ";
  for(unsigned int ii=0; ii<number_of_edge; ++ii)
    std::cout<<edge_idx[ii]<<'\t';
  std::cout<<" is generated. \n";
}


NodalBC_2D::~NodalBC_2D()
{}


// EOF
