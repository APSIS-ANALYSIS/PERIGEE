#include "NodalBC_3D.hpp"

NodalBC_3D::NodalBC_3D( const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_dir_nodes = 0;
  num_per_nodes = 0;

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D: No nodal BC is generated. \n";
}


NodalBC_3D::NodalBC_3D( const int &face_idx,
    const std::vector<std::vector<int> > &in_pt, const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;
  num_dir_nodes = in_pt[face_idx].size();

  dir_nodes.resize(num_dir_nodes);

  for(unsigned int ii=0; ii<num_dir_nodes; ++ii)
    dir_nodes[ii] = static_cast<unsigned int>( in_pt[face_idx][ii] );

  Create_ID(nFunc);

  std::cout<<"===> NodalBC_3D: Nodes from face "<<face_idx<<" is generated.\n";
}


NodalBC_3D::NodalBC_3D( const std::vector<int> &face_idx,
    const std::vector<std::vector<int> > &in_pt, const int &nFunc )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();
  num_per_nodes = 0;

  const unsigned int number_of_face = face_idx.size();

  for(unsigned int ii=0; ii<number_of_face; ++ii)
  {
    const unsigned int fa_idx = face_idx[ii];

    if( fa_idx >= in_pt.size() )
      SYS_T::print_fatal("Error: NodalBC_3D:: the face_idx is wrong.\n");

    const unsigned int face_ii_len = in_pt[fa_idx].size();
    for(unsigned int jj=0; jj<face_ii_len; ++jj)
      dir_nodes.push_back( static_cast<unsigned int>(in_pt[fa_idx][jj]) );
  }

  VEC_T::sort_unique_resize( dir_nodes ); // there may be redundant nodes

  num_dir_nodes = dir_nodes.size();

  Create_ID(nFunc);

  std::cout<<"===> NodalBC_3D: Nodes from faces ";
  for(unsigned int ii=0; ii<number_of_face; ++ii)
    std::cout<<face_idx[ii]<<'\t';
  std::cout<<" is generated. \n";
}


NodalBC_3D::~NodalBC_3D()
{}


NodalBC_3D::NodalBC_3D( const std::vector<int> &face_idx,
    const std::vector<std::vector<int> > &in_pt, const int &nFunc,
    const int &type )
{
  dir_nodes.clear();
  per_slave_nodes.clear();
  per_master_nodes.clear();

  num_dir_nodes = 0;
  num_per_nodes = 0;

  switch( type )
  {
    case 1:
      BC_type_1( face_idx, in_pt, nFunc );
      break;
    case 2:
      BC_type_2( face_idx, in_pt, nFunc );
      break;
    default:
      std::cerr<<"Error: NodalBC_3D with bc_type = "<<type<<" is not implemented. \n";
      exit(EXIT_FAILURE);
  }

  Create_ID( nFunc );
  std::cout<<"===> NodalBC_3D: Nodes from faces ";
  for(unsigned int ii=0; ii<face_idx.size(); ++ii)
    std::cout<<face_idx[ii]<<'\t';
  std::cout<<" with type = "<<type<<" is generated. \n";
}


void NodalBC_3D::BC_type_1( const std::vector<int> &face_idx,
    const std::vector<std::vector<int> > &inpt, const int &nFunc )
{
  const unsigned int num_file = face_idx.size();

  if( num_file != 2 ) SYS_T::print_fatal("Error: NodalBC_3D::BC_type_1 the number of input faces is wrong.\n");

  // The first face gives the Dirichlet nodes
  // conceptually this is the face where the arterial strip is fixed
  int fa_idx = face_idx[0];
  const unsigned int num_dn = inpt[fa_idx].size();
  for(unsigned int ii=0; ii<num_dn; ++ii)
    dir_nodes.push_back( static_cast<unsigned int>( inpt[fa_idx][ii] ) );

  num_dir_nodes = dir_nodes.size();
  
  // The second face gives the traction force applied face
  fa_idx = face_idx[1];
  const unsigned int num_pn = inpt[fa_idx].size(); 

  for(unsigned int ii=1; ii<num_pn; ++ii)
  {
    per_slave_nodes.push_back( static_cast<unsigned int>( inpt[fa_idx][ii]) );
    per_master_nodes.push_back( static_cast<unsigned int>( inpt[fa_idx][0]) );
  }

  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> Dirichlet nodes from face "<<face_idx[0];
  std::cout<<" and Master-slave from face "<<face_idx[1]<<" with 0th node master.\n";
}


void NodalBC_3D::BC_type_2( const std::vector<int> &face_idx,
    const std::vector<std::vector<int> > &inpt, const int &nFunc )
{
  const unsigned int num_file = face_idx.size();

  for( unsigned int ii = 0; ii<num_file; ++ii )
  {
    const unsigned int fa_idx = face_idx[ii];
    const unsigned int num_pn = inpt[fa_idx].size();

    for(unsigned int jj=1; jj<num_pn; ++jj)
    {
      per_slave_nodes.push_back( static_cast<unsigned int>( inpt[fa_idx][jj]) );
      per_master_nodes.push_back( static_cast<unsigned int>( inpt[fa_idx][0]) );
    }
  }

  num_per_nodes = per_slave_nodes.size();

  std::cout<<"-----> Master-Slave from face ";
  for(unsigned int ii=0; ii<num_file; ++ii) std::cout<<face_idx[ii]<<'\t';
  std::cout<<" with 0th node as master on these faces. \n";
}

// EOF
