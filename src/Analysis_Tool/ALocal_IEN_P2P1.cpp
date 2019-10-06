#include "ALocal_IEN_P2P1.hpp"

ALocal_IEN_P2P1::ALocal_IEN_P2P1( const std::string &fileBaseName,
    const int &cpu_rank,
    const APart_Node * const &node_u,
    const APart_Node * const &node_p )
: ALocal_IEN(fileBaseName, cpu_rank)
{
  // Make sure the read-in IEN data is for 10-node tet
  SYS_T::print_fatal_if( nLocBas != 10,
      "Error: The input for APart_Node_P2P1 should be 10-node tet IEN.\n");

  // Copy the local_to_global from the APart_Node_P2P1
  std::vector<int> l2g;
  const int len_l2g = node_p -> get_nlocghonode();
  l2g.resize( len_l2g );
  for(int ii=0; ii<len_l2g; ++ii)
    l2g[ii] = node_p -> get_local_to_global(ii);

  // temporary container
  std::vector<int> temp;
  temp.resize( 4 * nlocalele );
  
  // temp stores the location of the LIEN in the new local_to_global array
  for(int ee=0; ee<nlocalele; ++ee)
  {
    for(int ii=0; ii<4; ++ii)
    {
      const int val = LIEN[ee*10 + ii];
      temp[ee*4+ii] = VEC_T::get_pos( l2g, node_u->get_local_to_global(val) );
    }
  }
  
  // update data
  nLocBas = 4;
  LIEN.clear();
  LIEN = temp;
}


ALocal_IEN_P2P1::~ALocal_IEN_P2P1()
{}


void ALocal_IEN_P2P1::print_info() const
{
  std::cout<<"ALocal_IEN_P2P1: \n";

  std::cout<<"nlocalele = "<<nlocalele<<'\n';
  std::cout<<"nLocBas = "<<nLocBas<<'\n';

  std::vector<int> temp; temp.resize(nLocBas);
  for(int ee=0; ee < nlocalele; ++ee)
  {
    get_LIEN_e(ee, &temp[0]);
    std::cout<<ee<<" : \t";
    VEC_T::print(temp);
    std::cout<<'\n';
  }
}


// EOF
