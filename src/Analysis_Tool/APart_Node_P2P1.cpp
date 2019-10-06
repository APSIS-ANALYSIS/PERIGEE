#include "APart_Node_P2P1.hpp"

APart_Node_P2P1::APart_Node_P2P1( const std::string &fileBaseName, 
    const int &rank, const ALocal_IEN * const &ien_p2 )
: APart_Node(fileBaseName, rank)
{
  // Make sure that the input IEN is for tet-10
  SYS_T::print_fatal_if(ien_p2->get_stride() != 10, 
      "Error: The input for APart_Node_P2P1 should be 10-node tet IEN.\n");

  const int nlocelem = ien_p2 -> get_nlocalele();

  // temporary container for p1 loc node and ghost node
  std::vector<int> node_loc_p1, node_ghost_p1;

  node_loc_p1.clear();
  node_ghost_p1.clear();

  // Loop over element's first 4 nodes and collect the local and ghost node
  for(int ee=0; ee<nlocelem; ++ee)
  {
    for(int ii=0; ii<4; ++ii)
    {
      const int loc = ien_p2 -> get_LIEN(ee, ii);
      if(loc >= nlocalnode) node_ghost_p1.push_back(local_to_global[loc]);
      else node_loc_p1.push_back(local_to_global[loc]);
    }
  }

  // Sort and unique node collections
  VEC_T::sort_unique_resize(node_loc_p1);
  VEC_T::sort_unique_resize(node_ghost_p1);
  
  // Assign them to the true containers
  local_to_global.clear();
  node_loc.clear();
  node_ghost.clear();

  node_loc = node_loc_p1;
  node_ghost = node_ghost_p1;

  local_to_global = node_loc;
  VEC_T::insert_end( local_to_global, node_ghost );

  // Update the number of nodes
  nlocalnode = static_cast<int>( node_loc.size() );
  nghostnode = static_cast<int>( node_ghost.size() );
  nlocghonode = nlocalnode + nghostnode;
  ntotalnode = nlocghonode - nbadnode;
}


APart_Node_P2P1::~APart_Node_P2P1()
{}


void APart_Node_P2P1::print_info() const
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
