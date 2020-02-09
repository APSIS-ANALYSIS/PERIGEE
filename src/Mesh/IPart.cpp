#include "IPart.hpp"

void IPart::GenPart( const int &nElem,
    const int &nFunc, const int &nLocBas, const int &cpu_rank,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const bool &isPrintInfo,
    std::vector<int> &elem_loc,
    int &nlocalele,
    std::vector<int> &node_loc,
    std::vector<int> &node_loc_original,
    std::vector<int> &node_ghost,
    std::vector<int> &local_to_global,
    int &nlocalnode,
    int &nghostnode,
    int &ntotalnode,
    int &nbadnode,
    int &nlocghonode,
    std::vector<double> &ctrlPts_x_loc,
    std::vector<double> &ctrlPts_y_loc,
    std::vector<double> &ctrlPts_z_loc,
    std::vector<std::vector<int> > &LIEN ) const
{
  // 1. Create local partition based on epart & npart
  elem_loc.clear(); node_loc.clear();

  for( int e=0; e<nElem; ++e )
  {
    if( gpart->get_epart(e) == cpu_rank ) elem_loc.push_back(e);
  }
  VEC_T::shrink2fit(elem_loc);
  nlocalele = static_cast<int>( elem_loc.size() );

  for( int n=0; n<nFunc; ++n )
  {
    if( gpart->get_npart(n) == cpu_rank )
    {
      node_loc.push_back(n);
      node_loc_original.push_back(n);
    }
  }
  VEC_T::shrink2fit(node_loc);
  VEC_T::shrink2fit(node_loc_original);
  nlocalnode = static_cast<int>( node_loc.size() );

  if(isPrintInfo)
  {
    std::cout<<"-- proc "<<cpu_rank<<" -- elem_loc & node_loc arrays generated. \n";
    std::cout<<"-- proc "<<cpu_rank<<" local element number: "<<elem_loc.size()<<std::endl;
  }

  // 2. Reorder node_loc indices
  for( int ii=0; ii<nlocalnode; ++ii )
    node_loc[ii] = mnindex->get_old2new( node_loc[ii] );

  // 3. Generate node_tot, all the nodes needed by the subdomain
  std::vector<int> node_tot; node_tot.clear();
  for( int e=0; e<nlocalele; ++e )
  {
    for( int ii=0; ii<nLocBas; ++ii )
    {
      int temp_node = IEN->get_IEN(elem_loc[e], ii);
      temp_node = mnindex->get_old2new(temp_node);
      node_tot.push_back( temp_node );
    }
  }
  VEC_T::sort_unique_resize( node_tot );
  ntotalnode = static_cast<int>( node_tot.size() );

  // 4. node_ghost
  node_ghost.clear();
  for( int ii = 0; ii<ntotalnode; ++ii )
  {
    if( !VEC_T::is_invec(node_loc, node_tot[ii]) )
      node_ghost.push_back(node_tot[ii]);
  }
  VEC_T::shrink2fit(node_ghost);
  nghostnode = static_cast<int>( node_ghost.size() );

  nbadnode = 0;
  if( nghostnode + nlocalnode != ntotalnode )
  {
    std::vector<int> badnode;
    std::vector<int>::iterator badnode_it;
    for( int n=0; n<nlocalnode; ++n )
    {
      badnode_it = find( node_tot.begin(), node_tot.end(), node_loc[n] );
      if(badnode_it == node_tot.end()) badnode.push_back( node_loc[n] );
    }
    nbadnode = static_cast<int>( badnode.size() );
    if( nghostnode + nlocalnode != ntotalnode + nbadnode )
    {
      std::cerr<<"ERROR: The local node partition is wrong in proc "<<cpu_rank<<std::endl;
      std::cerr<<"       ghost node,  local node, total node, bad node: "<<std::endl;
      std::cerr<<"      "<<nghostnode<<'\t'<<nlocalnode<<'\t'<<ntotalnode<<'\t';
      std::cerr<<nbadnode<<std::endl;
      exit(EXIT_FAILURE);
    }
    std::cout<<"WARNING: The partition is poor for proecssor: "<<cpu_rank<<std::endl;
  }

  if( isPrintInfo )
  {
    std::cout<<"-- proc "<<cpu_rank;
    std::cout<<" -- ntotalnode: "<<ntotalnode;
    std::cout<<" -- nlocalnode: "<<nlocalnode;
    std::cout<<" -- nghostnode: "<<nghostnode;
    std::cout<<" -- nbadnode: "<<nbadnode<<std::endl;
  }

  // 5. local_to_global mapping
  local_to_global = node_loc;
  for( int n=0; n<nghostnode; ++n )
    local_to_global.push_back( node_ghost[n] );

  VEC_T::shrink2fit(local_to_global);
  nlocghonode = static_cast<int>( local_to_global.size() );

  if( isPrintInfo )
    std::cout<<"-- proc "<<cpu_rank<<" local_to_global generated. \n";

  // 6. LIEN
  LIEN.resize( nlocalele );
  for(int e=0; e<nlocalele; ++e) LIEN[e].resize(nLocBas);

  std::vector<int>::iterator lien_ptr;
  int global_index;
  for(int e=0; e<nlocalele; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      global_index = IEN->get_IEN(elem_loc[e], i);
      global_index = mnindex->get_old2new(global_index);
      lien_ptr = find( local_to_global.begin(), local_to_global.end(),
          global_index);
      LIEN[e][i] = lien_ptr - local_to_global.begin();
    }
  }
  if(isPrintInfo) std::cout<<"-- proc "<<cpu_rank<<" LIEN generated. \n";

  // 7. Control points coordinates
  ctrlPts_x_loc.resize(nlocghonode);
  ctrlPts_y_loc.resize(nlocghonode);
  ctrlPts_z_loc.resize(nlocghonode);

  for(int ii=0; ii<nlocghonode; ++ii)
  {
    int aux_index = local_to_global[ii]; // new global index
    aux_index = mnindex->get_new2old(aux_index); // back to old global index
    ctrlPts_x_loc[ii] = ctrlPts[3*aux_index + 0];
    ctrlPts_y_loc[ii] = ctrlPts[3*aux_index + 1];
    ctrlPts_z_loc[ii] = ctrlPts[3*aux_index + 2];
  }

  VEC_T::shrink2fit(ctrlPts_x_loc);
  VEC_T::shrink2fit(ctrlPts_y_loc);
  VEC_T::shrink2fit(ctrlPts_z_loc);

  if(isPrintInfo)
    std::cout<<"-- proc "<<cpu_rank<<" Local control points generated. \n";
}

// EOF
