#include "Part_FSI_PV.hpp"

Part_FSI_PV::Part_FSI_PV( const IMesh * const &mesh_p,
    const IMesh * const &mesh_v,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const Map_Node_Index * const &mnindex_p,
    const Map_Node_Index * const &mnindex_v,
    const IIEN * const &IEN_p,
    const IIEN * const &IEN_v,
    const std::vector<double> &ctrlPts,
    const std::vector<int> &phytag,
    const std::vector<int> &node_f,
    const std::vector<int> &node_s,
    const int &in_start_idx_p, const int &in_start_idx_v,
    const int &in_cpu_rank, const int &in_cpu_size,
    const int &in_elemType,
    const bool isPrintInfo )
: cpu_rank(in_cpu_rank), cpu_size(in_cpu_size),
  dual_edge_ncommon( gpart->get_dual_edge_ncommon() ),
  nElem( mesh_p->get_nElem() ), probDim(3), elemType(in_elemType),
  nFunc_p( mesh_p->get_nFunc() ), 
  sDegree_p( mesh_p->get_s_degree() ), tDegree_p( mesh_p->get_t_degree() ),
  uDegree_p( mesh_p->get_u_degree() ), nLocBas_p( mesh_p->get_nLocBas() ),
  nFunc_v( mesh_v->get_nFunc() ), 
  sDegree_v( mesh_v->get_s_degree() ), tDegree_v( mesh_v->get_t_degree() ),
  uDegree_v( mesh_v->get_u_degree() ), nLocBas_v( mesh_v->get_nLocBas() )
{
  SYS_T::print_fatal_if( nElem != mesh_v -> get_nElem(), "Error: nElem in pressure mesh and velocity mesh should be the same.\n" );
  
  // 1. create local element partition based on epart
  elem_loc.clear();

  for(int ee=0; ee<nElem; ++ee)
  {
    if( gpart->get_epart(ee) == cpu_rank ) elem_loc.push_back(ee);
  }

  VEC_T::shrink2fit(elem_loc);
  nlocalele = static_cast<int>( elem_loc.size() );

  elem_phy_tag.resize( nlocalele );
  for(int ii=0; ii<nlocalele; ++ii) elem_phy_tag[ii] = phytag[ elem_loc[ii] ];

  // 2. Local node partition based on npart
  node_loc_p.clear(); node_loc_v.clear();
  node_loc_original_p.clear(); node_loc_original_v.clear();
  for(int nn=0; nn<nFunc_p; ++nn)
  {
    if( gpart->get_npart(nn) == cpu_rank )
    {
      node_loc_p.push_back(nn);
      node_loc_original_p.push_back(nn);
    }
  }

  VEC_T::shrink2fit(node_loc_p); VEC_T::shrink2fit(node_loc_original_p);
  nlocalnode_p = static_cast<int>( node_loc_p.size() );

  for(int nn=0; nn<nFunc_v; ++nn)
  {
    if( gpart->get_npart(nn + nFunc_p) == cpu_rank )
    {
      node_loc_v.push_back(nn);
      node_loc_original_v.push_back(nn);
    }
  }
  VEC_T::shrink2fit(node_loc_v); VEC_T::shrink2fit(node_loc_original_v);
  nlocalnode_v = static_cast<int>( node_loc_v.size() );

  if(isPrintInfo)
  {
    std::cout<<"-- proc "<<cpu_rank<<" -- elem_loc & node_loc arrays generated. \n";
    std::cout<<"-- proc "<<cpu_rank<<" local element number: "<<elem_loc.size()<<std::endl;
    std::cout<<"-- proc "<<cpu_rank<<" nlocalnode_p "<<nlocalnode_p<<" nlocalnode_v "<<nlocalnode_v<<std::endl;
  }

  for(int ii=0; ii<nlocalnode_p; ++ii)
    node_loc_p[ii] = mnindex_p->get_old2new( node_loc_p[ii] );

  for(int ii=0; ii<nlocalnode_v; ++ii)
    node_loc_v[ii] = mnindex_v->get_old2new( node_loc_v[ii] );

  std::vector<int> node_tot_p; node_tot_p.clear();
  std::vector<int> node_tot_v; node_tot_v.clear();
  for(int ee=0; ee<nlocalele; ++ee)
  {
    for(int ii=0; ii<nLocBas_p; ++ii)
    {
      int temp_node = IEN_p ->get_IEN( elem_loc[ee], ii );
      node_tot_p.push_back( mnindex_p->get_old2new(temp_node) );
    }
    for(int ii=0; ii<nLocBas_v; ++ii)
    {
      int temp_node = IEN_v -> get_IEN(elem_loc[ee], ii);
      node_tot_v.push_back( mnindex_v->get_old2new(temp_node) );
    }
  }
  VEC_T::sort_unique_resize( node_tot_p );
  VEC_T::sort_unique_resize( node_tot_v );
  ntotalnode_p = static_cast<int>( node_tot_p.size() );
  ntotalnode_v = static_cast<int>( node_tot_v.size() );

  // collect ghost node from two fields separately
  node_ghost_p.clear();
  for(int ii=0; ii<ntotalnode_p; ++ii)
  {
    if( !VEC_T::is_invec(node_loc_p, node_tot_p[ii]) )
      node_ghost_p.push_back(node_tot_p[ii]);
  }

  node_ghost_v.clear();
  for(int ii=0; ii<ntotalnode_v; ++ii)
  {
    if( !VEC_T::is_invec(node_loc_v, node_tot_v[ii] ) )
      node_ghost_v.push_back(node_tot_v[ii]);
  }
  VEC_T::shrink2fit(node_ghost_p); VEC_T::shrink2fit(node_ghost_v);
  nghostnode_p = static_cast<int>( node_ghost_p.size() );
  nghostnode_v = static_cast<int>( node_ghost_v.size() );

}


Part_FSI_PV::~Part_FSI_PV()
{}



// EOF
