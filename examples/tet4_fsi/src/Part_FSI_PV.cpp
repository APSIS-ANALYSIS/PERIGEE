#include "Part_FSI_PV.hpp"

Part_FSI_PV::Part_FSI_PV( const IMesh * const &mesh_p,
    const IMesh * const &mesh_v,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex_p,
    const Map_Node_Index * const &mnindex_v,
    const IIEN * const &IEN_p,
    const IIEN * const &IEN_v,
    const std::vector<double> &ctrlPts,
    const std::vector<int> &phytag,
    const std::vector<int> &p_node_f,
    const std::vector<int> &p_node_s,
    const std::vector<int> &v_node_f,
    const std::vector<int> &v_node_s,
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
    std::cout<<"-- proc "<<cpu_rank<<" elem_loc & node_loc arrays generated. \n";
    std::cout<<"-- proc "<<cpu_rank<<" local element number: "<<elem_loc.size()<<std::endl;
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

  SYS_T::print_fatal_if( nghostnode_p + nlocalnode_p != ntotalnode_p,
      "Error: the pressure node partition may contain bad node. \n");

  SYS_T::print_fatal_if( nghostnode_v + nlocalnode_v != ntotalnode_v,
      "Error: the velocity node partition may contain bad node. \n");

  if( isPrintInfo )
  {
    std::cout<<"-- proc "<<cpu_rank<<" ntotalnode_p "<<ntotalnode_p<<" ntotalnode_v "<<ntotalnode_v<<std::endl;
    std::cout<<"-- proc "<<cpu_rank<<" nlocalnode_p "<<nlocalnode_p<<" nlocalnode_v "<<nlocalnode_v<<std::endl;
    std::cout<<"-- proc "<<cpu_rank<<" nghostnode_p "<<nghostnode_p<<" nghostnode_v "<<nghostnode_v<<std::endl;
  }

  // collect pressure and velocity nodes in local node vector
  node_loc_p_fluid.clear(); node_loc_p_solid.clear();
  
  for(int ii=0; ii<nlocalnode_p; ++ii)
  {
    if( VEC_T::is_invec(p_node_f, node_loc_original_p[ii]) ) node_loc_p_fluid.push_back(ii);

    if( VEC_T::is_invec(p_node_s, node_loc_original_p[ii]) ) node_loc_p_solid.push_back(ii);
  }
  
  nlocalnode_p_fluid = static_cast<int>( node_loc_p_fluid.size() );
  nlocalnode_p_solid = static_cast<int>( node_loc_p_solid.size() );

  node_loc_v_fluid.clear(); node_loc_v_solid.clear();

  for(int ii=0; ii<nlocalnode_v; ++ii)
  {
    if( VEC_T::is_invec(v_node_f, node_loc_original_v[ii]) ) node_loc_v_fluid.push_back(ii);

    if( VEC_T::is_invec(v_node_s, node_loc_original_v[ii]) ) node_loc_v_solid.push_back(ii);
  }

  nlocalnode_v_fluid = static_cast<int>( node_loc_v_fluid.size() );
  nlocalnode_v_solid = static_cast<int>( node_loc_v_solid.size() );

  // local to global mapping
  local_to_global_p = node_loc_p;
  VEC_T::insert_end(local_to_global_p, node_ghost_p);

  local_to_global_v = node_loc_v;
  VEC_T::insert_end(local_to_global_v, node_ghost_v);

  VEC_T::shrink2fit(local_to_global_p);
  VEC_T::shrink2fit(local_to_global_v);

  nlocghonode_p = static_cast<int>( local_to_global_p.size() );
  nlocghonode_v = static_cast<int>( local_to_global_v.size() );

  if( isPrintInfo ) std::cout<<"-- proc "<<cpu_rank<<" local_to_global generated. \n";

  // LIEN
  pLIEN.resize(nlocalele * nLocBas_p);
  vLIEN.resize(nlocalele * nLocBas_v);

  for(int ee=0; ee<nlocalele; ++ee)
  {
    for(int ii=0; ii<nLocBas_p; ++ii)
    {
      int global_index = IEN_p -> get_IEN( elem_loc[ee], ii );
      global_index = mnindex_p -> get_old2new( global_index );
      const auto lien_ptr = find( local_to_global_p.begin(),local_to_global_p.end(),
          global_index );
      
      SYS_T::print_fatal_if( lien_ptr == local_to_global_p.end(), "Error: unable to locate an index.\n");
      
      pLIEN[ee*nLocBas_p + ii] = lien_ptr - local_to_global_p.begin();
    }

    for(int ii=0; ii<nLocBas_v; ++ii)
    {
      int global_index = IEN_v -> get_IEN( elem_loc[ee], ii );
      global_index = mnindex_v -> get_old2new(global_index);
      const auto lien_ptr = find( local_to_global_v.begin(), local_to_global_v.end(),
          global_index );
      
      SYS_T::print_fatal_if( lien_ptr == local_to_global_p.end(), "Error: unable to locate an index.\n");
      
      vLIEN[ee*nLocBas_v + ii] = lien_ptr - local_to_global_v.begin();
    }
  }

  if(isPrintInfo) std::cout<<"-- proc "<<cpu_rank<<" LIEN generated. \n";

  // 7. local copy of control points based on velocity mesh partition
  ctrlPts_x_loc.resize(nlocghonode_v);
  ctrlPts_y_loc.resize(nlocghonode_v);
  ctrlPts_z_loc.resize(nlocghonode_v);

  for(int ii=0; ii<nlocghonode_v; ++ii)
  {
    int aux_index = local_to_global_v[ii]; // new global index
    aux_index = mnindex_v->get_new2old(aux_index); // back to old global index
    ctrlPts_x_loc[ii] = ctrlPts[3*aux_index + 0];
    ctrlPts_y_loc[ii] = ctrlPts[3*aux_index + 1];
    ctrlPts_z_loc[ii] = ctrlPts[3*aux_index + 2];
  }

  VEC_T::shrink2fit(ctrlPts_x_loc);
  VEC_T::shrink2fit(ctrlPts_y_loc);
  VEC_T::shrink2fit(ctrlPts_z_loc);

  if(isPrintInfo) std::cout<<"-- proc "<<cpu_rank<<" Local control points generated. \n";
}


Part_FSI_PV::~Part_FSI_PV()
{}



// EOF
