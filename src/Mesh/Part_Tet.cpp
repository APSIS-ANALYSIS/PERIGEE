#include "Part_Tet.hpp"

Part_Tet::Part_Tet(
    const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const int &in_cpu_rank, const int &in_cpu_size,
    const int &in_dofNum, const int &in_elemType,
    const bool isPrintInfo )
: nElem( mesh->get_nElem() ), nFunc( mesh->get_nFunc() ),
  sDegree( mesh->get_s_degree() ), tDegree( mesh->get_t_degree() ),
  uDegree( mesh->get_u_degree() ), nLocBas( mesh->get_nLocBas() ),
  probDim(3), dofNum( in_dofNum ), dofMat(in_dofNum), 
  elemType(in_elemType)
{
  // Initialize group 3 data
  cpu_rank = in_cpu_rank;
  cpu_size = in_cpu_size;
  isMETIS           = gpart->get_isMETIS();
  part_isdual       = gpart->get_isDual();
  dual_edge_ncommon = gpart->get_dual_edge_ncommon();

  // Check the cpu info
  SYS_T::print_exit_if(cpu_size < 1, "Error: Part_Tet input cpu_size is wrong! \n");
  SYS_T::print_exit_if(cpu_rank >= cpu_size, "Error: Part_Tet input cpu_rank is wrong! \n");
  SYS_T::print_exit_if(cpu_rank < 0, "Error: Part_Tet input cpu_rank is wrong! \n");

  // Generate group 1, 2, 5, and 6.
  Generate_Partition( mesh, gpart, mnindex, IEN, ctrlPts, isPrintInfo );
}


Part_Tet::Part_Tet(
    const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const int &in_cpu_rank, const int &in_cpu_size,
    const int &in_dofNum, const int &in_dofMat,
    const int &in_elemType, const bool isPrintInfo )
: nElem( mesh->get_nElem() ), nFunc( mesh->get_nFunc() ),
  sDegree( mesh->get_s_degree() ), tDegree( mesh->get_t_degree() ),
  uDegree( mesh->get_u_degree() ), nLocBas( mesh->get_nLocBas() ),
  probDim(3), dofNum( in_dofNum ), dofMat( in_dofMat ), 
  elemType(in_elemType)
{
  // Initialize group 3 data
  cpu_rank = in_cpu_rank;
  cpu_size = in_cpu_size;
  isMETIS           = gpart->get_isMETIS();
  part_isdual       = gpart->get_isDual();
  dual_edge_ncommon = gpart->get_dual_edge_ncommon();

  // Check the cpu info
  SYS_T::print_exit_if(cpu_size < 1, "Error: Part_Tet input cpu_size is wrong! \n");
  SYS_T::print_exit_if(cpu_rank >= cpu_size, "Error: Part_Tet input cpu_rank is wrong! \n");
  SYS_T::print_exit_if(cpu_rank < 0, "Error: Part_Tet input cpu_rank is wrong! \n");

  // Generate group 1, 2, 5, and 6.
  Generate_Partition( mesh, gpart, mnindex, IEN, ctrlPts, isPrintInfo );
}


Part_Tet::~Part_Tet()
{
  for(int i=0; i<nlocalele; ++i) delete [] LIEN[i];
  delete [] LIEN;
}


void Part_Tet::Generate_Partition( const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const bool &isPrintInfo )
{
  // 1. Create local partition based on the epart & npart info
  elem_loc.clear(); node_loc.clear();

  for( s_int e=0; e<nElem; ++e )
  {
    if( gpart->get_epart(e) == cpu_rank )
      elem_loc.push_back(e);
  }
  VEC_T::shrink2fit(elem_loc);
  nlocalele = (int) elem_loc.size();

  for( s_int n=0; n<nFunc; ++n )
  {
    if( gpart->get_npart(n) == cpu_rank )
    {
      node_loc.push_back(n);
      node_loc_original.push_back(n);
    }
  }
  VEC_T::shrink2fit(node_loc);
  VEC_T::shrink2fit(node_loc_original);
  nlocalnode = (int) node_loc.size();

  if(isPrintInfo)
  {
    std::cout<<"-- proc "<<cpu_rank<<" -- elem_loc & node_loc arrays generated. \n";
    std::cout<<"-- proc "<<cpu_rank<<" local element number: "<<elem_loc.size()<<std::endl;
  }

  // 2. Reorder node_loc
  for( int ii=0; ii<nlocalnode; ++ii )
    node_loc[ii] = mnindex->get_old2new( node_loc[ii] );

  // 3. Generate node_tot, which stores the nodes needed by the subdomain
  std::vector<s_int> node_tot;
  node_tot.clear();
  for( int e=0; e<nlocalele; ++e )
  {
    for( int ii=0; ii<nLocBas; ++ii )
    {
      s_int temp_node = IEN->get_IEN(elem_loc[e], ii);
      temp_node = mnindex->get_old2new(temp_node);
      node_tot.push_back( temp_node );
    }
  }
  sort(node_tot.begin(), node_tot.end());
  std::vector<s_int>::iterator it = unique(node_tot.begin(), node_tot.end());
  node_tot.resize( it - node_tot.begin() );

  ntotalnode = (int) node_tot.size();

  // 4. generate node_ghost
  node_ghost.clear();
  for( int ii = 0; ii<ntotalnode; ++ii )
  {
    it = find(node_loc.begin(), node_loc.end(), node_tot[ii]);
    if( it == node_loc.end() )
      node_ghost.push_back(node_tot[ii]);
  }

  VEC_T::shrink2fit(node_ghost);
  nghostnode = (int) node_ghost.size();

  nbadnode = 0;
  if( nghostnode + nlocalnode != ntotalnode )
  {
    std::vector<s_int> badnode;
    std::vector<s_int>::iterator badnode_it;
    for( int n=0; n<nlocalnode; ++n )
    {
      badnode_it = find( node_tot.begin(), node_tot.end(), node_loc[n] );
      if(badnode_it == node_tot.end())
        badnode.push_back( node_loc[n] );
    }
    nbadnode = (int) badnode.size();
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
  local_to_global.clear();
  for( s_int n=0; n<nlocalnode; ++n )
    local_to_global.push_back( node_loc[n] );
  for( s_int n=0; n<nghostnode; ++n )
    local_to_global.push_back( node_ghost[n] );

  VEC_T::shrink2fit(local_to_global);
  nlocghonode = (int) local_to_global.size();

  if( isPrintInfo )
    std::cout<<"-- proc "<<cpu_rank<<" local_to_global generated. \n";

  // 6. LIEN
  LIEN = new int * [nlocalele];
  for(int e=0; e<nlocalele; ++e)
    LIEN[e] = new int [nLocBas];

  std::vector<int>::iterator lien_ptr;
  s_int global_index;
  for(int e=0; e<nlocalele; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      global_index = IEN->get_IEN(elem_loc[e], i);
      global_index = mnindex->get_old2new(global_index);
      lien_ptr = find( local_to_global.begin(), local_to_global.end(),
          global_index);

      if(lien_ptr == local_to_global.end())
      {
        std::cerr<<"ERROR: Failed to generate LIEN array for "<<global_index<<std::endl;
        exit(EXIT_FAILURE);
      }

      LIEN[e][i] = lien_ptr - local_to_global.begin();
    }
  }

  if(isPrintInfo)
    std::cout<<"-- proc "<<cpu_rank<<" LIEN generated. \n";

  // 7. local copy of control points
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


int Part_Tet::get_elemLocIndex(const int &gloindex) const
{
  std::vector<int>::const_iterator findindex;
  findindex = find(elem_loc.begin(), elem_loc.end(), gloindex);
  if(findindex == elem_loc.end())
    return -1;
  else
    return (findindex - elem_loc.begin());
}


void Part_Tet::write( const char * inputFileName ) const
{
  std::string fName( inputFileName );
  fName.append("_p");

  if( cpu_rank / 10 == 0 )
    fName.append("0000");
  else if( cpu_rank / 100 == 0 )
    fName.append("000");
  else if( cpu_rank / 1000 == 0 )
    fName.append("00");
  else if( cpu_rank / 10000 == 0 )
    fName.append("0");

  std::stringstream sstrm;
  sstrm<<cpu_rank;
  fName.append(sstrm.str());

  fName.append(".h5");

  hid_t file_id;
  file_id = H5Fcreate(fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  // group 1: local element
  hid_t group_id_1 = H5Gcreate(file_id, "/Local_Elem", H5P_DEFAULT, 
      H5P_DEFAULT, H5P_DEFAULT);
  
  h5w->write_intScalar( group_id_1, "nlocalele", nlocalele );
  h5w->write_intVector( group_id_1, "elem_loc", elem_loc );
    
  H5Gclose( group_id_1 );

  // group 2: local node
  hid_t group_id_2 = H5Gcreate( file_id, "/Local_Node", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  
  h5w->write_intScalar( group_id_2, "nlocalnode", nlocalnode );
  h5w->write_intScalar( group_id_2, "nghostnode", nghostnode );
  h5w->write_intScalar( group_id_2, "ntotalnode", ntotalnode );
  h5w->write_intScalar( group_id_2, "nbadnode", nbadnode );
  h5w->write_intScalar( group_id_2, "nlocghonode", nlocghonode );

  h5w->write_intVector( group_id_2, "node_loc", node_loc );
  h5w->write_intVector( group_id_2, "node_loc_original", node_loc_original );
  h5w->write_intVector( group_id_2, "local_to_global", local_to_global );
  if(nghostnode > 0)
    h5w->write_intVector( group_id_2, "node_ghost", node_ghost );
  
  H5Gclose( group_id_2 );

  // group 3: global mesh info
  hid_t group_id_3 = H5Gcreate(file_id, "/Global_Mesh_Info", 
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  h5w->write_intScalar( group_id_3, "nElem", nElem );
  h5w->write_intScalar( group_id_3, "nFunc", nFunc );
  
  std::vector<int> vdeg; vdeg.clear();
  vdeg.push_back(sDegree); vdeg.push_back(tDegree); vdeg.push_back(uDegree);
  
  h5w->write_intVector( group_id_3, "degree", vdeg );

  h5w->write_intScalar( group_id_3, "nLocBas", nLocBas );
  
  h5w->write_intScalar( group_id_3, "probDim", probDim );
  h5w->write_intScalar( group_id_3, "dofNum", dofNum );
  h5w->write_intScalar( group_id_3, "dofMat", dofMat );
  h5w->write_intScalar( group_id_3, "elemType", elemType );

  H5Gclose( group_id_3 );

  // group 4: part info
  hid_t group_id_4 = H5Gcreate( file_id, "/Part_Info", H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT ); 
  
  h5w->write_intScalar( group_id_4, "cpu_rank", cpu_rank );
  h5w->write_intScalar( group_id_4, "cpu_size", cpu_size );
  h5w->write_intScalar( group_id_4, "dual_edge_ncommon", dual_edge_ncommon );

  H5Gclose( group_id_4 );

  // group 5: LIEN
  hid_t group_id_5 = H5Gcreate(file_id, "/LIEN", H5P_DEFAULT, 
      H5P_DEFAULT, H5P_DEFAULT);

  std::vector<int> row_LIEN;

  row_LIEN.resize(nlocalele * nLocBas);

  for(int e=0; e<nlocalele; ++e)
  {
    for(int ii=0; ii<nLocBas; ++ii)
      row_LIEN[e*nLocBas + ii] = LIEN[e][ii];
  }

  h5w -> write_intMatrix( group_id_5, "LIEN", row_LIEN, nlocalele, nLocBas);

  H5Gclose( group_id_5 );
  
  // group 6: control points
  hid_t group_id_6 = H5Gcreate(file_id, "/ctrlPts_loc", H5P_DEFAULT, 
      H5P_DEFAULT, H5P_DEFAULT);

  h5w -> write_doubleVector( group_id_6, "ctrlPts_x_loc", ctrlPts_x_loc );
  h5w -> write_doubleVector( group_id_6, "ctrlPts_y_loc", ctrlPts_y_loc );
  h5w -> write_doubleVector( group_id_6, "ctrlPts_z_loc", ctrlPts_z_loc );

  H5Gclose( group_id_6 );

  // Finish writing, clean up
  delete h5w;
  H5Fclose(file_id);
}


void Part_Tet::print_part_ele() const
{
  std::cout<<"Proc: "<<cpu_rank<<" local elements: "<<std::endl;
  for(int e=0; e<nlocalele; ++e)
    std::cout<<elem_loc[e]<<'\t';
  std::cout<<std::endl;
}


void Part_Tet::print_part_node() const
{
  std::cout<<"Proc: "<<cpu_rank<<" local nodes: "<<std::endl;
  for(int n=0; n<nlocalnode; ++n)
    std::cout<<n<<'\t'<<node_loc_original[n]<<'\t'<<node_loc[n]<<'\n';
  std::cout<<std::endl;

  for(int n=0; n<nlocghonode; ++n)
  {
    std::cout<<n<<'\t'<<ctrlPts_x_loc[n]<<'\t'<<ctrlPts_y_loc[n]<<'\t';
    std::cout<<ctrlPts_z_loc[n]<<'\n';
  }
}


void Part_Tet::print_part_ghost_node() const
{
  std::cout<<"Proc: "<<cpu_rank<<" ghost nodes: "<<std::endl;
  for(int n=0; n<nghostnode; ++n)
    std::cout<<node_ghost[n]<<"("<<nlocalnode + n<<")"<<'\t';
  std::cout<<std::endl;
}


void Part_Tet::print_part_local_to_global() const
{
  std::cout<<"Proc: "<<cpu_rank<<" local_to_global array: "<<std::endl;
  for(int n=0; n<nlocghonode; ++n)
    std::cout<<local_to_global[n]<<"("<<n<<")"<<'\t';
  std::cout<<std::endl;
}


void Part_Tet::print_part_LIEN() const
{
  std::cout<<"Proc: "<<cpu_rank<<" LIEN: "<<std::endl;
  for(int ee=0; ee<nlocalele; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
      std::cout<<LIEN[ee][ii]<<'\t';
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
}


void Part_Tet::print_part_loadbalance_edgecut() const
{
  std::cout<<"Proc:"<<" "<<cpu_rank;
  std::cout<<" "<<"element ratio:"<<" "<<(double) nlocalele / (double) nElem;
  std::cout<<" "<<"node ratio:"<<" "<<(double) nlocalnode / (double) nFunc;
  std::cout<<" "<<"gho/loc ratio:"<<" "<<(double) nghostnode / (double) nlocalnode<<std::endl;
}

// EOF
