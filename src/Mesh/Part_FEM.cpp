#include "Part_FEM.hpp"

Part_FEM::Part_FEM(
    const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const int &in_cpu_rank, const int &in_cpu_size,
    const int &in_elemType, const Field_Property &fp )
: nElem( mesh->get_nElem() ), nFunc( mesh->get_nFunc() ),
  sDegree( mesh->get_s_degree() ), tDegree( mesh->get_t_degree() ),
  uDegree( mesh->get_u_degree() ), nLocBas( mesh->get_nLocBas() ),
  probDim(3), elemType(in_elemType),
  field_id( fp.get_id() ), dofNum( fp.get_dofNum() ),
  is_geo_field( fp.get_is_geo_field() ),
  field_name( fp.get_name() )
{
  // Initialize group 3 data
  cpu_rank = in_cpu_rank;
  cpu_size = in_cpu_size;

  // Check the cpu info
  SYS_T::print_fatal_if(cpu_size < 1, "Error: Part_FEM input cpu_size is wrong! \n");
  SYS_T::print_fatal_if(cpu_rank >= cpu_size, "Error: Part_FEM input cpu_rank is wrong! \n");
  SYS_T::print_fatal_if(cpu_rank < 0, "Error: Part_FEM input cpu_rank is wrong! \n");

  // Generate group 1, 2, and 5.
  Generate_Partition( mesh, gpart, mnindex, IEN, field_id );

  // Generate group 6, if the field is tagged as is_geo_field == true
  // local copy of control points
  if( is_geo_field == true )
  {
    ctrlPts_x_loc.resize(nlocghonode);
    ctrlPts_y_loc.resize(nlocghonode);
    ctrlPts_z_loc.resize(nlocghonode);

    PERIGEE_OMP_PARALLEL_FOR
    for(int ii=0; ii<nlocghonode; ++ii)
    {
      int aux_index = local_to_global[ii];         // new global index
      aux_index = mnindex->get_new2old(aux_index); // back to old global index
      ctrlPts_x_loc[ii] = ctrlPts[3*aux_index + 0];
      ctrlPts_y_loc[ii] = ctrlPts[3*aux_index + 1];
      ctrlPts_z_loc[ii] = ctrlPts[3*aux_index + 2];
    }

    VEC_T::shrink2fit(ctrlPts_x_loc);
    VEC_T::shrink2fit(ctrlPts_y_loc);
    VEC_T::shrink2fit(ctrlPts_z_loc);

    std::cout<<"-- proc "<<cpu_rank<<" Local control points generated. \n";
  }
  else
  {
    ctrlPts_x_loc.clear();
    ctrlPts_y_loc.clear();
    ctrlPts_z_loc.clear();
  }
}

Part_FEM::Part_FEM( const std::string &inputfileName, const int &in_cpu_rank )
{
  const std::string fName = SYS_T::gen_partfile_name( inputfileName, in_cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  // local elements 
  elem_loc = h5r->read_intVector( "Local_Elem", "elem_loc" );
  nlocalele = h5r->read_intScalar( "Local_Elem", "nlocalele" );

  // local node
  nlocalnode  = h5r->read_intScalar("Local_Node", "nlocalnode");
  nghostnode  = h5r->read_intScalar("Local_Node", "nghostnode");
  nbadnode    = h5r->read_intScalar("Local_Node", "nbadnode");
  nlocghonode = h5r->read_intScalar("Local_Node", "nlocghonode");
  ntotalnode  = h5r->read_intScalar("Local_Node", "ntotalnode");

  local_to_global = h5r->read_intVector("Local_Node", "local_to_global");
  if( nghostnode > 0)
    node_ghost = h5r->read_intVector("Local_Node", "node_ghost");
  else
    node_ghost.clear();

  node_loc = h5r->read_intVector("Local_Node", "node_loc");
  node_loc_original = h5r->read_intVector("Local_Node", "node_loc_original");

  // Part info
  cpu_rank = h5r->read_intScalar("Part_Info", "cpu_rank");

  SYS_T::print_fatal_if( cpu_rank != in_cpu_rank, "Error: Part_FEM::cpu_rank is inconsistent.\n");

  cpu_size = h5r->read_intScalar("Part_Info", "cpu_size");

  // global mesh info
  std::vector<int> vdeg = h5r -> read_intVector("Global_Mesh_Info", "degree");

  sDegree = vdeg[0]; tDegree = vdeg[1]; uDegree = vdeg[2];

  nElem    = h5r -> read_intScalar("Global_Mesh_Info", "nElem");
  nFunc    = h5r -> read_intScalar("Global_Mesh_Info", "nFunc");
  nLocBas  = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");
  probDim  = h5r -> read_intScalar("Global_Mesh_Info", "probDim");
  elemType = h5r -> read_intScalar("Global_Mesh_Info", "elemType");
  dofNum   = h5r -> read_intScalar("Global_Mesh_Info", "dofNum");
  field_id = h5r -> read_intScalar("Global_Mesh_Info", "field_id");
  field_name = h5r -> read_string("Global_Mesh_Info", "field_name" );
  
  const int temp = h5r -> read_intScalar("Global_Mesh_Info", "is_geo_field");
  if(temp == 1) is_geo_field = true;
  else is_geo_field = false;

  // LIEN
  int num_row, num_col;
  const std::vector<int> LIEN_vec = h5r -> read_intMatrix("LIEN", "LIEN", num_row, num_col);

  SYS_T::print_fatal_if( num_row != nlocalele, "Error: Part_FEM::LIEN size does not match the number of element. \n");

  SYS_T::print_fatal_if( num_col != nLocBas, "Error: Part_FEM::LIEN size does not match the value of nLocBas. \n");

  LIEN = new int * [nlocalele];
  for(int ee=0; ee<nlocalele; ++ee) LIEN[ee] = new int [nLocBas];

  for(int ee=0; ee<nlocalele; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii) LIEN[ee][ii] = LIEN_vec[ee*nLocBas + ii];
  }

  // control points
  ctrlPts_x_loc = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_x_loc");
  ctrlPts_y_loc = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_y_loc");
  ctrlPts_z_loc = h5r -> read_doubleVector("ctrlPts_loc", "ctrlPts_z_loc");

  delete h5r; H5Fclose( file_id );
}

Part_FEM::~Part_FEM()
{
  for(int ii=0; ii<nlocalele; ++ii) delete [] LIEN[ii];
  delete [] LIEN;
}

void Part_FEM::Generate_Partition( const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const int &field )
{
  // 1. Create local partition based on the epart & npart info
  elem_loc.clear(); node_loc.clear();

  for( int e=0; e<nElem; ++e )
  {
    if( gpart->get_epart(e) == cpu_rank ) elem_loc.push_back(e);
  }
  VEC_T::shrink2fit(elem_loc);
  nlocalele = VEC_T::get_size( elem_loc );

  for( int n=0; n<nFunc; ++n )
  {
    if( gpart->get_npart(n, field) == cpu_rank )
    {
      node_loc.push_back(n);
      node_loc_original.push_back(n);
    }
  }
  VEC_T::shrink2fit(node_loc); VEC_T::shrink2fit(node_loc_original);
  nlocalnode = VEC_T::get_size( node_loc );

  std::cout<<"-- proc "<<cpu_rank<<" -- elem_loc & node_loc arrays generated. \n";
  std::cout<<"-- proc "<<cpu_rank<<" local element number: "<<elem_loc.size()<<std::endl;

  // 2. Reorder node_loc
  PERIGEE_OMP_PARALLEL_FOR
  for( int ii=0; ii<nlocalnode; ++ii ) 
    node_loc[ii] = mnindex->get_old2new( node_loc[ii] );

  // 3. Generate node_tot, which stores the nodes needed by the elements in the subdomain
  std::vector<int> node_tot {};
  PERIGEE_OMP_PARALLEL
  {
    std::vector<int> temp_node_tot {};
    PERIGEE_OMP_FOR
    for( int e=0; e<nlocalele; ++e )
    {
      for( int ii=0; ii<nLocBas; ++ii )
      {
        int temp_node = IEN->get_IEN(elem_loc[e], ii);
        temp_node = mnindex->get_old2new(temp_node);
        temp_node_tot.push_back( temp_node );
      }
    }
    PERIGEE_OMP_CRITICAL
    VEC_T::insert_end(node_tot, temp_node_tot);
  }
  VEC_T::sort_unique_resize( node_tot );

  ntotalnode = VEC_T::get_size( node_tot );

  // 4. generate node_ghost
  node_ghost.clear();
  for( int ii = 0; ii<ntotalnode; ++ii )
  {
    if( !VEC_T::is_invec(node_loc, node_tot[ii]) )
      node_ghost.push_back(node_tot[ii]);
  }

  VEC_T::shrink2fit(node_ghost);
  nghostnode = VEC_T::get_size( node_ghost );

  nbadnode = 0;
  if( nghostnode + nlocalnode != ntotalnode )
  {
    std::vector<int> badnode {};
    for( int n=0; n<nlocalnode; ++n )
    {
      if( !VEC_T::is_invec(node_tot, node_loc[n]) )
        badnode.push_back( node_loc[n] );
    }
    nbadnode = VEC_T::get_size( badnode );
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

  std::cout<<"-- proc "<<cpu_rank;
  std::cout<<" -- ntotalnode: "<<ntotalnode;
  std::cout<<" -- nlocalnode: "<<nlocalnode;
  std::cout<<" -- nghostnode: "<<nghostnode;
  std::cout<<" -- nbadnode: "<<nbadnode<<std::endl;

  // 5. local_to_global mapping
  local_to_global.clear();
  for( int n=0; n<nlocalnode; ++n )
    local_to_global.push_back( node_loc[n] );
  for( int n=0; n<nghostnode; ++n )
    local_to_global.push_back( node_ghost[n] );

  VEC_T::shrink2fit(local_to_global);
  nlocghonode = VEC_T::get_size( local_to_global );

  std::cout<<"-- proc "<<cpu_rank<<" local_to_global generated. \n";

  // 6. LIEN
  LIEN = new int * [nlocalele];
  for(int ee=0; ee<nlocalele; ++ee) LIEN[ee] = new int [nLocBas];

  PERIGEE_OMP_PARALLEL_FOR
  for(int ee=0; ee<nlocalele; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int global_index = mnindex->get_old2new( IEN->get_IEN(elem_loc[ee], ii) );
      const auto lien_ptr = find( local_to_global.begin(), local_to_global.end(), global_index );

      if(lien_ptr == local_to_global.end())
      {
        std::cerr<<"ERROR: Failed to generate LIEN array for "<<global_index<<std::endl;
        exit(EXIT_FAILURE);
      }

      LIEN[ee][ii] = lien_ptr - local_to_global.begin();
    }
  }

  std::cout<<"-- proc "<<cpu_rank<<" LIEN generated. \n";
}

void Part_FEM::write( const std::string &inputFileName ) const
{
  const std::string fName = SYS_T::gen_partfile_name( inputFileName, cpu_rank );

  hid_t file_id = H5Fcreate(fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  // group 1: local element
  hid_t group_id_1 = H5Gcreate(file_id, "/Local_Elem", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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
  hid_t group_id_3 = H5Gcreate(file_id, "/Global_Mesh_Info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  h5w->write_intScalar( group_id_3, "nElem", nElem );
  h5w->write_intScalar( group_id_3, "nFunc", nFunc );

  const std::vector<int> vdeg { sDegree, tDegree, uDegree };

  h5w->write_intVector( group_id_3, "degree", vdeg );

  h5w->write_intScalar( group_id_3, "nLocBas", nLocBas );

  h5w->write_intScalar( group_id_3, "probDim", probDim );
  h5w->write_intScalar( group_id_3, "dofNum", dofNum );
  h5w->write_intScalar( group_id_3, "elemType", elemType );

  h5w->write_intScalar( group_id_3, "field_id", field_id );
  h5w->write_intScalar( group_id_3, "is_geo_field", (is_geo_field ? 1 : 0) );
  h5w->write_string( group_id_3, "field_name", field_name );

  H5Gclose( group_id_3 );

  // group 4: part info
  hid_t group_id_4 = H5Gcreate( file_id, "/Part_Info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ); 

  h5w->write_intScalar( group_id_4, "cpu_rank", cpu_rank );
  h5w->write_intScalar( group_id_4, "cpu_size", cpu_size );

  H5Gclose( group_id_4 );

  // group 5: LIEN
  hid_t group_id_5 = H5Gcreate(file_id, "/LIEN", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  std::vector<int> row_LIEN(nlocalele * nLocBas, -1);

  for(int e=0; e<nlocalele; ++e)
  {
    for(int ii=0; ii<nLocBas; ++ii)
      row_LIEN[e*nLocBas + ii] = LIEN[e][ii];
  }

  h5w -> write_intMatrix( group_id_5, "LIEN", row_LIEN, nlocalele, nLocBas);

  H5Gclose( group_id_5 );

  // group 6: control points
  if( is_geo_field == true )
  {
    hid_t group_id_6 = H5Gcreate(file_id, "/ctrlPts_loc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    h5w -> write_doubleVector( group_id_6, "ctrlPts_x_loc", ctrlPts_x_loc );
    h5w -> write_doubleVector( group_id_6, "ctrlPts_y_loc", ctrlPts_y_loc );
    h5w -> write_doubleVector( group_id_6, "ctrlPts_z_loc", ctrlPts_z_loc );

    H5Gclose( group_id_6 );
  }

  // Finish writing, clean up
  delete h5w;
  H5Fclose(file_id);
}

void Part_FEM::print_part_ele() const
{
  std::cout<<"Proc: "<<cpu_rank<<" local elements: "<<std::endl;
  for(int e=0; e<nlocalele; ++e)
    std::cout<<elem_loc[e]<<'\t';
  std::cout<<std::endl;
}

void Part_FEM::print_part_node() const
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

void Part_FEM::print_part_ghost_node() const
{
  std::cout<<"Proc: "<<cpu_rank<<" ghost nodes: "<<std::endl;
  for(int n=0; n<nghostnode; ++n)
    std::cout<<node_ghost[n]<<"("<<nlocalnode + n<<")"<<'\t';
  std::cout<<std::endl;
}

void Part_FEM::print_part_local_to_global() const
{
  std::cout<<"Proc: "<<cpu_rank<<" local_to_global array: "<<std::endl;
  for(int n=0; n<nlocghonode; ++n)
    std::cout<<local_to_global[n]<<"("<<n<<")"<<'\t';
  std::cout<<std::endl;
}

void Part_FEM::print_part_LIEN() const
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

void Part_FEM::print_part_loadbalance_edgecut() const
{
  std::cout<<"Proc:"<<" "<<cpu_rank;
  std::cout<<" "<<"element ratio:"<<" "<<(double) nlocalele / (double) nElem;
  std::cout<<" "<<"node ratio:"<<" "<<(double) nlocalnode / (double) nFunc;
  std::cout<<" "<<"gho/loc ratio:"<<" "<<(double) nghostnode / (double) nlocalnode<<std::endl;
}

// EOF
