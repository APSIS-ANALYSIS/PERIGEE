#include "Part_NURBS_1Patch_3D_METIS.hpp"

Part_NURBS_1Patch_3D_METIS::Part_NURBS_1Patch_3D_METIS( 
    const class IMesh * const &mesh,
    const class IGlobal_Part * const &gpart,
    const class Map_Node_Index * const &mnindex,
    const class IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const std::vector<NURBS_T::BezierElem*> &in_seg_x,
    const std::vector<NURBS_T::BezierElem*> &in_seg_y,
    const std::vector<NURBS_T::BezierElem*> &in_seg_z,
    const int &in_cpu_rank, const int &in_cpu_size,
    const int &in_dofNum, const int &in_elemType,
    const bool isPrintInfo )
{
  cpu_rank = in_cpu_rank;
  cpu_size = in_cpu_size;

  dofNum   = in_dofNum;
  elemType = in_elemType;
  probDim  = 3;

  isMETIS = gpart->get_isMETIS();
  part_isdual = gpart->get_isDual();
  dual_edge_ncommon = gpart->get_dual_edge_ncommon();

  if( cpu_size < 1 )
  {
    std::cerr<<"ERROR: given cpu_size "<<cpu_size<<" is wrong. \n";
    exit(1);
  }
  if( cpu_rank >= cpu_size )
  {
    std::cerr<<"ERROR: cpu_rank "<<cpu_rank<<" is greater than cpu_size "<<cpu_size<<std::endl;
    exit(1);
  }
  
  // Get basic global mesh info: degree, nFunc, nElem, nLocBas, hx_min,
  // hx_max.
  Get_global_meshinfo(mesh);

  // Create local partition info
  Generate_Partition(mesh, gpart, mnindex, IEN, ctrlPts, isPrintInfo);

  // store the pointers to seg_xyz
  pseg_x.clear(); pseg_x = in_seg_x;
  pseg_y.clear(); pseg_y = in_seg_y;
  pseg_z.clear(); pseg_z = in_seg_z;
  
  // Create indices for nonzero element
  nonzero_x.clear(); nonzero_y.clear(); nonzero_z.clear();
  int temp_index = 0;
  for(int ee = 0 ; ee<nElem_x; ++ee)
  {
    if(mesh->get_hx(ee) == 0.0)
      nonzero_x.push_back(-1);
    else
    {
      nonzero_x.push_back(temp_index);
      temp_index += 1;
    }
  }
  temp_index = 0;
  for(int ee=0; ee<nElem_y; ++ee)
  {
    if(mesh->get_hy(ee*nElem_x) == 0.0)
      nonzero_y.push_back(-1);
    else
    {
      nonzero_y.push_back(temp_index);
      temp_index += 1;
    }
  }
  temp_index = 0;
  for(int ee=0; ee<nElem_z; ++ee)
  {
    if(mesh->get_hz(ee*nElem_x*nElem_y) == 0.0)
      nonzero_z.push_back(-1);
    else
    {
      nonzero_z.push_back(temp_index);
      temp_index += 1;
    }
  }
  VEC_T::shrink2fit(nonzero_x);
  VEC_T::shrink2fit(nonzero_y);
  VEC_T::shrink2fit(nonzero_z);
}

Part_NURBS_1Patch_3D_METIS::~Part_NURBS_1Patch_3D_METIS()
{
  for(int i=0; i<nlocalele; ++i)
    delete [] LIEN[i];
  delete [] LIEN;
}

void Part_NURBS_1Patch_3D_METIS::Generate_Partition( const class IMesh * const &mesh,
    const class IGlobal_Part * const &gpart,
    const class Map_Node_Index * const &mnindex, 
    const class IIEN * const &IEN,
    const std::vector<double> &ctrlPts, const bool isPrintInfo )
{

  // 1. Create local partition based on the epart npart info
  elem_loc.clear(); node_loc.clear();

  for( int e=0; e<nElem; ++e )
  {
    if( gpart->get_epart(e) == cpu_rank )
      elem_loc.push_back(e);
  }
  VEC_T::shrink2fit(elem_loc);
  nlocalele = elem_loc.size();

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
  nlocalnode = node_loc.size();

  if(isPrintInfo)
  {
    std::cout<<"-- proc "<<cpu_rank<<" -- elem_loc & node_loc arrays generated. \n";
    std::cout<<"-- proc "<<cpu_rank<<" local element number: "<<elem_loc.size()<<std::endl;
  }
  // 2. reorder node_loc
  for( int ii=0; ii<nlocalnode; ++ii )
    node_loc[ii] = mnindex->get_old2new( node_loc[ii] );

  // 3. generate node_tot, which stores the nodes needed by the subdomain
  std::vector<int> node_tot;
  node_tot.clear();
  for( int e=0; e<nlocalele; ++e )
  {
    for( int ii=0; ii<nLocBas; ++ii )
    {
      int temp_node = IEN->get_IEN(elem_loc[e], ii);
      temp_node = mnindex->get_old2new(temp_node);
      node_tot.push_back( temp_node );
    }
  }
  sort(node_tot.begin(), node_tot.end());
  std::vector<int>::iterator it = unique(node_tot.begin(), node_tot.end());
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
    std::vector<int> badnode;
    std::vector<int>::iterator badnode_it;
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
      exit(1);
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
  for( int n=0; n<nlocalnode; ++n )
    local_to_global.push_back( node_loc[n] );
  for( int n=0; n<nghostnode; ++n )
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
  int global_index;
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
        exit(1);
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
  ctrlPts_w_loc.resize(nlocghonode);

  for(int ii=0; ii<nlocghonode; ++ii)
  {
    int aux_index = local_to_global[ii]; // new global index
    aux_index = mnindex->get_new2old(aux_index); // back to old global index
    ctrlPts_x_loc[ii] = ctrlPts[4*aux_index + 0];
    ctrlPts_y_loc[ii] = ctrlPts[4*aux_index + 1];
    ctrlPts_z_loc[ii] = ctrlPts[4*aux_index + 2];
    ctrlPts_w_loc[ii] = ctrlPts[4*aux_index + 3];
  }
 
  VEC_T::shrink2fit(ctrlPts_x_loc);
  VEC_T::shrink2fit(ctrlPts_y_loc);
  VEC_T::shrink2fit(ctrlPts_z_loc);
  VEC_T::shrink2fit(ctrlPts_w_loc);
  
  if(isPrintInfo)
    std::cout<<"-- proc "<<cpu_rank<<" Local control points generated. \n";

  // 8. hx hy hz for local elements
  for( int e=0; e<nlocalele; ++e )
  {
    int e_global = get_elem_loc(e);
    hx.push_back(mesh->get_hx(e_global));
    hy.push_back(mesh->get_hy(e_global));
    hz.push_back(mesh->get_hz(e_global));
  }
}

void Part_NURBS_1Patch_3D_METIS::Get_global_meshinfo( 
    const class IMesh * const &mesh )
{
  // get global mesh basic data
  sDegree = mesh->get_s_degree();
  tDegree = mesh->get_t_degree();
  uDegree = mesh->get_u_degree();

  nFunc = mesh->get_nFunc();
  nFunc_x = mesh->get_nFunc_x();
  nFunc_y = mesh->get_nFunc_y();
  nFunc_z = mesh->get_nFunc_z();

  nElem = mesh->get_nElem();
  nElem_x = mesh->get_nElem_x();
  nElem_y = mesh->get_nElem_y();
  nElem_z = mesh->get_nElem_z();

  nLocBas = mesh->get_nLocBas();

  hx_max = mesh->get_hx_max();
  hy_max = mesh->get_hy_max();
  hz_max = mesh->get_hz_max();

  hx_min = mesh->get_hx_min();
  hy_min = mesh->get_hy_min();
  hz_min = mesh->get_hz_min();
}

void Part_NURBS_1Patch_3D_METIS::print_info() const
{
  std::cout<<"=================================="<<std::endl;
  std::cout<<"isMETIS: "<<get_isMETIS()<<std::endl;
  std::cout<<"isDual: "<<get_part_isDual()<<std::endl;
  std::cout<<"dual_dege_ncommon: "<<get_dual_edge_ncommon()<<std::endl;
  std::cout<<"=================================="<<std::endl;
}

void Part_NURBS_1Patch_3D_METIS::print_part_ele() const
{
  std::cout<<"Proc: "<<cpu_rank<<" local elements: "<<std::endl;
  for(int e=0; e<nlocalele; ++e)
    std::cout<<get_elem_loc(e)<<'\t';
  std::cout<<std::endl;
}

void Part_NURBS_1Patch_3D_METIS::print_part_node() const
{
  std::cout<<"Proc: "<<cpu_rank<<" local nodes: "<<std::endl;
  for(int n=0; n<nlocalnode; ++n)
  {
    std::cout<<n<<'\t'<<node_loc_original[n]<<'\t'<<node_loc[n]<<'\t';
  }
  std::cout<<std::endl;
  
  for(int n=0; n<nlocghonode; ++n)
  {
    std::cout<<n<<'\t'<<ctrlPts_x_loc[n]<<'\t'<<ctrlPts_y_loc[n]<<'\t';
    std::cout<<ctrlPts_z_loc[n]<<'\t'<<ctrlPts_w_loc[n]<<'\n';
  }
}

void Part_NURBS_1Patch_3D_METIS::print_part_ghost_node() const
{
  std::cout<<"Proc: "<<cpu_rank<<" ghost nodes: "<<std::endl;
  for(int n=0; n<nghostnode; ++n)
    std::cout<<node_ghost[n]<<"("<<nlocalnode + n<<")"<<'\t';
  std::cout<<std::endl;
}

void Part_NURBS_1Patch_3D_METIS::print_part_local_to_global() const
{
  std::cout<<"Proc: "<<cpu_rank<<" local_to_global array: "<<std::endl;
  for(int n=0; n<nlocghonode; ++n)
    std::cout<<local_to_global[n]<<"("<<n<<")"<<'\t';
  std::cout<<std::endl;
}

void Part_NURBS_1Patch_3D_METIS::print_part_LIEN() const
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

void Part_NURBS_1Patch_3D_METIS::print_part_loadbalance_edgecut() const
{
  std::cout<<"Proc:"<<" "<<cpu_rank;
  std::cout<<" "<<"element ratio:"<<" "<<(double) nlocalele / (double) nElem;
  std::cout<<" "<<"node ratio:"<<" "<<(double) nlocalnode / (double) nFunc;
  std::cout<<" "<<"gho/loc ratio:"<<" "<<(double) nghostnode / (double) nlocalnode<<std::endl;
}

void Part_NURBS_1Patch_3D_METIS::print_part_bezier_ext_x(int elem) const
{
  int elem_x, elem_y, elem_z;
  SYS_T::get_xyz_index(elem, nElem_x, nElem_y, elem_x, elem_y, elem_z);
  
  std::vector<double> extractor;
  Get_bezier_ext_x(elem_x, extractor);

  VEC_T::print(extractor);
}

void Part_NURBS_1Patch_3D_METIS::print_part_bezier_ext_y(int elem) const
{
  int elem_x, elem_y, elem_z;
  SYS_T::get_xyz_index(elem, nElem_x, nElem_y, elem_x, elem_y, elem_z);
  
  std::vector<double> extractor;
  Get_bezier_ext_y(elem_y, extractor);

  VEC_T::print(extractor);
}

void Part_NURBS_1Patch_3D_METIS::print_part_bezier_ext_z(int elem) const
{
  int elem_x, elem_y, elem_z;
  SYS_T::get_xyz_index(elem, nElem_x, nElem_y, elem_x, elem_y, elem_z);
  
  std::vector<double> extractor;
  Get_bezier_ext_z(elem_z, extractor);

  VEC_T::print(extractor);
}

void Part_NURBS_1Patch_3D_METIS::Get_bezier_ext_x(const int &elem_x, std::vector<double> &ext_x) const
{
  ext_x.clear(); ext_x.resize((sDegree+1)*(sDegree+1));

  int nonzero_index = nonzero_x[elem_x];

  if(nonzero_index == -1)
  {
    for(int ii=0; ii<(sDegree+1)*(sDegree+1); ++ii)
      ext_x[ii] = 0.0;
  }
  else
  {
    for(int ii=0; ii<sDegree+1; ++ii)
    {
      for(int jj=0; jj<sDegree+1; ++jj)
        ext_x[ii*(sDegree+1) + jj] = pseg_x[nonzero_index]->coefficient(ii,jj);
    }
  }
}

void Part_NURBS_1Patch_3D_METIS::Get_bezier_ext_y(const int &elem_y, std::vector<double> &ext_y) const
{
  ext_y.clear(); ext_y.resize((tDegree+1)*(tDegree+1));
  int nonzero_index = nonzero_y[elem_y];

  if(nonzero_index == -1)
  {
    for(int ii=0; ii<(tDegree+1)*(tDegree+1); ++ii)
      ext_y[ii] = 0.0;
  }
  else
  {
    for(int ii=0; ii<tDegree+1; ++ii)
    {
      for(int jj=0; jj<tDegree+1; ++jj)
        ext_y[ii*(tDegree+1) +jj] = pseg_y[nonzero_index]->coefficient(ii,jj);
    }
  }
}

void Part_NURBS_1Patch_3D_METIS::Get_bezier_ext_z(const int &elem_z, std::vector<double> &ext_z) const
{
  ext_z.clear(); ext_z.resize((uDegree+1)*(uDegree+1));
  int nonzero_index = nonzero_z[elem_z];

  if(nonzero_index == -1)
  {
    for(int ii=0; ii<(uDegree+1)*(uDegree+1); ++ii)
      ext_z[ii] = 0.0;
  }
  else
  {
    for(int ii=0; ii<uDegree+1; ++ii)
    {
      for(int jj=0; jj<uDegree+1; ++jj)
        ext_z[ii*(uDegree+1) + jj] = pseg_z[nonzero_index]->coefficient(ii,jj);
    }
  }
}


void Part_NURBS_1Patch_3D_METIS::write( const char * inputFileName ) const
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

  // 1. local element
  write_local_elem(file_id);

  // 2. local node
  write_local_node(file_id);

  // 3. part info
  write_part_info(file_id);
  
  // 4. global mesh info  
  write_global_mesh_info(file_id);

  // 5. LIEN
  write_LIEN(file_id);
  
  // 6. local control points
  write_ctrlPts_loc(file_id);

  // 7. extraction operator
  write_extractor( file_id );

  // 8. hx hy hz
  write_hxyz( file_id );

  H5Fclose(file_id);
}

void Part_NURBS_1Patch_3D_METIS::write_local_elem( hid_t file_id ) const
{
  hid_t group_id, dataspace_id_1, dataspace_id_n, setid_nlocalele, setid_elem_loc;
  hsize_t dim1[1], dim2[1];
 
  // dataspace setup 
  dim1[0] = 1;
  dim2[0] = nlocalele;

  dataspace_id_1 = H5Screate_simple(1, dim1, NULL);
  dataspace_id_n = H5Screate_simple(1, dim2, NULL);

  // Create group
  group_id = H5Gcreate( file_id, "/Local_Elem", H5P_DEFAULT, H5P_DEFAULT,
     H5P_DEFAULT );
  
  setid_nlocalele = H5Dcreate( group_id, "nlocalele", H5T_NATIVE_INT,
      dataspace_id_1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_nlocalele, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nlocalele);
  H5Dclose( setid_nlocalele );

  setid_elem_loc = H5Dcreate( group_id, "elem_loc", H5T_NATIVE_INT,
      dataspace_id_n, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_elem_loc, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
      &elem_loc[0]);
  H5Dclose( setid_elem_loc );

  H5Sclose( dataspace_id_1 ); 
  H5Sclose( dataspace_id_n ); 
  H5Gclose( group_id );
}
  
void Part_NURBS_1Patch_3D_METIS::write_local_node( hid_t file_id ) const
{
  hid_t group_id, dataspace_id_1, dataspace_id_node_loc;
  hid_t dataspace_id_local_to_global;
  hsize_t dim1[1], dim2[1], dim4[1];

  dim1[0] = 1; dim2[0] = nlocalnode;
  dim4[0] = nlocghonode;
  dataspace_id_1               = H5Screate_simple(1, dim1, NULL);
  dataspace_id_node_loc        = H5Screate_simple(1, dim2, NULL);
  dataspace_id_local_to_global = H5Screate_simple(1, dim4, NULL);

  // Create group
  group_id = H5Gcreate( file_id, "/Local_Node", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  // scalar data
  hid_t setid_nlocalnode, setid_nghostnode, setid_ntotalnode, setid_nbadnode;
  hid_t setid_nlocghonode;

  setid_nlocalnode = H5Dcreate( group_id, "nlocalnode", H5T_NATIVE_INT,
     dataspace_id_1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_nlocalnode, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nlocalnode);
  H5Dclose( setid_nlocalnode );

  setid_nghostnode = H5Dcreate( group_id, "nghostnode", H5T_NATIVE_INT,
     dataspace_id_1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_nghostnode, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nghostnode);
  H5Dclose( setid_nghostnode );
  
  setid_ntotalnode = H5Dcreate( group_id, "ntotalnode", H5T_NATIVE_INT,
     dataspace_id_1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_ntotalnode, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ntotalnode);
  H5Dclose( setid_ntotalnode );

  setid_nbadnode = H5Dcreate( group_id, "nbadnode", H5T_NATIVE_INT,
     dataspace_id_1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_nbadnode, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nbadnode);
  H5Dclose( setid_nbadnode );

  setid_nlocghonode = H5Dcreate( group_id, "nlocghonode", H5T_NATIVE_INT,
     dataspace_id_1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_nlocghonode, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
      &nlocghonode);
  H5Dclose( setid_nlocghonode );

  // array data
  hid_t setid_node_loc, setid_node_loc_ori, setid_local_to_global;

  setid_node_loc = H5Dcreate( group_id, "node_loc", H5T_NATIVE_INT, 
     dataspace_id_node_loc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_node_loc, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &node_loc[0] );
  H5Dclose( setid_node_loc );

  setid_node_loc_ori = H5Dcreate( group_id, "node_loc_original", H5T_NATIVE_INT, 
     dataspace_id_node_loc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_node_loc_ori, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
      H5P_DEFAULT, &node_loc_original[0] );
  H5Dclose( setid_node_loc_ori );


  // if nghostnode > 0, write in node_ghost
  if(nghostnode > 0)
  {
    hid_t setid_node_ghost, dataspace_id_node_ghost;
    hsize_t dim3[1];
    dim3[0] = nghostnode;
    dataspace_id_node_ghost = H5Screate_simple(1, dim3, NULL);
    setid_node_ghost = H5Dcreate( group_id, "node_ghost", H5T_NATIVE_INT, 
        dataspace_id_node_ghost, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( setid_node_ghost, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
        H5P_DEFAULT, &node_ghost[0] );
    H5Dclose( setid_node_ghost );
    H5Sclose( dataspace_id_node_ghost );
  }


  setid_local_to_global = H5Dcreate( group_id, "local_to_global", H5T_NATIVE_INT, 
      dataspace_id_local_to_global, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_local_to_global, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
      H5P_DEFAULT, &local_to_global[0] );
  H5Dclose( setid_local_to_global );

  H5Sclose( dataspace_id_1 );
  H5Sclose( dataspace_id_node_loc );
  H5Sclose( dataspace_id_local_to_global );
  H5Gclose( group_id );
}

void Part_NURBS_1Patch_3D_METIS::write_part_info( hid_t file_id ) const
{
  hid_t group_id, dataspace_id;
  hsize_t dim[1];

  dim[0] = 1;
  dataspace_id = H5Screate_simple(1, dim, NULL);

  group_id = H5Gcreate( file_id, "/Part_Info", H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT );

  // int type
  hid_t setid_cpu_rank, setid_cpu_size, setid_dual_edge_ncommon; 

  setid_cpu_rank = H5Dcreate( group_id, "cpu_rank", H5T_NATIVE_INT,
      dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_cpu_rank, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &cpu_rank );
  H5Dclose( setid_cpu_rank );

  setid_cpu_size = H5Dcreate( group_id, "cpu_size", H5T_NATIVE_INT,
      dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_cpu_size, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &cpu_size );
  H5Dclose( setid_cpu_size );

  setid_dual_edge_ncommon = H5Dcreate( group_id, "dual_edge_ncommon", H5T_NATIVE_INT,
      dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_dual_edge_ncommon, H5T_NATIVE_INT, H5S_ALL, 
      H5S_ALL, H5P_DEFAULT, &dual_edge_ncommon );
  H5Dclose( setid_dual_edge_ncommon );

  H5Sclose( dataspace_id ); 
  H5Gclose( group_id );
}

void Part_NURBS_1Patch_3D_METIS::write_global_mesh_info( hid_t file_id ) const
{
  hid_t group_id, dataspace_id_4, dataspace_id_3, dataspace_id_1;
  hid_t setid_nElem, setid_nFunc, setid_degree, setid_hmax, setid_hmin, setid_lbas;
  hsize_t dim1[1], dim2[1], dim3[1];

  // Create group
  group_id = H5Gcreate(file_id, "/Global_Mesh_Info", H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT);

  // scalar type
  dim1[0] = 4;
  dataspace_id_4 = H5Screate_simple(1,dim1, NULL);

  dim2[0] = 3;
  dataspace_id_3 = H5Screate_simple(1, dim2, NULL);

  dim3[0] = 1;
  dataspace_id_1 = H5Screate_simple(1, dim3, NULL);

  // element dataset crate and write
  setid_nElem = H5Dcreate( group_id, "nElem", H5T_NATIVE_INT, dataspace_id_4,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  std::vector<int> vec_temp;
  vec_temp.clear();
  vec_temp.push_back(nElem); vec_temp.push_back(nElem_x);
  vec_temp.push_back(nElem_y); vec_temp.push_back(nElem_z);
  H5Dwrite( setid_nElem, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec_temp[0] );
  H5Dclose( setid_nElem );

  // function number dataset create and write
  setid_nFunc = H5Dcreate( group_id, "nFunc", H5T_NATIVE_INT, dataspace_id_4,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  vec_temp.clear();
  vec_temp.push_back(nFunc); vec_temp.push_back(nFunc_x);
  vec_temp.push_back(nFunc_y); vec_temp.push_back(nFunc_z);
  H5Dwrite( setid_nFunc, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec_temp[0] );
  H5Dclose( setid_nFunc );

  // degree dataset
  setid_degree = H5Dcreate( group_id, "degree", H5T_NATIVE_INT, dataspace_id_3,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  vec_temp.clear();
  vec_temp.push_back(sDegree); vec_temp.push_back(tDegree); vec_temp.push_back(uDegree);
  H5Dwrite( setid_degree, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec_temp[0] );
  H5Dclose( setid_degree );

  // hx_max
  setid_hmax = H5Dcreate( group_id, "h_max", H5T_NATIVE_DOUBLE, dataspace_id_3,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  std::vector<double> vec_temp2;
  vec_temp2.clear();
  vec_temp2.push_back(hx_max); vec_temp2.push_back(hy_max); vec_temp2.push_back(hz_max);
  H5Dwrite( setid_hmax, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec_temp2[0] );
  H5Dclose( setid_hmax );

  // hx_min
  setid_hmin = H5Dcreate( group_id, "h_min", H5T_NATIVE_DOUBLE, dataspace_id_3,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  vec_temp2.clear();
  vec_temp2.push_back(hx_min); vec_temp2.push_back(hy_min); vec_temp2.push_back(hz_min);
  H5Dwrite( setid_hmin, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec_temp2[0] );
  H5Dclose( setid_hmin );

  // nLocBas
  setid_lbas = H5Dcreate( group_id, "nLocBas", H5T_NATIVE_INT, dataspace_id_1,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( setid_lbas, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nLocBas );
  H5Dclose( setid_lbas );
  
  // probDim & dofNum
  HDF5_Writer * hdf5writer = new HDF5_Writer(file_id);
  hdf5writer->write_intScalar(group_id, "probDim", probDim);
  hdf5writer->write_intScalar(group_id, "dofNum", dofNum);

  delete hdf5writer;
  H5Sclose( dataspace_id_4 );
  H5Sclose( dataspace_id_3 );
  H5Sclose( dataspace_id_1 );
  H5Gclose( group_id );
}

void Part_NURBS_1Patch_3D_METIS::write_LIEN( hid_t file_id ) const
{
  hid_t group_id, dataspace_id, setid;
  hsize_t dim[2];

  group_id = H5Gcreate(file_id, "/LIEN", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // dim dataspace
  dim[0] = nlocalele; dim[1] = nLocBas;

  dataspace_id = H5Screate_simple(2, dim, NULL);

  // dataset
  setid = H5Dcreate( group_id, "LIEN", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT );

  int * row_LIEN = new int [nlocalele * nLocBas];
  for(int e=0; e<nlocalele; ++e)
  {
    for(int ii=0; ii<nLocBas; ++ii)
      row_LIEN[e*nLocBas + ii] = LIEN[e][ii];
  } 

  H5Dwrite( setid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, row_LIEN );

  delete [] row_LIEN; 

  H5Dclose( setid );
  H5Sclose( dataspace_id );
  H5Gclose( group_id );
}

void Part_NURBS_1Patch_3D_METIS::write_ctrlPts_loc( hid_t file_id ) const
{
  hid_t group_id, dataspace_id, setid_x, setid_y, setid_z, setid_w;
  hsize_t dim[1];

  dim[0] = nlocghonode;
  dataspace_id = H5Screate_simple(1, dim, NULL);

  group_id = H5Gcreate(file_id, "/ctrlPts_loc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  setid_x = H5Dcreate( group_id, "ctrlPts_x_loc", H5T_NATIVE_DOUBLE, dataspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ); 
  setid_y = H5Dcreate( group_id, "ctrlPts_y_loc", H5T_NATIVE_DOUBLE, dataspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ); 
  setid_z = H5Dcreate( group_id, "ctrlPts_z_loc", H5T_NATIVE_DOUBLE, dataspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ); 
  setid_w = H5Dcreate( group_id, "ctrlPts_w_loc", H5T_NATIVE_DOUBLE, dataspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  H5Dwrite( setid_x, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &ctrlPts_x_loc[0] ); 
  H5Dwrite( setid_y, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &ctrlPts_y_loc[0] ); 
  H5Dwrite( setid_z, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &ctrlPts_z_loc[0] ); 
  H5Dwrite( setid_w, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      &ctrlPts_w_loc[0] );

  H5Dclose( setid_x );
  H5Dclose( setid_y );
  H5Dclose( setid_z );
  H5Dclose( setid_w );

  H5Sclose( dataspace_id );
  H5Gclose( group_id );
}

void Part_NURBS_1Patch_3D_METIS::write_extractor( hid_t file_id ) const
{
  hid_t group_id, dataspace_id_x, dataspace_id_y, dataspace_id_z;
  hid_t setid_x, setid_y, setid_z;
  hsize_t dim_x[2], dim_y[2], dim_z[2];

  group_id = H5Gcreate(file_id, "/Extraction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  dim_x[0] = nlocalele; dim_x[1] = (sDegree+1)*(sDegree+1); 
  dim_y[0] = nlocalele; dim_y[1] = (tDegree+1)*(tDegree+1); 
  dim_z[0] = nlocalele; dim_z[1] = (uDegree+1)*(uDegree+1); 

  dataspace_id_x = H5Screate_simple(2, dim_x, NULL);
  dataspace_id_y = H5Screate_simple(2, dim_y, NULL);
  dataspace_id_z = H5Screate_simple(2, dim_z, NULL);

  setid_x = H5Dcreate( group_id, "extractor_x", H5T_NATIVE_DOUBLE, dataspace_id_x,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_y = H5Dcreate( group_id, "extractor_y", H5T_NATIVE_DOUBLE, dataspace_id_y,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  setid_z = H5Dcreate( group_id, "extractor_z", H5T_NATIVE_DOUBLE, dataspace_id_z,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

  double * temp_ext_x = new double [nlocalele * (sDegree+1) * (sDegree+1)]; 
  double * temp_ext_y = new double [nlocalele * (tDegree+1) * (tDegree+1)]; 
  double * temp_ext_z = new double [nlocalele * (uDegree+1) * (uDegree+1)]; 

  for(int e=0; e<nlocalele; ++e)
  {
    int e_global = elem_loc[e];
    int e_global_x, e_global_y, e_global_z;
    SYS_T::get_xyz_index(e_global, nElem_x, nElem_y, e_global_x,
        e_global_y, e_global_z);
    std::vector<double> e_ext_x, e_ext_y, e_ext_z;

    Get_bezier_ext_x(e_global_x, e_ext_x);
    Get_bezier_ext_y(e_global_y, e_ext_y);
    Get_bezier_ext_z(e_global_z, e_ext_z);

    for(int ii=0; ii<(sDegree+1)*(sDegree+1); ++ii)
      temp_ext_x[e*(sDegree+1)*(sDegree+1)+ii] = e_ext_x[ii];

    for(int ii=0; ii<(tDegree+1)*(tDegree+1); ++ii)
      temp_ext_y[e*(tDegree+1)*(tDegree+1)+ii] = e_ext_y[ii];

    for(int ii=0; ii<(uDegree+1)*(uDegree+1); ++ii)
      temp_ext_z[e*(uDegree+1)*(uDegree+1)+ii] = e_ext_z[ii];
  }

  H5Dwrite( setid_x, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      temp_ext_x );

  H5Dwrite( setid_y, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      temp_ext_y );

  H5Dwrite( setid_z, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      temp_ext_z );

  H5Dclose(setid_x);
  H5Dclose(setid_y);
  H5Dclose(setid_z);

  delete [] temp_ext_x; delete [] temp_ext_y; delete [] temp_ext_z;

  HDF5_Writer * hdf5writer = new HDF5_Writer(file_id);
  hdf5writer->write_intScalar(group_id, "elemType", elemType);
   
  delete hdf5writer;
  H5Sclose( dataspace_id_x );
  H5Sclose( dataspace_id_y );
  H5Sclose( dataspace_id_z );
  H5Gclose(group_id);
}

void Part_NURBS_1Patch_3D_METIS::write_hxyz( hid_t file_id ) const
{
  HDF5_Writer * h5writer = new HDF5_Writer( file_id );
   
  hid_t group_id;

  group_id = H5Gcreate(file_id, "/Mesh_Size", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  h5writer->write_doubleVector( group_id, "hx", hx );
  h5writer->write_doubleVector( group_id, "hy", hy );
  h5writer->write_doubleVector( group_id, "hz", hz );

  H5Gclose(group_id);
  delete h5writer;
}

bool Part_NURBS_1Patch_3D_METIS::isElemInPart(int gloindex) const
{
  std::vector<int>::const_iterator findindex;
  findindex = find(elem_loc.begin(), elem_loc.end(), gloindex);
  return findindex != elem_loc.end();
}


int Part_NURBS_1Patch_3D_METIS::get_elemLocIndex(const int &gloindex) const
{
  std::vector<int>::const_iterator findindex;
  findindex = find(elem_loc.begin(), elem_loc.end(), gloindex);
  if(findindex == elem_loc.end())
    return -1;
  else
    return (findindex - elem_loc.begin());
}


bool Part_NURBS_1Patch_3D_METIS::isNodeInPart(int gloindex) const
{
  std::vector<int>::const_iterator findindex;
  findindex = find(node_loc.begin(), node_loc.end(), gloindex);
  return findindex != node_loc.end();
}


// EOF
