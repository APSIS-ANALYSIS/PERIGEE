#include "MeshPart_NURBS_multiPatch_3D.hpp"


MeshPart_NURBS_multiPatch_3D::MeshPart_NURBS_multiPatch_3D(
    const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const int &in_cpu_size,
    const int &in_dofNum, const int &in_elemType )
: numPat(mesh->get_num_patch()),
  cpu_size(in_cpu_size), dual_edge_ncommon(gpart->get_dual_edge_ncommon()),
  isMETIS(gpart->get_isMETIS()), part_isdual(gpart->get_isDual()),
  nElem(mesh->get_nElem()), nFunc(mesh->get_nFunc()),
  sDegree(mesh->get_s_degree()), tDegree(mesh->get_t_degree()),
  uDegree(mesh->get_u_degree()), nLocBas(mesh->get_nLocBas()),
  hx_max(mesh->get_hx_max()), hy_max(mesh->get_hy_max()),
  hz_max(mesh->get_hz_max()), hx_min(mesh->get_hx_min()),
  hy_min(mesh->get_hy_min()), hz_min(mesh->get_hz_min()),
  probDim(3), dofNum(in_dofNum), elemType(in_elemType)
{
  if(cpu_size<1)
  {
    std::cerr<<"Error: given cpu_size = "<<cpu_size<<" is negative. \n";
    exit(EXIT_FAILURE); 
  }

  // Calculate the nonzero index for each patch
  nonzero_x.resize(numPat);
  nonzero_y.resize(numPat);
  nonzero_z.resize(numPat);

  int temp_index = 0;
  const IMesh * patch_mesh = NULL;

  for(int ii=0; ii<numPat; ++ii)
  {
    nonzero_x[ii].clear();
    nonzero_y[ii].clear();
    nonzero_z[ii].clear();

    temp_index = 0;
    
    patch_mesh = mesh->get_patch_mesh(ii);

    const int nelemx = patch_mesh->get_nElem_x();
    const int nelemy = patch_mesh->get_nElem_y();
    const int nelemz = patch_mesh->get_nElem_z();
    
    for(int ee=0; ee<nelemx; ++ee)
    {
      if(patch_mesh->get_hx(ee) == 0.0)
        nonzero_x[ii].push_back(-1);
      else
      {
        nonzero_x[ii].push_back(temp_index);
        temp_index += 1;
      }
    }

    temp_index = 0;

    for(int ee=0; ee<nelemy; ++ee)
    {
      if(patch_mesh->get_hy(ee*nelemx) == 0.0)
        nonzero_y[ii].push_back(-1);
      else
      {
        nonzero_y[ii].push_back(temp_index);
        temp_index += 1;
      }
    }

    temp_index = 0;

    for(int ee=0; ee<nelemz; ++ee)
    {
      if(patch_mesh->get_hz(ee*nelemx*nelemy) == 0.0)
        nonzero_z[ii].push_back(-1);
      else
      {
        nonzero_z[ii].push_back(temp_index);
        temp_index += 1;
      }
    }

    VEC_T::shrink2fit(nonzero_x[ii]);
    VEC_T::shrink2fit(nonzero_y[ii]);
    VEC_T::shrink2fit(nonzero_z[ii]);
  }
  
  // rank dependent parameters are set to default values.
  // when calling mesh partition routine, they will be reset.
  cpu_rank = 0;
  
  elem_loc.clear(); nlocalele = 0;

  node_loc.clear(); node_loc_original.clear(); node_ghost.clear();
  local_to_global.clear();
  nlocalnode = 0; nghostnode = 0; ntotalnode = 0; nbadnode = 0;
  nlocghonode = 0;
  
  hx.clear(); hy.clear(); hz.clear();
  LIEN.clear();
  ctrlPts_x_loc.clear(); ctrlPts_y_loc.clear(); 
  ctrlPts_z_loc.clear(); ctrlPts_w_loc.clear();
}



MeshPart_NURBS_multiPatch_3D::~MeshPart_NURBS_multiPatch_3D()
{
  for(int ii=0; ii<numPat; ++ii)
  {
    VEC_T::clean(nonzero_x[ii]);
    VEC_T::clean(nonzero_y[ii]);
    VEC_T::clean(nonzero_z[ii]);
  }
  VEC_T::clean(nonzero_x);
  VEC_T::clean(nonzero_y);
  VEC_T::clean(nonzero_z);

  std::cout<<"-- MeshPart_NURBS_multiPatch_3D deleted. \n";
}


void MeshPart_NURBS_multiPatch_3D::set_part_info( const int &in_cpu_rank,
    const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const bool &isPrintInfo )
{
  cpu_rank = in_cpu_rank;

  if(cpu_rank >= cpu_size)
  {
    std::cerr<<"Error: cpu_rank = "<<cpu_rank<<" > cpu_size = "<<cpu_size<<std::endl;
    exit(EXIT_FAILURE);
  }

  elem_loc.clear(); node_loc.clear(); node_loc_original.clear();

  for( int ee =0; ee<nElem; ++ee )
  {
    if(gpart->get_epart(ee) == cpu_rank)
      elem_loc.push_back(ee);
  }
  VEC_T::shrink2fit(elem_loc);
  nlocalele = elem_loc.size();


  for( int ii=0; ii<nFunc; ++ii )
  {
    if(gpart->get_npart(ii) == cpu_rank)
    {
      node_loc.push_back(ii);
      node_loc_original.push_back(ii);
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

  for(int ii=0; ii<nlocalnode; ++ii)
    node_loc[ii] = mnindex->get_old2new(node_loc[ii]);

  std::vector<int> node_tot;
  node_tot.clear();
  int temp_node;
  for(int ee=0; ee<nlocalele; ++ee)
  {
    for(int ii=0; ii<nLocBas; ++ii)
    {
      temp_node = IEN->get_IEN(elem_loc[ee], ii);
      temp_node = mnindex->get_old2new(temp_node);
      node_tot.push_back(temp_node);
    }
  }
  VEC_T::sort_unique_resize(node_tot);
  ntotalnode = (int) node_tot.size();

  node_ghost.clear();
  std::vector<int>::const_iterator it;
  for(int ii=0; ii<ntotalnode; ++ii)
  {
    it = find(node_loc.begin(), node_loc.end(), node_tot[ii]);
    if(it == node_loc.end())
      node_ghost.push_back(node_tot[ii]);
  }
  VEC_T::shrink2fit(node_ghost);
  nghostnode = (int) node_ghost.size();

  nbadnode = 0;
  if( nghostnode + nlocalnode != ntotalnode )
  {
    std::vector<int> badnode;
    std::vector<int>::const_iterator badnode_it;
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


  local_to_global.clear();
  for( int n=0; n<nlocalnode; ++n )
    local_to_global.push_back( node_loc[n] );
  for( int n=0; n<nghostnode; ++n )
    local_to_global.push_back( node_ghost[n] );

  VEC_T::shrink2fit(local_to_global);
  nlocghonode = (int) local_to_global.size();

  if( isPrintInfo )
    std::cout<<"-- proc "<<cpu_rank<<" local_to_global generated. \n";


  LIEN.resize(nlocalele * nLocBas);

  std::vector<int>::const_iterator lien_ptr;
  int global_index;
  int counter = 0;
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

      LIEN[counter] = lien_ptr - local_to_global.begin();
      counter += 1;
    }
  }

  if(isPrintInfo)
    std::cout<<"-- proc "<<cpu_rank<<" LIEN generated. \n";


  ctrlPts_x_loc.resize(nlocghonode);
  ctrlPts_y_loc.resize(nlocghonode);
  ctrlPts_z_loc.resize(nlocghonode);
  ctrlPts_w_loc.resize(nlocghonode);

  int aux_index;
  for(int ii=0; ii<nlocghonode; ++ii)
  {
    aux_index = local_to_global[ii]; // new global index
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

  hx.clear(); hy.clear(); hz.clear();
  int e_global;
  for( int e=0; e<nlocalele; ++e )
  {
    e_global = get_elem_loc(e);
    hx.push_back(mesh->get_hx(e_global));
    hy.push_back(mesh->get_hy(e_global));
    hz.push_back(mesh->get_hz(e_global));
  }

}



void MeshPart_NURBS_multiPatch_3D::write( const std::string &inputFileName,
    const IMesh * const &mesh,
    const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_x,
    const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_y,
    const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_z ) const
{
  std::string fName = inputFileName;
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

  hid_t file_id, group_id;
  file_id = H5Fcreate(fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);  

  std::vector<int> temp_int;
  std::vector<double> temp_double;

  // 1. Local element
  group_id = H5Gcreate(file_id, "/Local_Elem", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5w->write_intScalar(group_id, "nlocalele", nlocalele);
  h5w->write_intVector(group_id, "elem_loc", elem_loc);
  H5Gclose(group_id);

  // 2. Local node
  group_id = H5Gcreate(file_id, "/Local_Node", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5w->write_intScalar(group_id, "nlocalnode", nlocalnode);
  h5w->write_intScalar(group_id, "nghostnode", nghostnode);
  h5w->write_intScalar(group_id, "ntotalnode", ntotalnode);
  h5w->write_intScalar(group_id, "nbadnode", nbadnode);
  h5w->write_intScalar(group_id, "nlocghonode", nlocghonode);
  h5w->write_intVector(group_id, "node_loc", node_loc);
  h5w->write_intVector(group_id, "node_loc_original", node_loc_original);
  if(nghostnode > 0)
    h5w->write_intVector(group_id, "node_ghost", node_ghost);
  h5w->write_intVector(group_id, "local_to_global", local_to_global);
  H5Gclose(group_id);

  // 3. Part info
  group_id = H5Gcreate(file_id, "/Part_Info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5w->write_intScalar(group_id, "cpu_rank", cpu_rank);
  h5w->write_intScalar(group_id, "cpu_size", cpu_size);
  h5w->write_intScalar(group_id, "dual_edge_ncommon", dual_edge_ncommon);
  H5Gclose(group_id);

  // 4. Global mesh info
  group_id = H5Gcreate(file_id, "/Global_Mesh_Info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5w->write_intScalar(group_id, "nElem", nElem);
  h5w->write_intScalar(group_id, "nFunc", nFunc);
  temp_int.clear(); 
  temp_int.push_back(sDegree);
  temp_int.push_back(tDegree);
  temp_int.push_back(uDegree);
  h5w->write_intVector(group_id, "degree", temp_int);
  temp_double.clear();
  temp_double.push_back(hx_max);
  temp_double.push_back(hy_max);
  temp_double.push_back(hz_max);
  h5w->write_doubleVector(group_id, "h_max",temp_double);
  temp_double.clear();
  temp_double.push_back(hx_min);
  temp_double.push_back(hy_min);
  temp_double.push_back(hz_min);
  h5w->write_doubleVector(group_id, "h_min",temp_double);
  h5w->write_intScalar(group_id, "nLocBas", nLocBas);
  h5w->write_intScalar(group_id, "probDim", probDim);
  h5w->write_intScalar(group_id, "dofNum", dofNum);
  H5Gclose(group_id);

  // 5. LIEN
  group_id = H5Gcreate(file_id, "/LIEN", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5w->write_intMatrix(group_id, "LIEN", LIEN, nlocalele, nLocBas ); 
  H5Gclose(group_id);

  // 6. Control points
  group_id = H5Gcreate(file_id, "/ctrlPts_loc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5w->write_doubleVector(group_id, "ctrlPts_x_loc", ctrlPts_x_loc);
  h5w->write_doubleVector(group_id, "ctrlPts_y_loc", ctrlPts_y_loc);
  h5w->write_doubleVector(group_id, "ctrlPts_z_loc", ctrlPts_z_loc);
  h5w->write_doubleVector(group_id, "ctrlPts_w_loc", ctrlPts_w_loc);
  H5Gclose(group_id);


  // 7. extraction operators
  std::vector<double> vec_ext_x, vec_ext_y, vec_ext_z;
  vec_ext_x.reserve(nlocalele * (sDegree+1) * (sDegree+1));
  vec_ext_y.reserve(nlocalele * (tDegree+1) * (tDegree+1));
  vec_ext_z.reserve(nlocalele * (uDegree+1) * (uDegree+1));

  vec_ext_x.clear(); vec_ext_y.clear(); vec_ext_z.clear();

  std::vector<double> e_ext_x, e_ext_y, e_ext_z;
  for(int ee=0; ee<nlocalele; ++ee)
  {
    int e_global = elem_loc[ee];

    Get_bezier_ext(e_global, bez_x, bez_y, bez_z, mesh, e_ext_x, e_ext_y, e_ext_z);

    vec_ext_x.insert(vec_ext_x.end(), e_ext_x.begin(), e_ext_x.end());
    vec_ext_y.insert(vec_ext_y.end(), e_ext_y.begin(), e_ext_y.end());
    vec_ext_z.insert(vec_ext_z.end(), e_ext_z.begin(), e_ext_z.end());
  }

  group_id = H5Gcreate(file_id, "/Extraction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5w->write_doubleMatrix(group_id, "extractor_x", vec_ext_x, nlocalele, (sDegree+1)*(sDegree+1));
  h5w->write_doubleMatrix(group_id, "extractor_y", vec_ext_y, nlocalele, (tDegree+1)*(tDegree+1));
  h5w->write_doubleMatrix(group_id, "extractor_z", vec_ext_z, nlocalele, (uDegree+1)*(uDegree+1));
  h5w->write_intScalar(group_id, "elemType", elemType);
  H5Gclose(group_id);

  // 8. hx hy hz
  group_id = H5Gcreate(file_id, "/Mesh_Size", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5w->write_doubleVector( group_id, "hx", hx );
  h5w->write_doubleVector( group_id, "hy", hy );
  h5w->write_doubleVector( group_id, "hz", hz );
  H5Gclose(group_id);

  delete h5w;
  H5Fclose(file_id);
}


void MeshPart_NURBS_multiPatch_3D::Get_bezier_ext( const int &elem,
    const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_x,
    const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_y,
    const std::vector< std::vector<NURBS_T::BezierElem*> > &bez_z,
    const IMesh * const &mesh, std::vector<double> &ext_x,
    std::vector<double> &ext_y, std::vector<double> &ext_z ) const
{
  ext_x.clear(); ext_x.resize((sDegree+1)*(sDegree+1));
  ext_y.clear(); ext_y.resize((tDegree+1)*(tDegree+1));
  ext_z.clear(); ext_z.resize((uDegree+1)*(uDegree+1));


  int p_index, loc_ee;
  mesh->get_locelem_index(elem, p_index, loc_ee);

  const int p_nelemx = mesh->get_patch_mesh(p_index)->get_nElem_x();
  const int p_nelemy = mesh->get_patch_mesh(p_index)->get_nElem_y();

  int loc_ee_x, loc_ee_y, loc_ee_z;
  SYS_T::get_xyz_index( loc_ee, p_nelemx, p_nelemy, loc_ee_x, loc_ee_y, loc_ee_z );

  const int nz_index_x = nonzero_x[p_index][loc_ee_x];
  const int nz_index_y = nonzero_y[p_index][loc_ee_y];
  const int nz_index_z = nonzero_z[p_index][loc_ee_z];

  if(nz_index_x == -1)
  {
    for(int ii=0; ii<(sDegree+1)*(sDegree+1); ++ii)
      ext_x[ii] = 0.0;
  }
  else
  {
    for(int ii=0; ii<sDegree+1; ++ii)
    {
      for(int jj=0; jj<sDegree+1; ++jj)
        ext_x[ii*(sDegree+1) + jj] = bez_x[p_index][nz_index_x]->coefficient(ii,jj);
    }
  }

  if(nz_index_y == -1)
  {
    for(int ii=0; ii<(tDegree+1)*(tDegree+1); ++ii)
      ext_y[ii] = 0.0;
  }
  else
  {
    for(int ii=0; ii<tDegree+1; ++ii)
    {
      for(int jj=0; jj<tDegree+1; ++jj)
        ext_y[ii*(tDegree+1)+jj] = bez_y[p_index][nz_index_y]->coefficient(ii,jj);
    }
  }


  if(nz_index_z == -1)
  {
    for(int ii=0; ii<(uDegree+1)*(uDegree+1); ++ii)
      ext_z[ii] = 0.0;
  }
  else
  {
    for(int ii=0; ii<uDegree+1; ++ii)
    {
      for(int jj=0; jj<uDegree+1; ++jj)
        ext_z[ii*(uDegree+1)+jj] = bez_z[p_index][nz_index_z]->coefficient(ii,jj);
    }
  }

}


inline bool MeshPart_NURBS_multiPatch_3D::isElemInPart(const int &eindex) const
{
  std::vector<int>::const_iterator finder = find(elem_loc.begin(), elem_loc.end(), eindex);
  return finder != elem_loc.end();
}


inline bool MeshPart_NURBS_multiPatch_3D::isNodeInPart(const int &nindex) const
{
  std::vector<int>::const_iterator finder = find(node_loc.begin(), node_loc.end(), nindex);
  return finder != node_loc.end();
}



inline int MeshPart_NURBS_multiPatch_3D::get_elemLocIndex(const int &pos) const
{
  std::vector<int>::const_iterator finder = find(elem_loc.begin(), elem_loc.end(), pos);
  if(finder == elem_loc.end())
    return -1;
  else
    return (finder - elem_loc.begin());
}



// EOF
