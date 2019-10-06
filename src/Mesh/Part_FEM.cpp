#include "Part_FEM.hpp"

Part_FEM::Part_FEM( const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const int &in_cpu_rank, const int &in_cpu_size,
    const int &in_dofNum, const int &in_dofMat,
    const int &in_probdim, const int &in_etype, 
    const bool &isPrintInfo )
: cpu_rank( in_cpu_rank ), cpu_size( in_cpu_size ),
  isMETIS( gpart->get_isMETIS() ),
  part_isdual( gpart->get_isDual() ),
  dual_edge_ncommon( gpart->get_dual_edge_ncommon() ),
  nElem( mesh->get_nElem() ), nFunc( mesh->get_nFunc() ),
  sDegree( mesh->get_s_degree() ), tDegree( mesh->get_t_degree() ),
  uDegree( mesh->get_u_degree() ), nLocBas( mesh->get_nLocBas() ),
  probDim( in_probdim ), dofNum( in_dofNum ), dofMat(in_dofMat),
  elemType( in_etype )
{
  SYS_T::print_exit_if(cpu_rank >= cpu_size || cpu_size < 0,
      "Error: Part_FEM input cpu rank is wrong! \n");

  IPart::GenPart( nElem, nFunc, nLocBas, cpu_rank,
      gpart, mnindex, IEN, ctrlPts, isPrintInfo,
      elem_loc, nlocalele, node_loc, node_loc_original,
      node_ghost, local_to_global, nlocalnode, nghostnode,
      ntotalnode, nbadnode, nlocghonode, ctrlPts_x_loc,
      ctrlPts_y_loc, ctrlPts_z_loc, LIEN );
}


Part_FEM::~Part_FEM()
{}


void Part_FEM::write( const char * inputFileName ) const
{
  const std::string bName( inputFileName );
  const std::string fName = SYS_T::gen_partfile_name( bName, cpu_rank );
  
  hid_t file_id = H5Fcreate(fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  HDF5_Writer * h5w = new HDF5_Writer( file_id );

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
