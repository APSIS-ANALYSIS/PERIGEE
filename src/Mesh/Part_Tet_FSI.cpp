#include "Part_Tet_FSI.hpp"

Part_Tet_FSI::Part_Tet_FSI( const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const std::vector<int> &phytag,
    const std::vector<int> &node_f,
    const std::vector<int> &node_s,
    const int &in_cpu_rank, const int &in_cpu_size,
    const int &in_dofNum, const int &in_dofMat,
    const int &in_elemType,
    const bool isPrintInfo )
: Part_Tet( mesh, gpart, mnindex, IEN, ctrlPts, in_cpu_rank,
    in_cpu_size, in_dofNum, in_dofMat, in_elemType, isPrintInfo )
{
  // Generate the local array tagging the element's property.
  elem_phy_tag.resize( nlocalele );
  for(int ii=0; ii<nlocalele; ++ii) elem_phy_tag[ii] = phytag[ elem_loc[ii] ];

  // Generate the node_loc_fluid/solid
  node_loc_fluid.clear();
  node_loc_solid.clear();
  
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    if( VEC_T::is_invec(node_f, node_loc_original[ii]) ) node_loc_fluid.push_back(ii);
    
    if( VEC_T::is_invec(node_s, node_loc_original[ii]) ) node_loc_solid.push_back(ii);
  } 

  nlocalnode_fluid = static_cast<int>( node_loc_fluid.size() );
  nlocalnode_solid = static_cast<int>( node_loc_solid.size() );
}


Part_Tet_FSI::~Part_Tet_FSI()
{}


void Part_Tet_FSI::write( const char * inputFileName ) const
{
  std::string inputName( inputFileName );
  std::string fName = SYS_T::gen_partfile_name( inputName, cpu_rank );

  hid_t file_id;
  file_id = H5Fcreate(fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  // group 1: local element
  hid_t group_id_1 = H5Gcreate(file_id, "/Local_Elem", H5P_DEFAULT, 
      H5P_DEFAULT, H5P_DEFAULT);
  
  h5w->write_intScalar( group_id_1, "nlocalele", nlocalele );
  h5w->write_intVector( group_id_1, "elem_loc", elem_loc );
  h5w->write_intVector( group_id_1, "elem_phy_tag", elem_phy_tag );

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
  
  h5w->write_intScalar( group_id_2, "nlocalnode_fluid", nlocalnode_fluid );
  h5w->write_intScalar( group_id_2, "nlocalnode_solid", nlocalnode_solid );

  if( nlocalnode_fluid > 0 )
    h5w->write_intVector( group_id_2, "node_loc_fluid", node_loc_fluid );

  if( nlocalnode_solid > 0 )
    h5w->write_intVector( group_id_2, "node_loc_solid", node_loc_solid );

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

// EOF
