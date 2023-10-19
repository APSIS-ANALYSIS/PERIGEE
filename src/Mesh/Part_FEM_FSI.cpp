#include "Part_FEM_FSI.hpp"

Part_FEM_FSI::Part_FEM_FSI( const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const std::vector<int> &phytag,
    const std::vector<int> &node_f,
    const std::vector<int> &node_s,
    const int &in_cpu_rank, 
    const int &in_cpu_size,
    const int &in_elemType,
    const int &field,
    const int &in_dof,
    const int &in_start_idx,
    const bool &in_is_geo_field ) 
: Part_FEM(), start_idx( in_start_idx ), is_geo_field(in_is_geo_field)
{
  nElem = mesh->get_nElem(); 
  nFunc = mesh->get_nFunc();
  sDegree = mesh->get_s_degree(); 
  tDegree = mesh->get_t_degree();
  uDegree = mesh->get_u_degree(); 
  nLocBas = mesh->get_nLocBas();
  probDim = 3; 
  elemType = in_elemType;
  cpu_rank = in_cpu_rank;
  cpu_size = in_cpu_size;

  // We set dofMat and dofNum to negative numbers, as they are not needed
  dofMat = -1;
  dofNum = in_dof;

  // Check the cpu info
  SYS_T::print_fatal_if(cpu_size < 1, "Error: Part_FEM_FSI input cpu_size is wrong! \n");
  SYS_T::print_fatal_if(cpu_rank >= cpu_size, "Error: Part_FEM_FSI input cpu_rank is wrong! \n");
  SYS_T::print_fatal_if(cpu_rank < 0, "Error: Part_FEM_FSI input cpu_rank is wrong! \n");

  // Generate group 1, 2, and 5.
  Generate_Partition( mesh, gpart, mnindex, IEN, field );

  // Generate group 6 if this field is geometry
  if( is_geo_field )
  {
    // local copy of control points
    ctrlPts_x_loc.resize(nlocghonode);
    ctrlPts_y_loc.resize(nlocghonode);
    ctrlPts_z_loc.resize(nlocghonode);
    PERIGEE_OMP_PARALLEL_FOR
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

    std::cout<<"-- proc "<<cpu_rank<<" Local control points generated. \n";
  }
  else
  {
    ctrlPts_x_loc.clear();
    ctrlPts_y_loc.clear();
    ctrlPts_z_loc.clear();
  }

  // Generate the local array tagging the element's property.
  elem_phy_tag.resize( nlocalele );
  for(int ii=0; ii<nlocalele; ++ii) elem_phy_tag[ii] = phytag[ elem_loc[ii] ];

  // Generate the node_loc_fluid/solid
  node_loc_fluid.clear();
  node_loc_solid.clear();
  PERIGEE_OMP_PARALLEL
  {
    std::vector<int> temp_node_loc_fluid {};
    std::vector<int> temp_node_loc_solid {};
    PERIGEE_OMP_FOR
    for(int ii=0; ii<nlocalnode; ++ii)
    {
      if( VEC_T::is_invec(node_f, node_loc_original[ii]) ) temp_node_loc_fluid.push_back(ii);

      if( VEC_T::is_invec(node_s, node_loc_original[ii]) ) temp_node_loc_solid.push_back(ii);
    }
    PERIGEE_OMP_CRITICAL
    {
      VEC_T::insert_end(node_loc_fluid, temp_node_loc_fluid);
      VEC_T::insert_end(node_loc_solid, temp_node_loc_solid);
    }
  }

  nlocalnode_fluid = static_cast<int>( node_loc_fluid.size() );
  nlocalnode_solid = static_cast<int>( node_loc_solid.size() );
}

Part_FEM_FSI::~Part_FEM_FSI()
{}

void Part_FEM_FSI::write( const std::string &inputFileName ) const
{
  const std::string fName = SYS_T::gen_partfile_name( inputFileName, cpu_rank );

  hid_t file_id = H5Fcreate(fName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  // group 1: local element
  hid_t group_id_1 = H5Gcreate(file_id, "/Local_Elem", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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
  hid_t group_id_3 = H5Gcreate(file_id, "/Global_Mesh_Info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  h5w->write_intScalar( group_id_3, "nElem", nElem );
  h5w->write_intScalar( group_id_3, "nFunc", nFunc );

  const std::vector<int> deg { sDegree, tDegree, uDegree };

  h5w->write_intVector( group_id_3, "degree", deg );

  h5w->write_intScalar( group_id_3, "nLocBas", nLocBas );

  h5w->write_intScalar( group_id_3, "probDim", probDim );
  h5w->write_intScalar( group_id_3, "dofNum",  dofNum );
  h5w->write_intScalar( group_id_3, "elemType", elemType );

  H5Gclose( group_id_3 );

  // group 4: part info
  hid_t group_id_4 = H5Gcreate( file_id, "/Part_Info", H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT ); 

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
  if( is_geo_field)
  {
    hid_t group_id_6 = H5Gcreate(file_id, "/ctrlPts_loc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    h5w -> write_doubleVector( group_id_6, "ctrlPts_x_loc", ctrlPts_x_loc );
    h5w -> write_doubleVector( group_id_6, "ctrlPts_y_loc", ctrlPts_y_loc );
    h5w -> write_doubleVector( group_id_6, "ctrlPts_z_loc", ctrlPts_z_loc );

    H5Gclose( group_id_6 );
  }

  // group 7: DOF mapper
  hid_t group_id_7 = H5Gcreate(file_id, "/DOF_mapper", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  h5w -> write_intScalar( group_id_7, "start_idx", start_idx );

  H5Gclose( group_id_7 );

  // Finish writing, clean up
  delete h5w;
  H5Fclose(file_id);
}

// EOF
