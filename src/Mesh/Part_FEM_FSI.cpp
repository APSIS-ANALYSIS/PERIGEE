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
    const int &in_start_idx,
    const Field_Property * const &fp ) 
: Part_FEM( mesh, gpart, mnindex, IEN, ctrlPts, in_cpu_rank, in_cpu_size, in_dof, in_elemType, fp ), 
  start_idx( in_start_idx ), is_geo_field(in_is_geo_field)
{
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

  nlocalnode_fluid = VEC_T::get_size( node_loc_fluid );
  nlocalnode_solid = VEC_T::get_size( node_loc_solid );
}

Part_FEM_FSI::~Part_FEM_FSI()
{}

void Part_FEM_FSI::write( const std::string &inputFileName ) const
{
  // ------------------------------------------------------
  // Call the base class writer to write the base class data
  Part_FEM::write( inputFileName );
  // ------------------------------------------------------

  const std::string fName = SYS_T::gen_partfile_name( inputFileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  // open group 1: local element
  hid_t group_id_1 = H5Gopen(file_id, "/Local_Elem", H5P_DEFAULT);

  h5w->write_intVector( group_id_1, "elem_phy_tag", elem_phy_tag );

  H5Gclose( group_id_1 );

  // group 2: local node
  hid_t group_id_2 = H5Gopen( file_id, "/Local_Node", H5P_DEFAULT );

  h5w->write_intScalar( group_id_2, "nlocalnode_fluid", nlocalnode_fluid );
  h5w->write_intScalar( group_id_2, "nlocalnode_solid", nlocalnode_solid );

  if( nlocalnode_fluid > 0 )
    h5w->write_intVector( group_id_2, "node_loc_fluid", node_loc_fluid );

  if( nlocalnode_solid > 0 )
    h5w->write_intVector( group_id_2, "node_loc_solid", node_loc_solid );

  H5Gclose( group_id_2 );

  // group 7: DOF mapper
  hid_t group_id_7 = H5Gcreate(file_id, "/DOF_mapper", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  h5w -> write_intScalar( group_id_7, "start_idx", start_idx );

  H5Gclose( group_id_7 );

  // Finish the writing of hdf5 file
  delete h5w;
  H5Fclose(file_id);
}

// EOF
