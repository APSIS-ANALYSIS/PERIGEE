#include "Part_FEM_Rotated.hpp"

Part_FEM_Rotated::Part_FEM_Rotated( const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const std::vector<int> &eletag,
    const std::vector<int> &node_f,
    const std::vector<int> &node_r,
    const int &in_cpu_rank, 
    const int &in_cpu_size,
    const FEType &in_elemType,
    const Field_Property &fp ) 
: Part_FEM( mesh, gpart, mnindex, IEN, ctrlPts, in_cpu_rank, in_cpu_size, in_elemType, fp )
{
  // Generate the local array tagging the element's property.
  elem_tag.resize( nlocalele );
  for(int ii=0; ii<nlocalele; ++ii) elem_tag[ii] = eletag[ elem_loc[ii] ];

  // Generate the node_loc_fixed/rotated
  node_loc_fixed.clear();
  node_loc_rotated.clear();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    if( VEC_T::is_invec(node_f, node_loc_original[ii]) ) node_loc_fixed.push_back(ii);

    if( VEC_T::is_invec(node_r, node_loc_original[ii]) ) node_loc_rotated.push_back(ii);
  }

  nlocalnode_fixed = VEC_T::get_size( node_loc_fixed );
  nlocalnode_rotated = VEC_T::get_size( node_loc_rotated );
}

void Part_FEM_Rotated::write( const std::string &inputFileName ) const
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

  h5w->write_intVector( group_id_1, "elem_phy_tag", elem_tag );

  H5Gclose( group_id_1 );

  // group 2: local node
  hid_t group_id_2 = H5Gopen( file_id, "/Local_Node", H5P_DEFAULT );

  h5w->write_intScalar( group_id_2, "nlocalnode_fixed", nlocalnode_fixed );
  h5w->write_intScalar( group_id_2, "nlocalnode_rotated", nlocalnode_rotated );

  if( nlocalnode_fixed > 0 )
    h5w->write_intVector( group_id_2, "node_loc_fixed", node_loc_fixed );

  if( nlocalnode_rotated > 0 )
    h5w->write_intVector( group_id_2, "node_loc_rotated", node_loc_rotated );

  H5Gclose( group_id_2 );

  // Finish the writing of hdf5 file
  delete h5w;
  H5Fclose(file_id);
}

// EOF
