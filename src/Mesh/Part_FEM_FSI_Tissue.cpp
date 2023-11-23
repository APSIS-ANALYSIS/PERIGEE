#include "Part_FEM_FSI_Tissue.hpp"

Part_FEM_FSI_Tissue::Part_FEM_FSI_Tissue( const IMesh * const &mesh,
    const IGlobal_Part * const &gpart,
    const Map_Node_Index * const &mnindex,
    const IIEN * const &IEN,
    const std::vector<double> &ctrlPts,
    const std::vector<int> &phytag,
    const std::vector<int> &node_f,
    const std::vector<int> &node_s,
    const std::vector<Vector_3> &basis_r,
    const std::vector<Vector_3> &basis_l,
    const std::vector<Vector_3> &basis_c,
    const int &in_cpu_rank, 
    const int &in_cpu_size,
    const int &in_elemType,
    const int &in_start_idx,
    const Field_Property &fp ) 
: Part_FEM_FSI( mesh, gpart, mnindex, IEN, ctrlPts, phytag, node_f, node_s, in_cpu_rank, in_cpu_size, in_elemType, in_start_idx, fp )
{
  // Generate the node_locgho_solid
  node_locgho_solid.resize(nlocghonode);
  loc_basis_r.clear();
  loc_basis_l.clear();
  loc_basis_c.clear();

  int index = 0;
  for(int ii=0; ii<nlocghonode; ++ii)
  {
    int aux_index = local_to_global[ii]; // new global index
    aux_index = mnindex->get_new2old(aux_index); // back to old global index
    if( VEC_T::is_invec(node_s, aux_index) )
    {
      node_locgho_solid[ii] = index; // record local solid node index
      index += 1;

      const int pos = VEC_T::get_pos( node_s, aux_index );
      loc_basis_r.push_back( basis_r[pos] );
      loc_basis_l.push_back( basis_l[pos] );
      loc_basis_c.push_back( basis_c[pos] );
    }
    else node_locgho_solid[ii] = -1;
  }

  nlocghonode_s = VEC_T::get_size( node_locgho_solid );
  std::cout<<"-- proc "<<cpu_rank<<" Local direction basis vectors generated. \n";
}

void Part_FEM_FSI_Tissue::write( const std::string &inputFileName ) const
{
  // ------------------------------------------------------
  // Call the base class writer to write the base class data
  Part_FEM_FSI::write( inputFileName );
  // ------------------------------------------------------

  const std::string fName = SYS_T::gen_partfile_name( inputFileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer(file_id);

  hid_t group_id_8 = H5Gcreate(file_id, "/directionBasis_loc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  h5w -> write_intScalar( group_id_8, "nlocghonode_s", nlocghonode_s );
  h5w -> write_intVector( group_id_8, "node_locgho_solid", node_locgho_solid );
  h5w -> write_Vector3_Vector( group_id_8, "loc_basis_r", loc_basis_r );
  h5w -> write_Vector3_Vector( group_id_8, "loc_basis_l", loc_basis_l );
  h5w -> write_Vector3_Vector( group_id_8, "loc_basis_c", loc_basis_c );

  H5Gclose( group_id_8 );

  // Finish the writing of hdf5 file
  delete h5w;
  H5Fclose(file_id);
}

// EOF
