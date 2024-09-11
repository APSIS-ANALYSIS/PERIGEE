#include "ALocal_Interface.hpp"

ALocal_Interface::ALocal_Interface( const std::string &fileBaseName, const int &cpu_rank)
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/sliding");

  num_itf = h5r -> read_intScalar( gname.c_str(), "num_interface" );
  SYS_T::print_fatal_if(num_itf < 1, "Error, ALocal_Interface: there is no interface in this geometric model.\n");

  num_fixed_ele = h5r -> read_intVector( gname.c_str(), "num_part_fixed_cell" );
  max_num_fixed_ele = h5r -> read_intVector( gname.c_str(), "max_num_fixed_cell" );
  num_tag = h5r -> read_intVector( gname.c_str(), "num_tag" );

  std::string groupbase(gname);
  groupbase.append("/interfaceid_");

  num_fixed_node.assign(num_itf, 0);
  num_rotated_node.assign(num_itf, 0);

  num_rotated_ele.resize(num_itf);

  fixed_ele_face_id.resize(num_itf);
  fixed_ele_tag.resize(num_itf);
  fixed_lien.resize(num_itf);
  fixed_pt_xyz.resize(num_itf);
  fixed_node_id.resize(num_itf);
  fixed_LID.resize(num_itf);

  rotated_lien.resize(num_itf);
  rotated_ele_face_id.resize(num_itf);
  rotated_pt_xyz.resize(num_itf);
  rotated_node_id.resize(num_itf);
  rotated_LID.resize(num_itf);

  const std::string mesh_info("/Global_Mesh_Info");
  nLocBas = h5r -> read_intScalar(mesh_info.c_str(), "nLocBas");

  for(int ii=0; ii<num_itf; ++ii)
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( std::to_string(ii) );

    fixed_ele_face_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_cell_face_id" );

    fixed_ele_tag[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_cell_tag" );

    fixed_lien[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_cell_ien" );

    fixed_pt_xyz[ii] = h5r -> read_doubleVector( subgroup_name.c_str(), "fixed_pt_xyz" );

    fixed_node_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_node_map" );

    fixed_LID[ii] = h5r -> read_intVector(  subgroup_name.c_str(), "fixed_LID" );

    num_fixed_node[ii] = VEC_T::get_size(fixed_node_id[ii]);

    rotated_lien[ii].resize(num_tag[ii]);
    rotated_ele_face_id[ii].resize(num_tag[ii]);
    num_rotated_ele[ii].resize(num_tag[ii]);
    
    std::string subsubgroupbase(subgroup_name);
    subsubgroupbase.append("/tag_");

    for(int jj=0; jj<num_tag[ii]; ++jj)
    {
      std::string subsubgroup_name(subsubgroupbase);
      subsubgroup_name.append( std::to_string(jj) );

      num_rotated_ele[ii][jj] = h5r -> read_intScalar( subsubgroup_name.c_str(), "num_rotated_cell" );

      if(num_rotated_ele[ii][jj] > 0)
      {
        rotated_lien[ii][jj] = h5r -> read_intVector( subsubgroup_name.c_str(), "rotated_cell_ien" );

        rotated_ele_face_id[ii][jj] = h5r -> read_intVector( subsubgroup_name.c_str(), "rotated_cell_face_id" );
      }
    }

    rotated_pt_xyz[ii] = h5r -> read_doubleVector( subgroup_name.c_str(), "rotated_pt_xyz" );

    rotated_node_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_node_map" );

    rotated_LID[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_LID" );

    num_rotated_node[ii] = VEC_T::get_size(rotated_node_id[ii]);
  }

  delete h5r; H5Fclose( file_id );
}

void ALocal_Interface::print_info() const
{
  SYS_T::commPrint("Interfaces: %d\n", num_itf);
}

void ALocal_Interface::get_fixed_ele_ctrlPts(const int &ii, const int &ee,
  double * const volctrl_x,  double * const volctrl_y,  double * const volctrl_z) const
{
  for(int nn{0}; nn < nLocBas; ++nn)
  {
    int node = get_fixed_lien(ii, nLocBas * ee + nn);

    volctrl_x[nn] = get_fixed_pt_xyz(ii, 3 * node);
    volctrl_y[nn] = get_fixed_pt_xyz(ii, 3 * node + 1);
    volctrl_z[nn] = get_fixed_pt_xyz(ii, 3 * node + 2);
  }
}

void ALocal_Interface::get_rotated_ele_ctrlPts(const int &ii, const int &tag, const int &ee,
  double * const volctrl_x, double * const volctrl_y, double * const volctrl_z) const
{
  for(int nn{0}; nn < nLocBas; ++nn)
  {
    int node = get_rotated_lien(ii, tag, nLocBas * ee + nn);

    volctrl_x[nn] = rotated_pt_xyz[ii][3 * node];
    volctrl_y[nn] = rotated_pt_xyz[ii][3 * node + 1];
    volctrl_z[nn] = rotated_pt_xyz[ii][3 * node + 2];
  }
}

// EOF
