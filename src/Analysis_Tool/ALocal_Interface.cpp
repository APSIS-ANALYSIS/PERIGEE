#include "ALocal_Interface.hpp"

ALocal_Interface::ALocal_Interface( const std::string &fileBaseName, const int &cpu_rank )
: num_itf {0}
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/sliding");

  num_itf = h5r -> read_intScalar( gname.c_str(), "num_interface" );
  SYS_T::print_fatal_if(num_itf < 1, "Error, ALocal_Interface: there is no interface in this geometric model.\n");

  num_fixed_ele = h5r -> read_intVector( gname.c_str(), "num_part_fixed_cell" );
  num_tag = h5r -> read_intVector( gname.c_str(), "num_tag" );

  std::string groupbase(gname);
  groupbase.append("/interfaceid_");

  num_fixed_node.assign(num_itf, 0);
  num_rotated_node.assign(num_itf, 0);

  num_rotated_ele.resize(num_itf);

  fixed_ele_face_id.resize(num_itf);
  fixed_ele_tag.resize(num_itf);
  fixed_layer_ien.resize(num_itf);
  fixed_node_xyz.resize(num_itf);
  fixed_node_id.resize(num_itf);

  rotated_layer_ien.resize(num_itf);
  rotated_layer_face_id.resize(num_itf);
  init_rotated_node_xyz.resize(num_itf);
  rotated_node_id.resize(num_itf);

  for(int ii=0; ii<num_itf; ++ii)
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( std::to_string(ii) );

    fixed_ele_face_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_cell_face_id" );

    fixed_ele_tag[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_cell_tag" );

    fixed_layer_ien[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_cell_ien" );

    fixed_node_xyz[ii] = h5r -> read_doubleVector( subgroup_name.c_str(), "fixed_node_xyz" );

    fixed_node_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_node_map" );

    rotated_layer_ien[ii].resize(num_tag[ii]);
    rotated_layer_face_id[ii].resize(num_tag[ii]);
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
        rotated_layer_ien[ii][jj] = h5r -> read_intVector( subsubgroup_name.c_str(), "rotated_cell_ien" );

        rotated_layer_face_id[ii][jj] = h5r -> read_intVector( subsubgroup_name.c_str(), "rotated_cell_face_id" );
      }
    }

    init_rotated_node_xyz[ii] = h5r -> read_doubleVector( subgroup_name.c_str(), "rotated_node_xyz" );

    rotated_node_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_node_map" );
  }

  const std::string mesh_info("/Global_Mesh_Info");
  nLocBas = h5r -> read_intScalar(mesh_info.c_str(), "nLocBas");

  delete h5r; H5Fclose( file_id );
}

void ALocal_Interface::print_info() const
{
  SYS_T::commPrint("Interfaces: %d\n", num_itf);
}

Vector_3 ALocal_Interface::get_curr_xyz(const int &ii, const int &node, const double &tt) const
{
  Vector_3 xyz (get_init_rotated_node_xyz(ii, 3 * node),
                get_init_rotated_node_xyz(ii, 3 * node + 1),
                get_init_rotated_node_xyz(ii, 3 * node + 2));

  // rotation around z-axis
  const double angular_velo = MATH_T::PI / 60;  // (rad/s)

  const double rr = std::sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2));

  double angle = MATH_T::get_angle_2d(xyz(1), xyz(2));

  angle += angular_velo * tt;

  xyz(1) = std::cos(angle) * rr;
  xyz(2) = std::sin(angle) * rr;

  return xyz;
}

void ALocal_Interface::get_fixed_ele_ctrlPts(const int &ii, const int &ee,
  double * const volctrl_x,  double * const volctrl_y,  double * const volctrl_z) const
{
  for(int nn{0}; nn < nLocBas; ++nn)
  {
    int node = get_fixed_layer_ien(ii, nLocBas * ee + nn);

    volctrl_x[nn] = get_fixed_node_xyz(ii, 3 * node);
    volctrl_y[nn] = get_fixed_node_xyz(ii, 3 * node + 1);
    volctrl_z[nn] = get_fixed_node_xyz(ii, 3 * node + 2);
  }
}

void ALocal_Interface::get_rotated_ele_ctrlPts(const int &ii, const int &tag, const int &ee, const double &tt,
  double * const volctrl_x, double * const volctrl_y, double * const volctrl_z) const
{
  for(int nn{0}; nn < nLocBas; ++nn)
  {
    int node = get_rotated_layer_ien(ii, tag, nLocBas * ee + nn);
    Vector_3 cuur_node_xyz = get_curr_xyz(ii, node, tt);

    volctrl_x[nn] = cuur_node_xyz(0);
    volctrl_y[nn] = cuur_node_xyz(1);
    volctrl_z[nn] = cuur_node_xyz(2);
  }
}

// EOF