#include "ALocal_Interface.hpp"

ALocal_Interface::ALocal_Interface( const std::string &fileBaseName, const int &cpu_rank,
    const double &angular, const Vector_3 &point_xyz, const Vector_3 &angular_direc )
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

  angular_velo = angular;
  
  direction_rotated = angular_direc;

  point_rotated = point_xyz;
}

void ALocal_Interface::print_info() const
{
  SYS_T::commPrint("Interfaces: %d\n", num_itf);
}

// Get the radius of rotation
Vector_3 ALocal_Interface::get_radius (const Vector_3 &coor) const
{ 
  // The vector from the rotation point to the input point
  const Vector_3 point_rotated_to_coor (coor.x() - point_rotated.x(), coor.y() - point_rotated.y(), coor.z() - point_rotated.z());

  const double projectd_length = Vec3::dot_product(point_rotated_to_coor, direction_rotated);
  
  // The projection point of the input point on the rotation axis
  const Vector_3 point_projected (point_rotated.x() +  projectd_length * direction_rotated.x(), point_rotated.y() +  projectd_length * direction_rotated.y(), point_rotated.z() +  projectd_length * direction_rotated.z());
  
  // The vector from the projection point to the input point
  return Vector_3 (coor.x()- point_projected.x(), coor.y()- point_projected.y(), coor.z()- point_projected.z());
}

void ALocal_Interface::get_currPts( const double * const &ept_x,
    const double * const &ept_y,
    const double * const &ept_z,
    const double &tt,
    double * const &currPt_x,
    double * const &currPt_y,
    double * const &currPt_z,
    const int &type) const
{
  double mag_angular_velo = 0.0; // (rad/s)

  for(int ii=0; ii<nLocBas; ++ii)
  {
    const Vector_3 ept_xyz (ept_x[ii], ept_y[ii], ept_z[ii]);
    const Vector_3 radius_ept = get_radius(ept_xyz);

    const double rr = radius_ept.norm2();
    
    double angle = 0.0;

      //case 0: x-axis, case 1: y-axis, case 2: z-axis
    switch(type) 
    {
      case 0:
        mag_angular_velo = angular_velo * direction_rotated.x();
        angle = MATH_T::get_angle_2d(ept_xyz(1), ept_xyz(2));        
        angle += mag_angular_velo * tt;
        currPt_x[ii] = ept_x[ii];
        currPt_y[ii] = std::cos(angle) * rr;
        currPt_z[ii] = std::sin(angle) * rr;            
        break;
      case 1: 
        mag_angular_velo = angular_velo * direction_rotated.y();
        angle = MATH_T::get_angle_2d(ept_xyz(2), ept_xyz(0));        
        angle += mag_angular_velo * tt;
        currPt_x[ii] = std::sin(angle) * rr;
        currPt_y[ii] = ept_y[ii];
        currPt_z[ii] = std::cos(angle) * rr;            
        break;            
      case 2: 
        mag_angular_velo = angular_velo * direction_rotated.z();
        angle = MATH_T::get_angle_2d(ept_xyz(0), ept_xyz(1));        
        angle += mag_angular_velo * tt;
        currPt_x[ii] = std::cos(angle) * rr;
        currPt_y[ii] = std::sin(angle) * rr;
        currPt_z[ii] = ept_z[ii];            
        break;             
      default:
        SYS_T::print_fatal("Error: ALocal_Interface::get_currPts: No such type of rotation axis. \n");
        break;        
    }
  }
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

void ALocal_Interface::get_rotated_ele_ctrlPts(const int &ii, const int &tag, const int &ee, const double &tt,
  double * const volctrl_x, double * const volctrl_y, double * const volctrl_z) const
{
  // std::vector<double> initPt_x(nLocBas, 0.0), initPt_y(nLocBas, 0.0), initPt_z(nLocBas, 0.0);

  for(int nn{0}; nn < nLocBas; ++nn)
  {
    int node = get_rotated_lien(ii, tag, nLocBas * ee + nn);

    // initPt_x[nn] = rotated_pt_xyz[ii][3 * node];
    // initPt_y[nn] = rotated_pt_xyz[ii][3 * node + 1];
    // initPt_z[nn] = rotated_pt_xyz[ii], [3 * node + 2];
    volctrl_x[nn] = rotated_pt_xyz[ii][3 * node];
    volctrl_y[nn] = rotated_pt_xyz[ii][3 * node + 1];
    volctrl_z[nn] = rotated_pt_xyz[ii][3 * node + 2];
  }

  // get_currPts(&initPt_x[0], &initPt_y[0], &initPt_z[0], tt, volctrl_x, volctrl_y, volctrl_z, 0);
}

// EOF
