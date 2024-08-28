#include "Interface_Partition.hpp"

Interface_Partition::Interface_Partition(const IPart * const &part,
  const Map_Node_Index * const &mnindex,
  const std::vector<Interface_pair> &interfaces,
  const std::vector<INodalBC *> &nbc_list) :
  cpu_rank{part->get_cpu_rank()},
  num_pair{VEC_T::get_size(interfaces)}
{
  num_tag.resize(num_pair);

  fixed_nlocalele.resize(num_pair);
  fixed_ele_face_id.resize(num_pair);
  fixed_lien.resize(num_pair);
  fixed_global_node.resize(num_pair);
  fixed_LID.resize(num_pair);
  fixed_pt_xyz.resize(num_pair);
  fixed_interval_tag.resize(num_pair);

  fixed_node_vol_part_tag.resize(num_pair);
  fixed_node_loc_pos.resize(num_pair);

  rotated_ele_face_id.resize(num_pair);
  rotated_lien.resize(num_pair);
  rotated_global_node.resize(num_pair);
  rotated_LID.resize(num_pair);
  rotated_pt_xyz.resize(num_pair);
  rotated_interval_tag.resize(num_pair);

  rotated_node_vol_part_tag.resize(num_pair);
  rotated_node_loc_pos.resize(num_pair);

  const int dof = 4;

  for(int ii=0; ii<num_pair; ++ii)
  {
    fixed_ele_face_id[ii] = std::vector<int> {};

    const std::vector<int> total_fixed_ien = interfaces[ii].get_fixed_vien();
    const std::vector<int> total_fixed_interval_tag = interfaces[ii].get_fixed_interval_tag();

    fixed_lien[ii] = std::vector<int> {};    
    fixed_global_node[ii] = interfaces[ii].get_fixed_global_node();
    fixed_interval_tag[ii] = std::vector<int> {};

    const int num_fixed_node = VEC_T::get_size(fixed_global_node[ii]);
    fixed_node_vol_part_tag[ii].resize(num_fixed_node);
    fixed_node_loc_pos[ii].resize(num_fixed_node);
    fixed_LID[ii].resize(dof * num_fixed_node);

    // convert the GlobalNodeID to new(mapped) global node id
    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_fixed_node; ++jj)
    {
      const int old_gid = fixed_global_node[ii][jj];
      const int new_gid = mnindex->get_old2new(old_gid);
      fixed_global_node[ii][jj] = new_gid;

      for(int kk=0; kk<dof; ++kk)
        fixed_LID[ii][kk * num_fixed_node + jj] = mnindex->get_old2new(nbc_list[kk]->get_ID(old_gid));

      if(part->isNodeInPart( new_gid ))
      {
        fixed_node_vol_part_tag[ii][jj] = cpu_rank;
        fixed_node_loc_pos[ii][jj] = part->get_nodeLocGhoIndex( new_gid );
      }
      else
      {
        fixed_node_vol_part_tag[ii][jj] = -1;
        fixed_node_loc_pos[ii][jj] = -1;
      }
    }

    fixed_pt_xyz[ii] = interfaces[ii].get_fixed_pt_xyz();

    // partition the fixed element according to the cpu_rank
    for(int ee=0; ee<interfaces[ii].get_num_fixed_ele(); ++ee)
    {
      if(interfaces[ii].get_fixed_cpu_rank(ee)==cpu_rank)
      {
        fixed_ele_face_id[ii].push_back(interfaces[ii].get_fixed_faceID(ee));

        for(int jj=0; jj<part->get_nLocBas(); ++jj)
          fixed_lien[ii].push_back(total_fixed_ien[ee * part->get_nLocBas() + jj]);

        fixed_interval_tag[ii].push_back(total_fixed_interval_tag[ee]);
      }
    }

    fixed_nlocalele[ii] = VEC_T::get_size(fixed_ele_face_id[ii]);

    rotated_global_node[ii] = interfaces[ii].get_rotated_global_node();

    const int num_rotated_node = VEC_T::get_size(rotated_global_node[ii]);
    rotated_node_vol_part_tag[ii].resize(num_rotated_node);
    rotated_node_loc_pos[ii].resize(num_rotated_node);
    rotated_LID[ii].resize(dof * num_rotated_node);

    // convert the GlobalNodeID to new(mapped) global node id
    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_rotated_node; ++jj)
    {
      const int old_gid = rotated_global_node[ii][jj];
      const int new_gid = mnindex->get_old2new(old_gid);
      rotated_global_node[ii][jj] = new_gid;

      for(int kk=0; kk<dof; ++kk)
        rotated_LID[ii][kk * num_rotated_node + jj] = mnindex->get_old2new(nbc_list[kk]->get_ID(old_gid));

      if(part->isNodeInPart( new_gid ))
      {
        rotated_node_vol_part_tag[ii][jj] = cpu_rank;
        rotated_node_loc_pos[ii][jj] = part->get_nodeLocGhoIndex( new_gid );
      }
      else
      {
        rotated_node_vol_part_tag[ii][jj] = -1;
        rotated_node_loc_pos[ii][jj] = -1;
      }
    }

    rotated_pt_xyz[ii] = interfaces[ii].get_rotated_pt_xyz();

    rotated_interval_tag[ii] = interfaces[ii].get_rotated_interval_tag();

    std::vector<int> tag = rotated_interval_tag[ii];
    VEC_T::sort_unique_resize(tag);

    // group the rotated element according to the interval tag
    int n_tag = VEC_T::get_size(tag);
    rotated_lien[ii].resize(n_tag);
    rotated_ele_face_id[ii].resize(n_tag);

    std::vector<int> total_rotated_lien = interfaces[ii].get_rotated_vien();
    std::vector<int> total_rotated_ele_face_id = interfaces[ii].get_rotated_faceID();

    for(int jj=0; jj<n_tag; ++jj)
    {
      rotated_lien[ii][jj] = std::vector<int> {};
      rotated_ele_face_id[ii][jj] = std::vector<int> {};
    }

    for(int ee=0; ee<VEC_T::get_size(rotated_interval_tag[ii]); ++ee)
    {
      int ee_tag = rotated_interval_tag[ii][ee];
      for(int jj=0; jj<part->get_nLocBas(); ++jj)
        rotated_lien[ii][ee_tag].push_back(total_rotated_lien[ee * part->get_nLocBas() + jj]);

      rotated_ele_face_id[ii][ee_tag].push_back(total_rotated_ele_face_id[ee]);
    }

    num_tag[ii] = n_tag;
  }
}

void Interface_Partition::write_hdf5(const std::string &FileName) const
{
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, "/sliding", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_intScalar( g_id, "num_interface", num_pair );

  h5w -> write_intVector( g_id, "num_part_fixed_cell", fixed_nlocalele );

  h5w -> write_intVector( g_id, "num_tag", num_tag );

  const std::string groupbase("interfaceid_");

  for(int ii=0; ii<num_pair; ++ii)
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( std::to_string(ii) );

    hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    h5w -> write_intScalar( group_id, "num_fixed_node", VEC_T::get_size(fixed_global_node[ii]) );

    h5w -> write_intVector( group_id, "fixed_cell_face_id", fixed_ele_face_id[ii] );

    h5w -> write_intVector( group_id, "fixed_cell_ien", fixed_lien[ii] );

    h5w -> write_intVector( group_id, "fixed_cell_tag", fixed_interval_tag[ii] );

    h5w -> write_intVector( group_id, "fixed_node_map", fixed_global_node[ii] );

    h5w -> write_intVector( group_id, "fixed_LID", fixed_LID[ii] );

    h5w -> write_doubleVector( group_id, "fixed_pt_xyz", fixed_pt_xyz[ii] );

    h5w -> write_intScalar( group_id, "num_rotated_node", VEC_T::get_size(rotated_global_node[ii]) );

    h5w -> write_intVector( group_id, "rotated_node_map", rotated_global_node[ii] );

    h5w -> write_intVector( group_id, "rotated_LID", rotated_LID[ii] );

    h5w -> write_doubleVector( group_id, "rotated_pt_xyz", rotated_pt_xyz[ii] );

    const std::string subgroupbase("tag_");

    for(int jj=0; jj<num_tag[ii]; ++jj)
    {
      std::string subsubgroup_name(subgroupbase);
      subsubgroup_name.append( std::to_string(jj) );

      hid_t subgroup_id = H5Gcreate(group_id, subsubgroup_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      h5w -> write_intScalar( subgroup_id, "num_rotated_cell", VEC_T::get_size(rotated_ele_face_id[ii][jj]) );

      h5w -> write_intVector( subgroup_id, "rotated_cell_ien", rotated_lien[ii][jj] );

      h5w -> write_intVector( subgroup_id, "rotated_cell_face_id", rotated_ele_face_id[ii][jj] );

      H5Gclose( subgroup_id );
    }

    H5Gclose( group_id );
  }
  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}

// EOF
