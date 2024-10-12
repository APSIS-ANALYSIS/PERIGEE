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
  fixed_ele_in_this_part.resize(num_pair);
  fixed_ele_face_id.resize(num_pair);
  fixed_lien.resize(num_pair);
  fixed_global_node.resize(num_pair);
  fixed_LID.resize(num_pair);
  fixed_pt_xyz.resize(num_pair);
  fixed_interval_tag.resize(num_pair);
  tagged_fixed_ele.resize(num_pair);
  num_tagged_fixed_ele.resize(num_pair);

  fixed_node_vol_part_tag.resize(num_pair);
  fixed_node_loc_pos.resize(num_pair);
  
  rotated_nlocalele.resize(num_pair);
  rotated_ele_in_this_part.resize(num_pair);
  rotated_ele_face_id.resize(num_pair);
  rotated_lien.resize(num_pair);
  rotated_global_node.resize(num_pair);
  rotated_LID.resize(num_pair);
  rotated_pt_xyz.resize(num_pair);
  rotated_interval_tag.resize(num_pair);
  tagged_rotated_ele.resize(num_pair);
  num_tagged_rotated_ele.resize(num_pair);

  rotated_node_vol_part_tag.resize(num_pair);
  rotated_node_loc_pos.resize(num_pair);

  L2_proj_fixed_node_pos.resize(num_pair);
  all_fixed_global_node = std::vector<int> {};
  all_fixed_inner_node = std::vector<int> {};

  L2_proj_rotated_node_pos.resize(num_pair);
  all_rotated_global_node = std::vector<int> {};
  all_rotated_inner_node = std::vector<int> {};

  const int dof = 4;

  for(int ii=0; ii<num_pair; ++ii)
  {
    // read info from interfaces
    fixed_ele_face_id[ii] = interfaces[ii].get_fixed_faceID();
    fixed_lien[ii] = interfaces[ii].get_fixed_vien();
    fixed_interval_tag[ii] = interfaces[ii].get_fixed_interval_tag();
    fixed_ele_in_this_part[ii] = std::vector<int> {};
    fixed_global_node[ii] = interfaces[ii].get_fixed_global_node();
    fixed_pt_xyz[ii] = interfaces[ii].get_fixed_pt_xyz();

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
    
    // for L2 projection
    VEC_T::insert_end(all_fixed_global_node, fixed_global_node[ii]);
    L2_proj_fixed_node_pos[ii].resize(num_fixed_node);

    std::vector<int> fixed_inner_node = interfaces[ii].get_fixed_inner_node();
    int num_inner_node = VEC_T::get_size(fixed_inner_node);
    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_inner_node; ++jj)
    {
      const int new_gid = mnindex->get_old2new(fixed_inner_node[jj]);
      fixed_inner_node[jj] = new_gid;
    }

    VEC_T::insert_end(all_fixed_inner_node, fixed_inner_node);

    // partition the fixed element according to the cpu_rank
    for(int ee=0; ee<interfaces[ii].get_num_fixed_ele(); ++ee)
    {
      if(interfaces[ii].get_fixed_cpu_rank(ee)==cpu_rank)
        fixed_ele_in_this_part[ii].push_back(ee);
    }

    fixed_nlocalele[ii] = VEC_T::get_size(fixed_ele_in_this_part[ii]);
    
    // read info from interfaces
    rotated_ele_face_id[ii] = interfaces[ii].get_rotated_faceID();
    rotated_lien[ii] = interfaces[ii].get_rotated_vien();
    rotated_interval_tag[ii] = interfaces[ii].get_fixed_interval_tag();
    rotated_ele_in_this_part[ii] = std::vector<int> {};
    rotated_global_node[ii] = interfaces[ii].get_rotated_global_node();
    rotated_pt_xyz[ii] = interfaces[ii].get_rotated_pt_xyz();

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

    // for L2 projection
    VEC_T::insert_end(all_rotated_global_node, rotated_global_node[ii]);
    L2_proj_rotated_node_pos[ii].resize(num_rotated_node);

    std::vector<int> rotated_inner_node = interfaces[ii].get_rotated_inner_node();
    num_inner_node = VEC_T::get_size(rotated_inner_node);
    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_inner_node; ++jj)
    {
      const int new_gid = mnindex->get_old2new(rotated_inner_node[jj]);
      rotated_inner_node[jj] = new_gid;
    }

    VEC_T::insert_end(all_rotated_inner_node, rotated_inner_node);

    // partition the rotated element according to the cpu_rank
    for(int ee=0; ee<interfaces[ii].get_num_rotated_ele(); ++ee)
    {
      if(interfaces[ii].get_rotated_cpu_rank(ee)==cpu_rank)
        rotated_ele_in_this_part[ii].push_back(ee);
    }

    rotated_nlocalele[ii] = VEC_T::get_size(rotated_ele_in_this_part[ii]);

    // sort local element indices by interval tag
    std::vector<int> tag = rotated_interval_tag[ii];
    VEC_T::sort_unique_resize(tag);
    int n_tag = VEC_T::get_size(tag);

    tagged_fixed_ele[ii].resize(n_tag);
    tagged_rotated_ele[ii].resize(n_tag);

    num_tagged_fixed_ele[ii].resize(n_tag);
    num_tagged_rotated_ele[ii].resize(n_tag);

    num_tag[ii] = n_tag;

    for(int tt=0; tt<n_tag; ++tt)
    {
      tagged_fixed_ele[ii][tt] = std::vector<int> {};
      for(int ee=0; ee<interfaces[ii].get_num_fixed_ele(); ++ee)
      {
        if(fixed_interval_tag[ii][ee] == tt)
          tagged_fixed_ele[ii][tt].push_back(ee);
      }
      num_tagged_fixed_ele[ii][tt] = VEC_T::get_size(tagged_fixed_ele[ii][tt]);

      tagged_rotated_ele[ii][tt] = std::vector<int> {};
      for(int ee=0; ee<interfaces[ii].get_num_rotated_ele(); ++ee)
      {
        if(rotated_interval_tag[ii][ee] == tt)
          tagged_rotated_ele[ii][tt].push_back(ee);
      }
      num_tagged_rotated_ele[ii][tt] = VEC_T::get_size(tagged_rotated_ele[ii][tt]);
    }
  }

  // for L2 projection
  VEC_T::sort_unique_resize(all_fixed_inner_node);
  VEC_T::sort_unique_resize(all_fixed_global_node);
  for(int ii=0; ii<num_pair; ++ii)
  {
    const int num_fixed_node = VEC_T::get_size(fixed_global_node[ii]);

    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_fixed_node; ++jj)
    {
      const int is_inner = VEC_T::get_pos(all_fixed_inner_node, fixed_global_node[ii][jj]);
      if(is_inner == -1)
        L2_proj_fixed_node_pos[ii][jj] = VEC_T::get_pos(all_fixed_global_node, fixed_global_node[ii][jj]);
      else
        L2_proj_fixed_node_pos[ii][jj] = -1;  // Inner node --> Dirichlet node
    }
  }

  for(int nn=0; nn<VEC_T::get_size(all_fixed_inner_node); ++nn)
  {
    const int inner_index = VEC_T::get_pos(all_fixed_global_node, all_fixed_inner_node[nn]);
    all_fixed_inner_node[nn] = inner_index;
  }

  num_all_fixed_node = VEC_T::get_size(all_fixed_global_node);

  VEC_T::sort_unique_resize(all_rotated_inner_node);
  VEC_T::sort_unique_resize(all_rotated_global_node);
  for(int ii=0; ii<num_pair; ++ii)
  {
    const int num_rotated_node = VEC_T::get_size(rotated_global_node[ii]);

    PERIGEE_OMP_PARALLEL_FOR
    for(int jj=0; jj<num_rotated_node; ++jj)
    {
      const int is_inner = VEC_T::get_pos(all_rotated_inner_node, rotated_global_node[ii][jj]);
      if(is_inner == -1)
        L2_proj_rotated_node_pos[ii][jj] = VEC_T::get_pos(all_rotated_global_node, rotated_global_node[ii][jj]);
      else
        L2_proj_rotated_node_pos[ii][jj] = -1;  // Inner node --> Dirichlet node
    }
  }

  for(int nn=0; nn<VEC_T::get_size(all_rotated_inner_node); ++nn)
  {
    const int inner_index = VEC_T::get_pos(all_rotated_global_node, all_rotated_inner_node[nn]);
    all_rotated_inner_node[nn] = inner_index;
  }

  num_all_rotated_node = VEC_T::get_size(all_rotated_global_node);
}

void Interface_Partition::write_hdf5(const std::string &FileName) const
{
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, "/sliding", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_intScalar( g_id, "num_interface", num_pair );

  h5w -> write_intVector( g_id, "num_tag", num_tag );

  h5w -> write_intVector( g_id, "num_local_fixed_cell", fixed_nlocalele );

  h5w -> write_intScalar( g_id, "num_all_fixed_node", num_all_fixed_node );

  h5w -> write_intVector( g_id, "all_fixed_inner_node", all_fixed_inner_node );

  h5w -> write_intVector( g_id, "num_local_rotated_cell", rotated_nlocalele );

  h5w -> write_intScalar( g_id, "num_all_rotated_node", num_all_rotated_node );

  h5w -> write_intVector( g_id, "all_rotated_inner_node", all_rotated_inner_node );

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

    h5w -> write_intVector( group_id, "L2_proj_fixed_node_pos", L2_proj_fixed_node_pos[ii] );

    h5w -> write_intVector( group_id, "fixed_LID", fixed_LID[ii] );

    h5w -> write_doubleVector( group_id, "fixed_pt_xyz", fixed_pt_xyz[ii] );

    h5w -> write_intVector(group_id, "local_fixed_cell", fixed_ele_in_this_part[ii] );

    h5w -> write_intVector(group_id, "num_tagged_fixed_cell", num_tagged_fixed_ele[ii]);

    h5w -> write_intScalar( group_id, "num_rotated_node", VEC_T::get_size(rotated_global_node[ii]) );

    h5w -> write_intVector( group_id, "rotated_cell_face_id", rotated_ele_face_id[ii] );

    h5w -> write_intVector( group_id, "rotated_cell_ien", rotated_lien[ii] );

    h5w -> write_intVector( group_id, "rotated_cell_tag", rotated_interval_tag[ii] );

    h5w -> write_intVector( group_id, "rotated_node_map", rotated_global_node[ii] );

    h5w -> write_intVector( group_id, "L2_proj_rotated_node_pos", L2_proj_rotated_node_pos[ii] );

    h5w -> write_intVector( group_id, "rotated_LID", rotated_LID[ii] );

    h5w -> write_doubleVector( group_id, "rotated_pt_xyz", rotated_pt_xyz[ii] );

    h5w -> write_intVector( group_id, "local_rotated_cell", rotated_ele_in_this_part[ii] );

    h5w -> write_intVector( group_id, "num_tagged_rotated_cell", num_tagged_rotated_ele[ii] );

    const std::string subgroupbase("tag_");

    for(int jj=0; jj<num_tag[ii]; ++jj)
    {
      std::string subsubgroup_name(subgroupbase);
      subsubgroup_name.append( std::to_string(jj) );

      hid_t subgroup_id = H5Gcreate(group_id, subsubgroup_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      h5w -> write_intVector( subgroup_id, "tagged_fixed_cell", tagged_fixed_ele[ii][jj] );

      h5w -> write_intVector( subgroup_id, "tagged_rotated_cell", tagged_rotated_ele[ii][jj] );

      H5Gclose( subgroup_id );
    }

    H5Gclose( group_id );
  }
  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}

// EOF
