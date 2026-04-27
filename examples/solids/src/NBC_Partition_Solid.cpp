#include "NBC_Partition_Solid.hpp"
#include "HDF5_Group.hpp"
#include "HDF5_Writer.hpp"

#include <algorithm>
#include <unordered_set>

NBC_Partition_Solid::NBC_Partition_Solid( const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const std::vector<NodalBC_Solid *> &solid_nbc_list_x,
    const std::vector<NodalBC_Solid *> &solid_nbc_list_y,
    const std::vector<NodalBC_Solid *> &solid_nbc_list_z )
: cpu_rank( part->get_cpu_rank() )
{
  constexpr int dof = 4;

  std::vector<std::vector<unsigned int>> dir_nodes(dof);
  std::vector<std::vector<unsigned int>> disp_nodes(dof);
  std::vector<std::unordered_set<unsigned int>> dir_node_sets(dof);

  const std::vector<std::vector<NodalBC_Solid *>> solid_nbc_lists
  {
    {},
    solid_nbc_list_x,
    solid_nbc_list_y,
    solid_nbc_list_z
  };

  for(int field=1; field<dof; ++field)
  {
    for(const auto * const nbc : solid_nbc_lists[field])
    {
      for(unsigned int ii=0; ii<nbc->get_num_dir_nodes(); ++ii)
      {
        const unsigned int node = nbc->get_dir_nodes(ii);
        dir_nodes[field].push_back(node);
        dir_node_sets[field].insert(node);
        if(nbc->get_is_disp_driven()) disp_nodes[field].push_back(node);
      }
    }
  }

  for(int field=1; field<dof; ++field)
  {
    std::sort(dir_nodes[field].begin(), dir_nodes[field].end());
    dir_nodes[field].erase(
        std::unique(dir_nodes[field].begin(), dir_nodes[field].end()),
        dir_nodes[field].end());

    std::sort(disp_nodes[field].begin(), disp_nodes[field].end());
    disp_nodes[field].erase(
        std::unique(disp_nodes[field].begin(), disp_nodes[field].end()),
        disp_nodes[field].end());
  }

  LID.clear(); LDN.clear(); LPSN.clear(); LPMN.clear();
  LocalMaster.clear(); LocalMasterSlave.clear();
  Num_LD.assign(dof, 0);
  Num_LPS.assign(dof, 0);
  Num_LPM.assign(dof, 0);
  LDN_is_disp_driven.clear();

  for(int field=0; field<dof; ++field)
  {
    std::unordered_set<unsigned int> disp_node_set(
        disp_nodes[field].begin(), disp_nodes[field].end());

    for(const auto node_old : dir_nodes[field])
    {
      const unsigned int node_new = mnindex->get_old2new(node_old);

      if(part->isNodeInPart(node_new))
      {
        LDN.push_back(node_new);
        LDN_is_disp_driven.push_back(disp_node_set.count(node_old) > 0 ? 1 : 0);
        Num_LD[field] += 1;
      }
    }
  }

  LDN.shrink_to_fit();
  LDN_is_disp_driven.shrink_to_fit();

  const int totnode = part->get_nlocghonode();
  LID.resize(totnode * dof);

  PERIGEE_OMP_PARALLEL_FOR
  for(int field=0; field<dof; ++field)
  {
    for(int jj=0; jj<totnode; ++jj)
    {
      const int new_index = part->get_local_to_global(jj);
      const int old_index = mnindex->get_new2old(new_index);

      if(dir_node_sets[field].count(static_cast<unsigned int>(old_index)) > 0)
        LID[field * totnode + jj] = -1;
      else
        LID[field * totnode + jj] = new_index;
    }
  }

  LID.shrink_to_fit();
}

void NBC_Partition_Solid::write_hdf5( const std::string &FileName,
    const std::string &GroupName ) const
{
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );
  auto h5writer = SYS_T::make_unique<HDF5_Writer>(fName, H5F_ACC_RDWR);
  const hid_t file_id = h5writer->get_file_id();

  auto nbc_group = HDF5_Group::create( file_id, GroupName );

  h5writer->write_intVector( nbc_group.id(), "LID", LID );
  if( LDN.size() > 0 ) h5writer->write_intVector( nbc_group.id(), "LDN", LDN );
  h5writer->write_intVector( nbc_group.id(), "Num_LD", Num_LD );
  h5writer->write_intVector( nbc_group.id(), "Num_LPS", Num_LPS );
  h5writer->write_intVector( nbc_group.id(), "Num_LPM", Num_LPM );
  h5writer->write_intVector( nbc_group.id(), "LDN_is_disp_driven", LDN_is_disp_driven );
}

// EOF
