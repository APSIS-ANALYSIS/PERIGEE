#ifndef INTERFACE_PARTITION_HPP
#define INTERFACE_PARTITION_HPP
// ============================================================================
// Interface_Partition.hpp
// 
// The partition for several interface-pairs.
// 
// Author: Xuanming Huag
// Date: Jun 24 2024
// ============================================================================
#include "IPart.hpp"
#include "Map_Node_Index.hpp"
#include "Interface_pair.hpp"
#include "INodalBC.hpp"

class Interface_Partition
{
  public:
    Interface_Partition( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const std::vector<Interface_pair> &interfaces,
        const std::vector<INodalBC *> &nbc_list);

    virtual ~Interface_Partition(){};

    // write the data to hdf5 file in folder /sliding
    virtual void write_hdf5( const std::string &FileName ) const;

    virtual std::vector<std::vector<int>> get_fixed_node_vol_part_tag() const
    {return fixed_node_vol_part_tag;}

    virtual std::vector<std::vector<int>> get_fixed_node_loc_pos() const
    {return fixed_node_loc_pos;}

    virtual std::vector<std::vector<int>> get_rotated_node_vol_part_tag() const
    {return rotated_node_vol_part_tag;}

    virtual std::vector<std::vector<int>> get_rotated_node_loc_pos() const
    {return rotated_node_loc_pos;}

    virtual int get_fixed_nlocalele(const int &ii) const
    {return fixed_nlocalele[ii];}

    virtual int get_rotated_nlocalele(const int &ii) const
    {return rotated_nlocalele[ii];}

  private:
    const int cpu_rank;

    // the number of interface-pairs
    const int num_pair;

    // the number of intervals of each interface-pair
    std::vector<int> num_tag;

    // the number of the local elements of the fixed interfaces
    std::vector<int> fixed_nlocalele;

    // stores the local fixed element index in this partition
    std::vector<std::vector<int>> fixed_ele_in_this_part;

    // stores the face id of the volume element of the fixed interface
    std::vector<std::vector<int>> fixed_ele_face_id;

    // stores ien of the volume element of the fixed layer
    std::vector<std::vector<int>> fixed_lien;

    // stores the fixed layer nodes indices, converted by get_old2new()
    std::vector<std::vector<int>> fixed_global_node;

    // stores the fixed layer nodes ID array
    std::vector<std::vector<int>> fixed_LID;

    // stores the fixed layer nodes' coordinates
    std::vector<std::vector<double>> fixed_pt_xyz;

    // stores the fixed layer nodes' volume partition tag & local position
    std::vector<std::vector<int>> fixed_node_vol_part_tag;

    std::vector<std::vector<int>> fixed_node_loc_pos;

    // stores the interval tag of each element of the fixed interface
    std::vector<std::vector<int>> fixed_interval_tag;

    // classified by interval tag
    std::vector<std::vector<std::vector<int>>> tagged_fixed_ele;

    std::vector<std::vector<int>> num_tagged_fixed_ele;

    // the number of the local elements of the rotated interfaces
    std::vector<int> rotated_nlocalele;

    // stores the local rotated element index in this partition
    std::vector<std::vector<int>> rotated_ele_in_this_part;

    // stores the face id of all the volume element of the rotated layer
    std::vector<std::vector<int>> rotated_ele_face_id;

    // stores ien of all the volume element of the rotated layer
    std::vector<std::vector<int>> rotated_lien;

    // stores the rotated layer nodes indices, converted by get_old2new()
    std::vector<std::vector<int>> rotated_global_node;

    // stores the rotated layer nodes ID array
    std::vector<std::vector<int>> rotated_LID;

    // stores the rotated layer nodes' coordinates
    std::vector<std::vector<double>> rotated_pt_xyz;

    // stores the rotated layer nodes' volume partition tag & local position
    std::vector<std::vector<int>> rotated_node_vol_part_tag;

    std::vector<std::vector<int>> rotated_node_loc_pos;

    // stores the interval tag of each element of the rotated interface
    std::vector<std::vector<int>> rotated_interval_tag;

    // classified by interval tag
    std::vector<std::vector<std::vector<int>>> tagged_rotated_ele;

    std::vector<std::vector<int>> num_tagged_rotated_ele;
};

#endif
