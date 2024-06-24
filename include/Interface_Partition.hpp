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
#include "HDF5_Writer.hpp"
#include "Interface_pair.hpp"

class Interface_Partition
{
  public:
    Interface_Partition( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const std::vector<Interface_pair> &interfaces );

    virtual ~Interface_Partition(){};

    // write the data to hdf5 file in folder /sliding
    virtual void write_hdf5( const std::string &FileName ) const;

  private:
    const int cpu_rank;

    // the number of interface-pairs
    const int num_pair;

    // the number of intervals of each interface-pair
    std::vector<int> num_tag;

    // the number of the local elements of the fixed interfaces
    std::vector<int> num_fixed_part_ele;

    // stores the face id of the volume element of the fixed interface in this part
    std::vector<std::vector<int>> fixed_ele_face_id;

    // stores ien of the volume element of the fixed layer in this part
    std::vector<std::vector<int>> fixed_layer_ien;

    // stores the fixed layer nodes indices, converted by get_old2new()
    std::vector<std::vector<int>> fixed_layer_global_node;

    // stores the fixed layer nodes' coordinates
    std::vector<std::vector<double>> fixed_layer_pt_xyz;

    // stores the interval tag of each element of the fixed interface
    std::vector<std::vector<int>> fixed_interval_tag;

    // stores the face id of all the volume element of the rotated layer
    std::vector<std::vector<std::vector<int>>> rotated_ele_face_id;

    // stores ien of all the volume element of the rotated layer
    std::vector<std::vector<std::vector<int>>> rotated_layer_ien;

    // stores the rotated layer nodes indices, converted by get_old2new()
    std::vector<std::vector<int>> rotated_layer_global_node;

    // stores the rotated layer nodes' coordinates
    std::vector<std::vector<double>> rotated_layer_pt_xyz;

    // // stores the interval tag of each element of the rotated interface
    std::vector<std::vector<int>> rotated_interval_tag;
};

#endif