#ifndef PART_FEM_ROTATED_HPP
#define PART_FEM_ROTATED_HPP
// ============================================================================
// Part_FEM_Rotated.hpp
//
// Object: Partition 3D tet mesh into subdomains. 
//         The mesh is partitioned into sub-domains that are tagged with a 
//         element tag. In the partitioned file, there will be a element tag,
//         which is an integer array with the length equaling to the number of 
//         local elements.
//
// Date Created: August 8 2024
// ============================================================================
#include "Part_FEM.hpp"

class Part_FEM_Rotated : public Part_FEM
{
  public:
    Part_FEM_Rotated( const IMesh * const &mesh,
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
        const Field_Property &in_fp );

    virtual ~Part_FEM_Rotated() = default;

    virtual void write( const std::string &inputFileName ) const;

  protected:
    std::vector<int> elem_tag {};
    
    // The location in the node_loc for the local fixed/rotated nodes
    // e.g. node_loc_fixed.size() gives the number of nodes belonging
    // to a fixed domain in this partition, and node_loc_fixed[ii] gives
    // the location in node_loc for that node's index.
    // Storing the index in node_loc array makes things more flexible,
    // because sometimes, we manage arrays using local partition's numbering
    // sometimes, we directly set entries using the global numbering.
    std::vector<int> node_loc_fixed {};
    std::vector<int> node_loc_rotated {}; 

    int nlocalnode_fixed, nlocalnode_rotated;
};

#endif
